import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import tqdm

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools, consts, utils
from variance import bulk_get_valid_cpgs_dataset


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--orig_meth_cells', required=True,
                        help='Path to file containing the list of cell to calculate original cells')
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--percentage_of_extreme_cells', required=False, type=int, default=15,
                        help='Percentage of extreme cells to take from top and bottom, default is 15% from '
                             'each side')
    parser.add_argument('--coverage_perc_cpgs_to_remove', type=int, default=5, required=False,
                        help='Percentage of cells to remove based on extreme coverage(high and low)')
    parser.add_argument('--min_cells_threshold', help='minimum samples per cell group to be a valid '
                                                      'cpg', required=False, default=5, type=int)
    args = parser.parse_args()
    return args


def get_cells_to_calculate_variance(file_path):
    cells_info_data = pd.read_csv(file_path)
    nc_cells = list(cells_info_data[cells_info_data["variance_cells"] == 1]["sample"])
    return [i.strip() for i in nc_cells]


def main():
    args = parse_input()

    all_files = files_tools.get_files_to_work(args.methylation_folder, consts.BULK_FILE_FORMAT % "*")
    orig_meth_cells = bulk_get_valid_cpgs_dataset.get_cells_to_calculate_original_meth(args.orig_meth_cells)
    variance_cells = get_cells_to_calculate_variance(args.orig_meth_cells)

    chromosomes_list = []
    for chromosome_file in tqdm.tqdm(all_files):
        chromosome, pmd_df = bulk_get_valid_cpgs_dataset.get_pmd_df(chromosome_file)
        orig_meth_values = bulk_get_valid_cpgs_dataset.get_cpgs_orig_methylated(pmd_df, orig_meth_cells)
        df = pmd_df.filter(items=variance_cells, axis=0).astype(np.float64)
        # df = pmd_df.astype(np.float64)

        # Remove extreme cells and cells with not enough coverage
        df = utils.remove_extreme_cpgs_by_coverage(df, args.coverage_perc_cpgs_to_remove)
        threshold_filter = df.notnull().sum(axis=0) > args.min_cells_threshold
        df = df.filter(items=threshold_filter[threshold_filter].index, axis=1)
        pmd_sample_mean = df.mean(axis=1)  # Mean of each sample in pmd

        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = df.columns
        chromosome_df["chromosome"] = chromosome
        chromosome_df = chromosome_df.set_index("location")
        # chromosome_df["pmd_index"] = df["pmd_index"]
        chromosome_df["meth"] = df.mean()
        chromosome_df["pearson_corr"] = df.corrwith(pmd_sample_mean)
        chromosome_df["orig_meth"] = orig_meth_values[chromosome_df.index]
        # chromosome_df["sequence"] = sequence_tools.get_sequences_for_cpgs(df.columns, str(chromosome))
        chromosomes_list.append(chromosome_df.reset_index())

    all_chromosome_df = pd.concat(chromosomes_list)
    # all_chromosome_df.plot.scatter("meth", "pearson_corr")
    # plt.show()
    all_chromosome_df.to_pickle(os.path.join(args.output_folder, "valid_cpg_zhou_corr.pkl"))


if __name__ == '__main__':
    main()
