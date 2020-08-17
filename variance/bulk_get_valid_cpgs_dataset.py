"""
Take all the methylation data created from the bulk pipeline and create a dataset which combine all the
chromosome to one, next we will use this to create the data for a NN.
We provide information about the pmd, the avg meth, the var in meth, the avg meth in normal cells and allow
to filter based on list of cells to use for calculation, number of cells cpg need to be valid and a ways to
remove percentage of top and bottom coverage cpgs
"""

import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd
import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import consts, sequence_tools, files_tools, utils
from format_files import handle_pmds


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--orig_meth_cells', required=True,
                        help='Path to file containing the list of cell to calculate original cells')
    parser.add_argument('--variance_meth_cells', help='Path to file containing the cells to calculate the '
                                                      'average methylation and variance', required=True)
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


def get_pmd_df(chromosome_file):
    """
    Get the data frame with information about the chromosome methylation only for pmd
    :param chromosome_file: The file with the information of the cpg in a specific chromosome
    :return: The data frame with information only for PMD
    :return: The chromosome we are working on and the data frame of this chromosome
    :rtype: tuple(str, pd.DataFrame)
    """
    chromosome = consts.CHR_FULL_NAME_RE.findall(os.path.basename(chromosome_file))[0]
    df = pd.read_pickle(chromosome_file)

    filtered_df = handle_pmds.filtered_out_non_pmd(df, chromosome, pmd_file=consts.PMD_FILE_LOCAL_DROR)
    filtered_df["pmd_index"] = handle_pmds.get_pmd_index(filtered_df, chromosome,
                                                         pmd_file=consts.PMD_FILE_LOCAL_DROR)
    return chromosome, filtered_df


def get_cpgs_orig_methylated(df, cells_to_use):
    """
    Get a list of cpg which where methylated to begin with using a specific set of cells
    :param df: The df to work on
    :type df: pd.DataFrame
    :param cells_to_use: A list of cells to use, should match the rows names of the df
    :return: The indexes of CpG to use
    """
    df_orig = df.filter(items=cells_to_use, axis=0)
    mean_values = df_orig.mean()
    return mean_values


def get_cells_to_calculate_original_meth(file_path):
    """
    Get a list of cells to define the original methylation level of each cpg
    The file format can change but in our file is csv with cell name called "sample" and a column
    "nc_cells" with 0 or 1 if we need to use this cell to calculate the original methylation
    :param file_path: The path for the file
    :return: A list of cells to use to calculate the original methylation.
    """
    cells_info_data = pd.read_csv(file_path)
    nc_cells = list(cells_info_data[cells_info_data["nc_cells"] == 1]["sample"])
    return [i.strip() for i in nc_cells]


def get_cells_to_calculate_variance(file_path, percentage_of_extreme_cells=15):
    """
    Get a list of cells to use to calculate the methylation and variance
    In this case we have a list of all cells that might be possible and the avg methylation and variance of
    them and we want to take the percentage to take top and bottom cells
    :param file_path:
    :return:
    """
    with open(file_path, "rb") as f:
        all_cells = pickle.load(f)

    all_cells.sort(key=lambda x: x[1])
    low = percentage_of_extreme_cells * len(all_cells) / 100 + 1
    top = len(all_cells) - percentage_of_extreme_cells * len(all_cells) / 100 - 1
    cells = all_cells[:int(low)] + all_cells[int(top):]
    return [i[0] for i in cells]


def main():
    args = parse_input()

    all_files = files_tools.get_files_to_work(args.methylation_folder, consts.BULK_FILE_FORMAT % "*")
    orig_meth_cells = get_cells_to_calculate_original_meth(args.orig_meth_cells)
    variance_cells = get_cells_to_calculate_variance(args.variance_meth_cells,
                                                     args.percentage_of_extreme_cells)

    chromosomes_list = []
    for chromosome_file in tqdm.tqdm(all_files):
        chromosome, pmd_df = get_pmd_df(chromosome_file)
        orig_meth_values = get_cpgs_orig_methylated(pmd_df, orig_meth_cells)
        df = pmd_df.filter(items=variance_cells, axis=0)

        # Remove extreme cells and cells with not enough coverage
        df = utils.remove_extreme_cpgs_by_coverage(df, args.coverage_perc_cpgs_to_remove)
        df = df[np.sum(~pd.isnull(df), axis=0) >= args.min_cells_threshold]
        pmd_sample_mean = df.mean(axis=1)  # Mean of each sample in pmd

        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = df.columns
        chromosome_df["chromosome"] = chromosome
        chromosome_df = chromosome_df.set_index("location")
        chromosome_df["pmd_index"] = handle_pmds.get_pmd_index(chromosome_df, chromosome)
        chromosome_df["meth"] = df.mean()
        chromosome_df["var"] = df.var()
        chromosome_df["pearson_corr"] = df.corrwith(pmd_sample_mean)
        covariance = [df.iloc[:, i].cov(pmd_sample_mean) for i in range(df.shape[1])]
        chromosome_df["covariance"] = covariance
        chromosome_df["orig_meth"] = orig_meth_values[chromosome_df.index]
        chromosome_df["sequence"] = sequence_tools.get_sequences_for_cpgs(df.columns, str(chromosome))
        chromosomes_list.append(chromosome_df.reset_index())

    all_chromosome_df = pd.concat(chromosomes_list)
    all_chromosome_df.to_pickle(os.path.join(args.output_folder, "valid_cpg_zhou.pkl"))


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
