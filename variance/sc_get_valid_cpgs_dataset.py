"""
Take all the methylation data created from the scWGBS pipeline and create a dataset which combine all the
patients and chromosome to one, next we will use this to create the data for a NN
We provide information about the pmd, the avg meth, the var in meth, the avg meth in normal cells for this
patient and across all patient and allow to filter based on another set of boundaries (not only pmd),
to provide list of cells to use and minimum number of samples per group of cells to be a valid cpg and a
ways to remove percentage of top and bottom coverage cpgs
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import files_tools, consts, utils, sequence_tools
from format_files import handle_pmds

PATIENTS = ["CRC01", "CRC11", "CRC13", "CRC10"]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take',
                        required=False, default=None)
    parser.add_argument('--nc_files', help='Path to nc files', required=True)
    parser.add_argument('--cells_to_use', help='Path to cells to use file, default is all', required=False,
                        type=str, default=None)
    parser.add_argument('--min_cells_threshold', help='minimum samples per cell group to be a valid '
                                                      'cpg', required=False, default=5, type=int)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--coverage_perc_cpgs_to_remove', type=int, default=5, required=False,
                        help='Percentage of cells to remove based on extreme coverage(high and low)')
    args = parser.parse_args()
    return args


def get_chromosome_df(file_path, cells_to_use=None):
    patient, chromosome = consts.DATA_FILE_SCWGBS_RE.findall(file_path)[0]
    if patient not in PATIENTS:
        return None, None, None

    if cells_to_use is not None and patient not in cells_to_use:
        return None, None, None

    return patient, chromosome, pd.read_pickle(file_path)


def filter_chromosome_df(df, patient, chromosome, boundaries_data=None, min_cells_threshold=5,
                         cells_to_use=None, perc_of_cpg_to_remove_based_on_coverage=5):
    """
    - filter out non-pmd cpg (and maybe add index)
    - Filtered out again based on boundaries
    - Only take non normal cells
    :param df: The df
    :param patient: name of the patient
    :param chromosome: name of the chromosome
    :param boundaries_data: The boundaries data, can be empty
    :param cells_to_use: A dictionary mapping patient to a dictionary of groups of cells to use
    {p1: {group1:[c1,c2,c3], group2:[c5,c8,c1]}, p2:...}
    :param min_cells_threshold: A min number of cells needed to be a valid cpg in each cells to use
    :param perc_of_cpg_to_remove_based_on_coverage: Percentage of cpg to remove from top and low coverage
    :return: A new dictionary in the same format of the input but filtered by what is required
    """
    # This need to be true
    filtered_df = handle_pmds.filtered_out_non_pmd(df, chromosome, pmd_file=consts.PMD_FILE_LOCAL_DROR)
    filtered_df = utils.filter_df_based_on_region_name(filtered_df, region_name=utils.NOT_NC)

    if boundaries_data:
        filtered_df = utils.filter_df_based_on_tuple_list(filtered_df, boundaries_data)

    if cells_to_use:
        for cell_group in cells_to_use[patient]:
            cells = cells_to_use[patient][cell_group]

            group_coverage = np.sum(~pd.isnull(filtered_df.loc[cells]), axis=0)
            filtered_df = filtered_df[group_coverage >= min_cells_threshold]

    else:
        group_coverage = np.sum(~pd.isnull(filtered_df), axis=0)
        filtered_df = filtered_df.loc[:, group_coverage >= min_cells_threshold]

    # filtered out low and high coverage cpg
    filtered_df = utils.remove_extreme_cpgs_by_coverage(
        filtered_df, top_low_level_to_remove=perc_of_cpg_to_remove_based_on_coverage)

    return filtered_df


def main():
    args = parse_input()
    cells_to_use = files_tools.load_compressed_pickle(args.cells_to_use) if args.cells_to_use else None
    min_cells_threshold = args.min_cells_threshold

    all_files = files_tools.get_files_to_work(os.path.join(args.methylation_folder, "*"), pattern="*.pkl.zip")
    boundaries_data = files_tools.load_compressed_pickle(args.windows_file)

    df_list = []
    for file_path in tqdm.tqdm(all_files):
        patient, chromosome, chromosome_df = get_chromosome_df(file_path, cells_to_use=cells_to_use)
        if not patient:
            continue

        # TODO: boundaries data should be with chr1 and not 1
        filtered_df = filter_chromosome_df(chromosome_df, patient=patient, chromosome=chromosome,
                                           boundaries_data=boundaries_data[int(chromosome[3:])],
                                           cells_to_use=cells_to_use,
                                           min_cells_threshold=min_cells_threshold,
                                           perc_of_cpg_to_remove_based_on_coverage=args.coverage_perc_cpgs_to_remove)

        # Create df with location, chromosome and sequence
        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = filtered_df.columns
        chromosome_df["chromosome"] = chromosome
        chromosome_df["patient"] = patient
        chromosome_df = chromosome_df.set_index("location")

        chromosome_df["sequence"] = sequence_tools.get_sequences_for_cpgs(chromosome_df.index, chromosome)
        chromosome_df["pmd_index"] = handle_pmds.get_pmd_index(chromosome_df, chromosome,
                                                               pmd_file=consts.PMD_FILE_LOCAL_DROR)

        original_meth = utils.cpg_meth_in_cells(patient, chromosome, chromosome_df.index,
                                                args.methylation_folder, sublineage_name=utils.ONLY_NC)
        chromosome_df["orig_meth"] = original_meth.mean()

        pmd_sample_mean = filtered_df.mean(axis=1)

        nc_meth_avg = utils.get_nc_avg(chromosome, chromosome_df.index, args.nc_files)
        chromosome_df.loc[nc_meth_avg.index, "orig_meth_avg"] = nc_meth_avg

        chromosome_df["meth"] = filtered_df.mean()
        chromosome_df["var"] = filtered_df.var()
        chromosome_df["pearson_corr"] = filtered_df.corrwith(pmd_sample_mean)
        covariance = [filtered_df.iloc[:, i].cov(pmd_sample_mean) for i in range(filtered_df.shape[1])]
        chromosome_df["coveriance"] = covariance

        df_list.append(chromosome_df.reset_index())

    all_chromosome_df = pd.concat(df_list).reset_index()
    all_chromosome_df.to_pickle(os.path.join(args.output_folder, "scwgbs_valid_cpg_df.pkl"))

if __name__ == '__main__':
    main()
