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

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import files_tools, consts, utils
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
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--coverage_perc_cpgs_to_remove', type=int, default=5, required=False,
                        help='Percentage of cells to remove based on extreme coverage(high and low)')
    args = parser.parse_args()
    return args


def get_chromosome_df_dict(all_file_paths, cells_to_use=None):
    chromosome_df_dict = {}
    for file_path in all_file_paths:
        patient, chromosome = consts.DATA_FILE_SCWGBS_RE.findall(file_path)[0]
        if patient not in PATIENTS:
            continue

        if cells_to_use is not None and patient not in cells_to_use:
            continue

        elif chromosome not in chromosome_df_dict:
            chromosome_df_dict[chromosome] = []

        df = pd.read_pickle(file_path)

        chromosome_df_dict[chromosome].append((patient, df))

    return chromosome_df_dict


def filter_chromosome_df_dict(chromosome_df_dict, boundaries_data=None, min_cells_threshold=5,
                              cells_to_use=None, perc_of_cpg_to_remove_based_on_coverage=5):
    """
    - For each patient and each file inside it (chromosome file) ->
    - filter out non-pmd cpg (and maybe add index)
    - Filtered out again based on boundaries
    - Only take non normal cells
    :param chromosome_df_dict: The dictionary which matchs chromosome to all patients and df related to it
    :param boundaries_data: The boundaries data, can be empty
    :param cells_to_use: A dictionary mapping patient to a dictionary of groups of cells to use
    {p1: {group1:[c1,c2,c3], group2:[c5,c8,c1]}, p2:...}
    :param min_cells_threshold: A min number of cells needed to be a valid cpg in each cells to use
    :param perc_of_cpg_to_remove_based_on_coverage: Percentage of cpg to remove from top and low coverage
    :return: A new dictionary in the same format of the input but filtered by what is required
    """

    filtered_dict = {chromosome: [] for chromosome in chromosome_df_dict}
    for chromosome in chromosome_df_dict:
        patient, df = chromosome_df_dict[chromosome]
        filtered_df = handle_pmds.filtered_out_non_pmd(df, chromosome, add_pmd_index=True)
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
            filtered_df = filtered_df[group_coverage >= min_cells_threshold]

        # filtered out low and high coverage cpg
        filtered_df = utils.remove_extreme_cpgs_by_coverage(
            filtered_df, top_low_level_to_remove=perc_of_cpg_to_remove_based_on_coverage)
        filtered_dict[chromosome].append(patient, filtered_df)

    return filtered_dict


def main():
    args = parse_input()
    cells_to_use = files_tools.load_compressed_pickle(args.cells_to_use) if args.cells_to_use else None
    min_cells_threshold = args.min_cells_threshold

    all_files = files_tools.get_files_to_work(os.path.join(args.methylation_folder, "*"), pattern="*.pkl.zip")
    boundaries_data = files_tools.load_compressed_pickle(args.windows_file)

    chromosome_df_dict = get_chromosome_df_dict(all_files, cells_to_use=cells_to_use)
    chromosome_df_dict = filter_chromosome_df_dict(
        chromosome_df_dict, boundaries_data, min_cells_threshold=min_cells_threshold,
        cells_to_use=cells_to_use, perc_of_cpg_to_remove_based_on_coverage=args.coverage_perc_cpgs_to_remove)

    chromosomes_list = []
    for chromosome in chromosome_df_dict:
        cpg_locations = utils.get_all_indexes_from_patients_list(chromosome_df_dict[chromosome])

        # Create df with location, chromosome and sequence
        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = cpg_locations
        chromosome_df["chromosome"] = chromosome
        # chromosome_df["sequence"] = sequence_tools.get_sequences_for_cpgs(cpg_locations, chromosome)

        chromosome_df = chromosome_df.set_index("location")

        nc_meth_avg = utils.get_nc_avg(chromosome, chromosome_df.index, args.nc_files)
        chromosome_df.loc[nc_meth_avg.index, "orig_meth_avg"] = nc_meth_avg

        for df_tuple in chromosome_df_dict[chromosome]:
            patient, patient_df = df_tuple
            patient_num = patient[-2:]
            cpgs_indexes = patient_df.columns
            pmd_sample_mean = patient_df.mean(axis=1)

            chromosome_df.loc[cpgs_indexes, "pmd_index"] = patient_df["pmd_index"]

            chromosome_df.loc[cpgs_indexes, "meth%s" % patient_num] = patient_df.mean()
            chromosome_df.loc[cpgs_indexes, "var%s" % patient_num] = patient_df.var()
            chromosome_df["pearson_corr"] = patient_df.corrwith(pmd_sample_mean)
            covariance = [patient_df.iloc[:, i].cov(pmd_sample_mean) for i in range(patient_df.shape[1])]
            chromosome_df["coveriance"] = covariance

            original_meth = utils.cpg_meth_in_cells(patient, chromosome, cpgs_indexes,
                                                    args.methylation_folder,
                                                    sublineage_name=utils.ONLY_NC)
            chromosome_df.loc[original_meth.index, "orig_meth%s" % patient_num] = original_meth.mean()

        chromosomes_list.append(chromosome_df.reset_index())

    all_chromosome_df = pd.concat(chromosomes_list)
    all_chromosome_df.to_pickle(os.path.join(args.output_folder, "valid_cpg.pkl"))


if __name__ == '__main__':
    main()
