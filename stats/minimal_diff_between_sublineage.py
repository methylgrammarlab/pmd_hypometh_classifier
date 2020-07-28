import argparse
import glob
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

HIGH_T = 0.5
LOW_T = 0.25
MID = 0.5

sys.path.append(os.getcwd())
from format_files.format_cpg_context_map import NUMBER_OF_ORPH_PER_INDEX
from commons import consts
from format_files import handle_pmds, format_sublineage_info

ORPH_COLS = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"

BEDGRAPH_FILE_NAME_RE = re.compile(".*_chr_(\d+)_.*")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--sublineage_cov_folder', help='Path to folder or file of parsed scWGBS',
                        required=False)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def format_args():
    """
    Format the args for this script
    :return: The path of the files and the output directory
    """
    args = parse_input()
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])
    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    covariance_dict = get_covariance_dict(args.sublineage_cov_folder)

    return all_cpg_format_file_paths, covariance_dict, output


def get_covariance_dict(covariance_path):
    covariance_path = os.path.join(covariance_path, "*.bedgraph")
    d = {}
    bedgraph_files = glob.glob(covariance_path)
    for f in bedgraph_files:
        name = f.split("and")[1].split(".")[0]
        d[name] = f

    return d


def minimal_diff_in_sub_lineage_meth_levels(df, patient, chromosome):
    sublineage_info = format_sublineage_info.get_sublineage_info(consts.SUBLINEAGE_FILE_LOCAL_DROR)
    patient_info = sublineage_info[patient]
    dfs = []
    for sublineage in patient_info:
        if sublineage == "NC":
            continue

        region_cell_ids = []
        for sample in patient_info[sublineage]:
            region_cell_ids.extend([cell_id for cell_id in df.index if cell_id.startswith(sample)])

        region_df = df.loc[region_cell_ids, :].mean(axis=0, skipna=True)
        region_df[region_df >= HIGH_T] = 1
        region_df[region_df <= LOW_T] = 0
        dfs.append(region_df)

    accuracy = 0
    for couple in list(itertools.combinations(dfs, 2)):
        df1 = couple[0]
        df2 = couple[1]
        index = np.logical_and(np.logical_or(df1 == 0, df1 == 1), np.logical_or(df2 == 0, df2 == 1))
        diff = np.sum(np.abs(df1[index] - df2[index]) == 0)
        accuracy = max(diff / np.sum(index == 1), accuracy)

    print(accuracy)


def minimal_diff_in_sub_lineage_cov(covariance_dict, patient, chromosome):
    ### Unfinished ###

    sublineage_info = format_sublineage_info.get_sublineage_info(consts.SUBLINEAGE_FILE_LOCAL_DROR)
    patient_info = sublineage_info[patient]
    dfs = {}

    for sublineage in patient_info:
        if sublineage == "NC":
            continue

        covariance_pmd_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(covariance_dict[sublineage],
                                                                               chromosome)

    accuracy = 0
    for couple in list(itertools.combinations(dfs, 2)):
        df1 = couple[0]
        df2 = couple[1]
        index = np.logical_and(np.logical_or(df1 == 0, df1 == 1), np.logical_or(df2 == 0, df2 == 1))
        diff = np.sum(np.abs(df1[index] - df2[index]) == 0)
        accuracy = max(diff / np.sum(index == 1), accuracy)

    print(accuracy)


def main():
    input_files, covariance_dict, output_dir = format_args()

    for file_path in tqdm(input_files):
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.filtered_out_non_pmd(df, chromosome)

        # minimal_diff_in_sub_lineage_meth_levels(pmd_df, patient, chromosome)
        minimal_diff_in_sub_lineage_cov(covariance_dict, patient, chromosome)


if __name__ == '__main__':
    main()
