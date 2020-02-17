import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

HIGH_T = 0.5
LOW_T = 0.25
MID = 0.5

KEPT_METH = 1
LOST_METH = 0

sys.path.append(os.getcwd())
from format_files.format_cpg_context_map import NUMBER_OF_ORPH_PER_INDEX
from format_files import format_cpg_context_map, handle_pmds

ORPH_COLS = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
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

    return all_cpg_format_file_paths, output


def predict_df(df, chr_info, chromosome, is_filtred_to_pmd=True):
    if not is_filtred_to_pmd:
        df = handle_pmds.get_pmd_df(df, chromosome)

    # weak\strong
    weak = format_cpg_context_map.get_weak_column(chr_info)
    strong = format_cpg_context_map.get_strong_column(chr_info)
    weak_or_strong = np.logical_or(weak, strong)

    # solo\ not solo
    orph35_col = format_cpg_context_map.get_orph_35_column(chr_info)
    orph35_col[orph35_col > 0] = 1
    not_solo = orph35_col.astype(np.bool)
    is_solo = np.logical_not(not_solo)

    nc_index = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    pt_index = [cell_id for cell_id in df.index if cell_id.startswith('PT')]

    nc_values = df.loc[nc_index, :].mean(axis=0, skipna=True)
    pt_values = df.loc[pt_index, :].mean(axis=0, skipna=True)

    pt_values[pt_values >= HIGH_T] = 1
    pt_values[pt_values <= LOW_T] = 0
    nc_values[nc_values <= LOW_T] = 0
    nc_values[nc_values >= HIGH_T] = 1

    valid_indexes = np.logical_and(np.logical_or(pt_values == 1, pt_values == 0), nc_values == 1)
    prediction = []

    for index in range(valid_indexes.size):
        if not valid_indexes.iloc[index]:
            continue

        if not_solo[index]:
            prediction.append(KEPT_METH)

        elif is_solo[index]:
            if weak[index]:
                prediction.append(LOST_METH)

            else:
                prediction.append(KEPT_METH)

    prediction = np.array(prediction)
    diff = prediction - pt_values[valid_indexes]
    accuracy = np.sum(diff == 0) / np.sum(valid_indexes == 1)

    print(accuracy)


def main():
    input_files, output_dir = format_args()
    pmd_context_map = handle_pmds.get_pmd_context_map()

    for file_path in tqdm(input_files):
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.get_pmd_df(df, chromosome)
        predict_df(pmd_df, pmd_context_map[chromosome], chromosome)


if __name__ == '__main__':
    main()
