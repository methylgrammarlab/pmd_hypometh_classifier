import argparse
import glob
import json
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

NUMBER_OF_ORPH_PER_INDEX = [1, 2, 3, 4, 5, 10, 20, 35, 50, 75, 100, 150, 200]

sys.path.append(os.getcwd())
sys.path.append("/cs/usr/drordod/Desktop/project/proj_scwgbs")
from commons import files_tools

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_context_as_int_for_chr(chr_info):
    return chr_info[:, -3]


def collect_data(df, chr_info):
    orph_coloms = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]

    df.reset_index(drop=True)
    sampels = df.axes[0]
    pt_index = [i for i in range(len(sampels)) if sampels[i].startswith("PT")]
    nc_index = [i for i in range(len(sampels)) if sampels[i].startswith("NC")]

    nc_values = df.iloc[nc_index,:].mean(axis=0, skipna=True)
    pt_values = df.iloc[pt_index,:].mean(axis=0, skipna=True)

    ones_or_zeros_pt = np.where(np.logical_or(pt_values == 1, pt_values == 0))[0]
    ones_or_zeros_nc = np.where(np.logical_or(nc_values == 1, nc_values == 0))[0]
    together = set(ones_or_zeros_pt) & set(ones_or_zeros_nc)
    nc_values = nc_values.iloc[list(together)]
    pt_values = pt_values.iloc[list(together)]

    index = pt_values.axes[0]._values
    location = chr_info[:, 0]
    location = np.in1d(location,index)
    context = chr_info[location, -3]
    orph = chr_info[location, 1:14]

    # Combine the different tables to one table and than convert to df
    final_table = np.hstack((index[:, None], context[:, None], orph, nc_values[:,None], pt_values[:,None]))

    end_df = pd.DataFrame(final_table, index=index, columns=["locations", "context"] + orph_coloms + ["nc", "pt"])
    return end_df


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


def main():
    input_files, output_dir = format_args()

    for file_path in tqdm(input_files):
        cpg_dict = files_tools.get_all_cpg_locations_across_chr(full_name=True, full_data=True)
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]

        df = pd.read_pickle(file_path)
        data = collect_data(df, cpg_dict[chromosome])

        save_output(data, output_dir, file_path)


def save_output(data, output, data_file_path):
    """
    Save the data
    :param data: The main data of the script
    :param json_data: The json data
    :param output: The output folder
    :param data_file_path: The original file path
    """
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(data_file_path)[0]
    output_csv = os.path.join(output, patient, "%s_all_data.csv" % chromosome)

    if not os.path.exists(os.path.join(output, patient)):
        os.mkdir(os.path.join(output, patient))
    data.to_csv(output_csv)


if __name__ == '__main__':
    main()
