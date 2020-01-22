import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

HIGH_T = 0.75
LOW_T = 0.25

sys.path.append(os.getcwd())
from commons import files_tools
from format_files.format_chr_cpg_seq import NUMBER_OF_ORPH_PER_INDEX

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def collect_data(df, chr_info):
    orph_coloms = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]

    df.reset_index(drop=True)
    sampels = df.axes[0]
    pt_index = [i for i in range(len(sampels)) if sampels[i].startswith("PT")]
    nc_index = [i for i in range(len(sampels)) if sampels[i].startswith("NC")]

    nc_values = df.iloc[nc_index, :].mean(axis=0, skipna=True)
    pt_values = df.iloc[pt_index, :].mean(axis=0, skipna=True)

    pt_values[pt_values >= HIGH_T] = 1
    pt_values[pt_values <= LOW_T] = 0
    nc_values[nc_values <= LOW_T] = 0
    nc_values[nc_values >= HIGH_T] = 1

    ones_or_zeros_pt = np.where(np.logical_or(pt_values == 1, pt_values == 0))[0]
    ones_or_zeros_nc = np.where(np.logical_or(nc_values == 1, nc_values == 0))[0]

    pt_values[~np.logical_or(pt_values == 1, pt_values == 0)] = 2
    nc_values[~np.logical_or(nc_values == 1, nc_values == 0)] = 2
    location = chr_info[:, 0]
    context = chr_info[:, -3]
    orph = chr_info[:, 1:14]

    # Combine the different tables to one table and than convert to df
    final_table = np.hstack(
        (location[:, None], context[:, None], orph, nc_values[:, None], pt_values[:, None]))

    end_df = pd.DataFrame(final_table, columns=["locations", "context"] + orph_coloms + ["nc", "pt"])
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
        cpg_dict = files_tools.get_cpg_context_map(drop_chr_prefix=True, get_full_mapping=True)
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]

        df = pd.read_pickle(file_path)
        data = collect_data(df, cpg_dict[chromosome])

        save_output(data, output_dir, file_path)


def save_output(data, output, data_file_path):
    """
    Save the data
    :param data: The main data of the script
    :param output: The output folder
    :param data_file_path: The original file path
    """
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(data_file_path)[0]
    output_csv = os.path.join(output, patient, "%s_all_data.csv.gzip" % chromosome)

    if not os.path.exists(os.path.join(output, patient)):
        os.mkdir(os.path.join(output, patient))

    data.to_csv(output_csv, compression='gzip')


if __name__ == '__main__':
    main()
