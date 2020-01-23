"""
Collect the following data: "locations", "avg_rate", "is_weak", "is_strong", "context", "is_solo(35)",
                                   "number_of_samples_list", "cpg_avg_coverage", "window_avg_coverage",
                                   "cpg_median_coverage", "window_median_coverage",
                                   "window_start", "window_end", "window_size"
"""

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

sys.path.append(os.path.dirname(os.getcwd()))
import commons.slurm_tools
from commons.data_tools import mean_of_counter_obj, median_of_counter_obj, extend_lists, numpy_counter_to_dict
from stats import coverage_between_pairs
from format_files import format_cpg_context_map
from commons import files_tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

WINDOWS_SIZE = 5000


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def save_output(data, json_data, output, data_file_path):
    """
    Save the data
    :param data: The main data of the script
    :param json_data: The json data
    :param output: The output folder
    :param data_file_path: The original file path
    """
    patient, chromosome = coverage_between_pairs.CPG_FORMAT_FILE_RE.findall(data_file_path)[0]
    output_csv = os.path.join(output, patient, "%s_all_data.csv" % chromosome)
    output_json = os.path.join(output, patient, "%s_json_coverage.pickle.zip" % chromosome)

    if not os.path.exists(os.path.join(output, patient)):
        os.mkdir(os.path.join(output, patient))
    data.to_csv(output_csv)
    json_data.to_pickle(output_json)


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


def collect_data(df, chromosome):
    """
    Collect all the data we wanted and save it in two tables
    :param df: The data frame
    :param chromosome: The name of the chr
    :return: Two tables with all the data to be saved
    """
    # Collect global sequence information from the cpg dict
    cpg_dict = files_tools.get_cpg_context_map()
    chr_info = cpg_dict[chromosome]

    context_as_chr = format_cpg_context_map.get_context_as_str_for_chr(chr_info)
    is_weak = format_cpg_context_map.get_weak_column(chr_info)
    is_strong = format_cpg_context_map.get_strong_column(chr_info)
    orph_35 = format_cpg_context_map.get_orph_35_column(chr_info)

    # Start to work on generating more data per window

    # Init too many lists to keep the data
    number_of_samples_list, col_coverage_avg_list, col_coverage_median_list, col_coverage_json_list, \
    window_start_list, window_end_list, window_size, window_coverage_avg_list, window_coverage_median_list, \
    window_coverage_json_list = [], [], [], [], [], [], [], [], [], []

    # Get the position of each CpG and calculate the avg rate for every CpG
    # The median is commented out because we only have values between 0-4 and it didn't gave something interesting +
    # it takes very long to run this part
    positions = df.columns._values
    avg_rate = df.mean(axis=0, skipna=True)
    # median_rate = np.nanmedian(df._values, axis=0)

    # Convert all nans to 0 and all values to 1, this make calculating the common coverage easier
    converted_matrix = np.where(~np.isnan(df), 1, 0)
    chr_size = converted_matrix.shape[1]

    for i in range(0, converted_matrix.shape[1], WINDOWS_SIZE):
        window_counter = Counter()
        window = converted_matrix[:, i:i + min(WINDOWS_SIZE, chr_size)]

        window_start_list.append(positions[i])
        window_end_list.append(positions[i + window.shape[1] - 1])
        window_size.append(window.shape[1])

        # Go over each CpG
        for index in range(window.shape[1]):
            cpg_counter = Counter()

            # Calculate the coverage between the CpG and all other pairs
            col_coverage = coverage_between_pairs.compare_matrix(window, index)

            # Extract the number of samples for this CpG, and the coverage with all others CpGs
            # Update the counters for this CpG and this window
            number_of_samples, col_coverage = col_coverage[0], col_coverage[1:]
            number_of_samples_list.append(number_of_samples)
            cpg_counter.update(col_coverage)
            window_counter.update(col_coverage)

            # Calculate and update the lists for average, median for CpG and also save the data as json
            col_coverage_avg_list.append(np.average(col_coverage))
            col_coverage_median_list.append(np.median(col_coverage))
            col_coverage_json_list.append(json.dumps(numpy_counter_to_dict(cpg_counter)))

        # Calculate and update the lists for average, median for window and also save the data as json
        window_coverage_avg_list.append(mean_of_counter_obj(window_counter))
        window_coverage_median_list.append(median_of_counter_obj(window_counter))
        window_coverage_json_list.append(json.dumps(numpy_counter_to_dict(window_counter)))

    # Extend the different lists to fill the correct size
    window_start_list = extend_lists(window_start_list, window_size)
    window_end_list = extend_lists(window_end_list, window_size)
    window_size_list = extend_lists(window_size, window_size)
    window_coverage_avg_list = extend_lists(window_coverage_avg_list, window_size)
    window_coverage_median_list = extend_lists(window_coverage_median_list, window_size)
    window_coverage_json_list = extend_lists(window_coverage_json_list, window_size)

    # Combine the different tables to one table and than convert to df
    final_table = np.vstack((positions, avg_rate, is_weak, is_strong, context_as_chr, orph_35,
                             number_of_samples_list, col_coverage_avg_list, window_coverage_avg_list,
                             col_coverage_median_list, window_coverage_median_list,
                             window_start_list, window_end_list, window_size_list))

    end_df = pd.DataFrame(final_table.T, index=positions,
                          columns=["locations", "avg_rate", "is_weak", "is_strong", "context", "is_solo(35)",
                                   "number_of_samples_list", "cpg_avg_coverage", "window_avg_coverage",
                                   "cpg_median_coverage", "window_median_coverage",
                                   "window_start", "window_end", "window_size"
                                   ]
                          )

    json_table = np.vstack((positions, col_coverage_json_list, window_coverage_json_list,
                            window_start_list, window_end_list, window_size_list))
    json_df = pd.DataFrame(json_table.T, index=positions,
                           columns=["locations", "col_coverage_json_list", "window_coverage_json_list",
                                    "window_start", "window_end", "window_size"]
                           )

    return end_df, json_df


def main():
    input_files, output_dir = format_args()

    for file_path in tqdm(input_files):
        patient, chromosome = coverage_between_pairs.CPG_FORMAT_FILE_RE.findall(file_path)[0]
        output_csv = os.path.join(output_dir, patient, "%s_all_data.csv" % chromosome)

        if os.path.exists(output_csv):
            continue

        df = pd.read_pickle(file_path)
        data, json_data = collect_data(df, chromosome)

        save_output(data, json_data, output_dir, file_path)


if __name__ == '__main__':
    commons.slurm_tools.init_slurm(main)
    # main()
