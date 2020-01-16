import datetime
import argparse
import glob
import json
import os
import re
import statistics
import sys
from collections import Counter
from tqdm import tqdm

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
from format_files import coverage_between_pairs
from format_files import format_chr_cpg_seq
from commons import tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

WINDOWS_SIZE = 5000


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def mean_of_counter(counter):
    sum_of_numbers = 0
    for i in counter.items():
        sum_of_numbers += i[0] * i[1]

    count = sum(counter.values())
    mean = sum_of_numbers / count
    return mean


def median_of_counter(counter):
    l = []
    for item in counter.items():
        l += [item[1]] * int(item[0])

    return statistics.median(l)


def extend_lists(value_list, each_value_size):
    l = []
    for i in range(len(value_list)):
        l += [value_list[i]] * each_value_size[i]

    return l


def counter_to_dict(counter):
    return {str(item[0]): item[1] for item in counter.items()}


def collect_data(file_path, chromosome):
    window_start_list, window_end_list, window_size = [], [], []
    number_of_samples_list, col_coverage_avg_list, col_coverage_median_list, col_coverage_json_list \
        = [], [], [], []
    window_coverage_avg_list, window_coverage_median_list, window_coverage_json_list = [], [], []
    df = pd.read_pickle(file_path)
    # df = df.iloc[:, :12000]

    locations = df.columns._values
    avg_rate = df.mean(axis=0, skipna=True)
    # median_rate = np.nanmedian(df._values, axis=0)

    converted_matrix = np.where(~np.isnan(df), 1, 0)
    max_size = converted_matrix.shape[1]
    # start_time = datetime.datetime.now()
    # current_index = 0
    for i in range(0, converted_matrix.shape[1], WINDOWS_SIZE):
        # if current_index % 5 == 0:
        #     print("%s / %s, previous took %s" % (
        #         current_index, converted_matrix.shape[1] / WINDOWS_SIZE, datetime.datetime.now() - start_time))
        #     start_time = datetime.datetime.now()
        # current_index += 1

        window = converted_matrix[:, i:i + min(WINDOWS_SIZE, max_size)]
        window_counter = Counter()

        window_start_list.append(locations[i])
        window_end_list.append(locations[i + window.shape[1] - 1])
        window_size.append(window.shape[1])

        for index in range(window.shape[1]):
            cpg_counter = Counter()
            col_coverage = coverage_between_pairs.compare_matrix(window, index)
            number_of_samples, col_coverage = col_coverage[0], col_coverage[1:]

            number_of_samples_list.append(number_of_samples)

            cpg_counter.update(col_coverage)
            window_counter.update(col_coverage)

            col_coverage_avg_list.append(np.average(col_coverage))
            col_coverage_median_list.append(np.median(col_coverage))
            col_coverage_json_list.append(json.dumps(counter_to_dict(cpg_counter)))

        window_coverage_avg_list.append(mean_of_counter(window_counter))
        window_coverage_median_list.append(median_of_counter(window_counter))
        window_coverage_json_list.append(json.dumps(counter_to_dict(window_counter)))

    window_start_list = extend_lists(window_start_list, window_size)
    window_end_list = extend_lists(window_end_list, window_size)
    window_size_list = extend_lists(window_size, window_size)
    window_coverage_avg_list = extend_lists(window_coverage_avg_list, window_size)
    window_coverage_median_list = extend_lists(window_coverage_median_list, window_size)
    window_coverage_json_list = extend_lists(window_coverage_json_list, window_size)

    # Info from dict
    cpg_dict = tools.get_all_cpg_locations_across_chr(full_name=True, full_data=True)
    chr_info = cpg_dict[chromosome]

    context_as_int = chr_info[:, -3]
    context_as_chr = [format_chr_cpg_seq.convert_context_int_to_chr(i) for i in context_as_int]
    is_weak = chr_info[:, -2]
    is_strong = chr_info[:, -1]
    orph_35 = chr_info[:, -9]

    final_table = np.vstack((locations, avg_rate, is_weak, is_strong, context_as_chr, orph_35,
                             col_coverage_avg_list, window_coverage_avg_list,
                             col_coverage_median_list, window_coverage_median_list,
                             window_start_list, window_end_list, window_size_list))

    end_df = pd.DataFrame(final_table.T, index=locations,
                          columns=["locations", "avg_rate", "is_weak", "is_strong", "context", "is_solo(35)",
                                   "cpg_avg_coverage", "window_avg_coverage",
                                   "cpg_median_coverage", "window_median_coverage",
                                   "window_start", "window_end", "window_size"
                                   ]
                          )

    json_table = np.vstack((locations, col_coverage_json_list, window_coverage_json_list,
                            window_start_list, window_end_list, window_size_list))
    json_df = pd.DataFrame(json_table.T, index=locations,
                           columns=["locations", "col_coverage_json_list", "window_coverage_json_list",
                                    "window_start", "window_end", "window_size"]
                           )

    return end_df, json_df


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    for file_path in tqdm(all_cpg_format_file_paths):
        patient, chromosome = coverage_between_pairs.CPG_FORMAT_FILE_RE.findall(file_path)[0]
        output_csv = os.path.join(output, patient, "%s_all_data.csv" % chromosome)
        output_json = os.path.join(output, patient, "%s_json_coverage.pickle.zip" % chromosome)


        data, json_data = collect_data(file_path, chromosome)

        if not os.path.exists(os.path.join(output, patient)):
            os.mkdir(os.path.join(output, patient))

        data.to_csv(output_csv)
        json_data.to_pickle(output_json)


if __name__ == '__main__':
    tools.init_slurm(main)
    # main()
