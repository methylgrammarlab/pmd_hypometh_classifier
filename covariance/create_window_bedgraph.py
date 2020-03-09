"""
Convert a covariance bedgraph file which contains information per position to a bedgraph with the same
value per window
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))

CSV_FILE = "common_cpg_in_cov_matrix_%s_chr_%s_region_%s.csv"

BEDGRPH_FILE_FORMAT = "*.bedgraph"
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr_(\d+).*")
OUTPUT_FILE_FORMAT = "smooth_%s_chr%s_window_%s.bedgraph"
OUTPUT_FILE_FORMAT_NORM = "norm_smooth_%s_chr%s_window_%s.bedgraph"


def get_files_to_work(files):
    if os.path.isdir(files):
        file_path = os.path.join(files, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [files]

    return all_file_paths


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Path to bedgrph files or folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--window_size', help='The window size, default is 500', default=5000,
                        required=False, type=int)

    args = parser.parse_args()
    return args


def create_window_bedgraph(file_path, output_folder, window_size):
    patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
    output_path = os.path.join(output_folder, OUTPUT_FILE_FORMAT % (patient, chromosome, window_size))
    output_path_norm = os.path.join(output_folder, OUTPUT_FILE_FORMAT_NORM % (patient, chromosome,
                                                                              window_size))
    input_file = pd.read_csv(file_path, sep="\t", header=None, names=["chr", "start", "end", "cov"])
    number_of_lines = input_file.shape[0]
    for i in range(0, number_of_lines, window_size):
        window = input_file[i:min(number_of_lines, i + window_size - 1)]
        window_value = window["cov"].median()
        input_file.loc[i:min(number_of_lines, i + window_size - 1), "cov"] = window_value

    input_file.to_csv(output_path, sep="\t", header=None, index=False)

    norm_cov = input_file["cov"] - input_file["cov"].min()
    norm_cov /= norm_cov.max()
    input_file["cov"] = norm_cov
    input_file.to_csv(output_path_norm, sep="\t", header=None, index=False)


def main():
    args = parse_input()

    all_file_paths = get_files_to_work(args.input)

    for file_path in all_file_paths:
        create_window_bedgraph(file_path, args.output_folder, args.window_size)


if __name__ == '__main__':
    main()
