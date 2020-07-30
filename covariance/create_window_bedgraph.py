"""
Read a bedgraph file and with information per positions and a window based file and create 2 files:
1. bedgraph file with the values normalized based on the values we have in this bedgraph
2. bedgraph file with the median of each window for all positions in this window
"""

import argparse
import os
import sys

import numpy as np
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools, consts

BEDGRPH_FILE_FORMAT = "*.bedgraph"
SMOOTH_OUTPUT = "smooth_%s_%s_window_%s.bedgraph"
NORMED_OUTPUT = "norm_smooth_%s_%s_window_%s.bedgraph"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Path to bedgraph files or folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--window_size', help='The window size, default is 5000', default=5000,
                        required=False, type=int)
    parser.add_argument('--window_boundaries', help='File with the window boundries', required=False)

    args = parser.parse_args()
    return args


def create_window_bedgraph(bedgraph_path, output_folder, window_size, window_boundaries_path):
    """
    Read a bedgraph file and with information per positions and a window based file and create 2 files:
    1. bedgraph file with the values normalized based on the values we have in this bedgraph
    2. bedgraph file with the median of each window for all positions in this window
    :param bedgraph_path: The path for the bedgraph file
    :param output_folder: The output folder path
    :param window_size: The window size
    :param window_boundaries_path: The boundaries of the windows
    """
    patient, chromosome = consts.PATIENT_CHR_NAME_RE.findall(bedgraph_path)[0]

    output_path = os.path.join(output_folder, "smooth", SMOOTH_OUTPUT % (patient, chromosome, window_size))
    output_path_norm = os.path.join(output_folder, "norm", NORMED_OUTPUT % (patient, chromosome, window_size))

    input_file = files_tools.load_bedgraph(bedgraph_path)
    window_boundaries = files_tools.load_compressed_pickle(window_boundaries_path)

    for tup in window_boundaries[chromosome]:
        indices = np.logical_and(tup[0] <= input_file.start, input_file.start <= tup[1])
        window = input_file[indices]
        window_value = window["coverage"].median()
        input_file.loc[indices, "coverage"] = window_value

    input_file.to_csv(output_path, sep="\t", header=None, index=False)

    input_file["cov"] = np.interp(input_file["coverage"], (input_file["cov"].min(), input_file["cov"].max()),
                                  (0, 1))
    input_file.to_csv(output_path_norm, sep="\t", header=None, index=False)


def main():
    args = parse_input()

    all_file_paths = files_tools.get_files_to_work(args.input, pattern=BEDGRPH_FILE_FORMAT)

    smooth_path = os.path.join(args.output_folder, "smooth")
    norm_path = os.path.join(args.output_folder, "norm")
    if not os.path.exists(smooth_path):
        os.mkdir(smooth_path)
    if not os.path.exists(norm_path):
        os.mkdir(norm_path)

    for file_path in tqdm(all_file_paths):
        create_window_bedgraph(file_path, args.output_folder, args.window_size, args.window_boundaries)


if __name__ == '__main__':
    main()
