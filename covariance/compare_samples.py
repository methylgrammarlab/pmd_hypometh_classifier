import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('seaborn')

import numpy as np

HIGH_T = 0.5
LOW_T = 0.25
MID = 0.5

sys.path.append(os.getcwd())
from format_files.format_cpg_context_map import NUMBER_OF_ORPH_PER_INDEX
from format_files import handle_pmds

ORPH_COLS = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"

BEDGRAPH_FILE_NAME_RE = re.compile(".*_chr_(\d+)_.*")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_file', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--second_file', help='Path to folder or file of parsed scWGBS', required=True)
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

    first_file = args.first_file
    second_file = args.second_file

    assert BEDGRAPH_FILE_NAME_RE.findall(first_file)[0] == BEDGRAPH_FILE_NAME_RE.findall(second_file)[0]

    return first_file, second_file, output, BEDGRAPH_FILE_NAME_RE.findall(first_file)[0]


def compare_samples_covariance(first_file, second_file, chromosome):
    first_df = handle_pmds.get_covariance_pmd_df(first_file, chromosome)
    second_df = handle_pmds.get_covariance_pmd_df(second_file, chromosome)

    first_df = first_df[first_df.coverage < 0.1]
    first_df = first_df[first_df.coverage > -0.03]

    second_df = second_df[second_df.coverage < 0.1]
    second_df = second_df[second_df.coverage > -0.03]

    indexes = first_df.index.intersection(second_df.index)
    first_df = first_df.filter(indexes, axis=0)
    second_df = second_df.filter(indexes, axis=0)

    vo = np.hstack((first_df._values, second_df._values))
    df = pd.DataFrame(data=vo, columns=["x", "y"])
    plt.plot('x', 'y', data=df, linestyle='', marker='o', markersize=0.2)
    plt.savefig("compare_2")


def main():
    first_file, second_file, output_dir, chromosome = format_args()
    chromosome = "chr%s" % chromosome
    # compare_samples_covariance(first_file, second_file, chromosome)


if __name__ == '__main__':
    main()
