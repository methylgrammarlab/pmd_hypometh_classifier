import argparse
import glob
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s_mini.dummy.pkl.zip"  # TODO remove the mini


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def count_similar(loc_1, loc_2):
    loc_1_val = np.where(~np.isnan(loc_1))
    loc_2_val = np.where(~np.isnan(loc_2))
    return np.intersect1d(loc_1_val, loc_2_val).size


def create_pairwise_coverage(cpg_format_file):
    """
    Creates a matrix where each cell holds the number of cells both locations covered
    note: The correlation function used for the calculations returns 1 on th diagonal, but the diagonal isn't used so can be ignored.
    :param cpg_format_file:
    :return:
    """
    df = pd.read_pickle(cpg_format_file)
    return df.corr(method=count_similar)


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    chr = args.chr

    cpg_format_file_path = os.path.join(args.cpg_format_folder, CPG_FORMAT_FILE_FORMAT % '*')
    if chr:
        cpg_format_file_path = os.path.join(args.cpg_format_folder, CPG_FORMAT_FILE_FORMAT % chr)

    all_cpg_format_file_paths = glob.glob(cpg_format_file_path)
    for file in all_cpg_format_file_paths:
        coverage_matrix = create_pairwise_coverage(file)
        pairwise_coverage = coverage_matrix.where(
            np.triu(np.ones(coverage_matrix.shape), k=1).astype(bool)).stack().reset_index()
        pairwise_coverage.loc[:, 0].plot.hist(bins=27)
        plt.xlabel('number of cells covering both location pairs')
        plt.savefig('histogram.png')
        plt.show()
        pass


if __name__ == '__main__':
    main()
