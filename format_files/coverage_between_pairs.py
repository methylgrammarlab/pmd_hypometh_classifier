import argparse
import glob
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s_30mini.dummy.pkl.zip"  # TODO remove the mini
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+)_30mini.dummy.pkl.zip")
HISTOGRAM_FORMAT = "histogram_%s_%s_part_%s"  # TODO remove the mini
SPLIT = 10


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


def create_histogram(series, patient, chromosome, num_of_bins, part, output):
    series.plot.hist(bins=num_of_bins)
    plt.xlabel('number of cells covering both location pairs')
    plt.xticks(range(num_of_bins))
    plt.title(HISTOGRAM_FORMAT % (patient, chromosome, part))
    plt.savefig(os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome, part)))
    plt.show()


def create_pairwise_coverage(cpg_format_file, output):
    """
    Creates a matrix where each cell holds the number of cells both locations covered
    note: The correlation function used for the calculations returns 1 on th diagonal, but the diagonal isn't used so can be ignored.
    :param cpg_format_file:
    :return:
    """
    df = pd.read_pickle(cpg_format_file)
    total_hist = pd.Series()
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
    for i in range(0, df.shape[1], SPLIT):
        coverage_matrix = df.iloc[:, i:i + SPLIT].corr(method=count_similar)
        pairwise_coverage = coverage_matrix.where(
            np.triu(np.ones(coverage_matrix.shape), k=1).astype(bool)).stack().reset_index()
        total_hist = pd.concat([total_hist, pairwise_coverage.loc[:, 0]], ignore_index=True)
        create_histogram(pairwise_coverage.loc[:, 0], patient, chromosome, df.shape[0], str(int(i / SPLIT)), output)
    create_histogram(total_hist, patient, chromosome, df.shape[0], 'all', output)


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
        create_pairwise_coverage(file, output)


if __name__ == '__main__':
    main()


"""venv/bin/python format_files/coverage_between_pairs.py --cpg_format_folder format_files/example_files/ --output_folder format_files/example_files/"""