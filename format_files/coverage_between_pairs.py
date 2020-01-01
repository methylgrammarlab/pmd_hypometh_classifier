#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python

import argparse
import glob
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s_mini.dummy.pkl.zip"  # TODO remove the mini
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+)_mini.dummy.pkl.zip")
HISTOGRAM_FORMAT = "histogram_%s_%s_part_%s"
SPLIT = 10000


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_histogram(series, patient, chromosome, num_of_bins, part, output):
    series.plot.hist(bins=num_of_bins)
    plt.xlabel('number of cells covering both location pairs')
    plt.xticks(range(num_of_bins))
    plt.title(HISTOGRAM_FORMAT % (patient, chromosome, part))
    plt.savefig(os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome, part)))
    # plt.show()
    plt.close()


def count_similar(loc_1, loc_2):
    """
    Using np.where, for each series, returns the indexes of locations that aren't nan. Next, using np.intersect1d counts
    how many indices appeared in both series (because indices are unique counts how many times both indices had a value.
    :param loc_1: first series
    :param loc_2: second series
    :return: number of locations both series had a value(not nan).
    """
    loc_1_value_inds = np.where(~np.isnan(loc_1))
    loc_2_value_inds = np.where(~np.isnan(loc_2))
    return np.intersect1d(loc_1_value_inds, loc_2_value_inds).size


def compare_matrix(series, matrix):
    return matrix.apply(lambda col: count_similar(col, series), axis=0)
    pass
    # loc_1_val = np.where(~np.isnan(series))
    # loc_2_val = np.where(~np.isnan(matrix))


def create_pairwise_coverage(cpg_format_file, output):
    """
    Creates a matrix where each cell holds the number of cells both locations covered
    note: The correlation function used for the calculations returns 1 on th diagonal, but the diagonal isn't used so can be ignored.
    :param cpg_format_file:
    :return:
    """
    tqdm.pandas()
    df = pd.read_pickle(cpg_format_file)
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
    coverage_matrix = df.progress_apply(lambda col: compare_matrix(col, df), axis=0)
    pairwise_coverage = coverage_matrix.where(
        np.triu(np.ones(coverage_matrix.shape), k=1).astype(bool)).stack().reset_index()
    create_histogram(pairwise_coverage.loc[:, 0], patient, chromosome, df.shape[0], 'all', output)

    # total_hist = pd.Series()
    # patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
    # count_similar(df.iloc[:, 0], df.iloc[:, 1])
    # for i in tqdm(range(0, df.shape[1], SPLIT), desc='section of file'):
    #     coverage_matrix = df.iloc[:, i:i + SPLIT].corr(method=count_similar)
    #     pairwise_coverage = coverage_matrix.where(
    #         np.triu(np.ones(coverage_matrix.shape), k=1).astype(bool)).stack().reset_index()
    #     total_hist = pd.concat([total_hist, pairwise_coverage.loc[:, 0]], ignore_index=True)
    #     create_histogram(pairwise_coverage.loc[:, 0], patient, chromosome, df.shape[0], str(int(i / SPLIT)), output)
    # create_histogram(total_hist, patient, chromosome, df.shape[0], 'all', output)


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
    for file in tqdm(all_cpg_format_file_paths, desc='files'):
        create_pairwise_coverage(file, output)


if __name__ == '__main__':
    print("hey")
    main()

"""/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin/python /cs/usr/liorf/PycharmProjects/proj_scwgbs/format_files/coverage_between_pairs.py --cpg_format_folder /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/cpg_format/threshold/CRC09/ --output_folder /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/stats/coverage_histograms/CRC09/"""
