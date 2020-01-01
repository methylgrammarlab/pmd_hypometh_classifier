#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python

import argparse
import glob
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from tqdm import tqdm, trange
import datetime
import collections
import numba


sys.path.append(os.path.dirname(os.getcwd()))
from commons import tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"  # TODO remove the mini
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
HISTOGRAM_FORMAT = "histogram_%s_%s"
SPLIT = 10000


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_histogram(series, patient, chromosome, num_of_bins, output):
    series.plot.hist(bins=num_of_bins)
    plt.xlabel('number of cells covering both location pairs')
    plt.xticks(range(num_of_bins))
    plt.title(HISTOGRAM_FORMAT % (patient, chromosome))
    plt.savefig(os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome)))
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
    series_val_ind = np.where(series == 1)
    masked_matrix = matrix[series_val_ind[0], :]
    return np.sum(masked_matrix, axis=0)


def work_counter(col_coverage, counter, amount_of_samples):
    for i in range(amount_of_samples + 1):
        counter[i] += np.where(col_coverage == i)[0].size
    return counter


def create_pairwise_coverage(cpg_format_file, output):
    """
    Creates a matrix where each cell holds the number of cells both locations covered
    note: The correlation function used for the calculations returns 1 on th diagonal, but the diagonal isn't used so can be ignored.
    :param cpg_format_file:
    :return:
    """
    tqdm.pandas()
    df = pd.read_pickle(cpg_format_file)
    converted_matrix = np.where(~np.isnan(df), 1, 0)
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
    counter = {}
    amount_of_samples = converted_matrix.shape[0]
    for i in range(amount_of_samples + 1):
        counter[i] = 0

    for col in trange(converted_matrix.shape[1], desc='location progress'):
        col_coverage = compare_matrix(converted_matrix[:, col], converted_matrix)
        counter = work_counter(col_coverage, counter, amount_of_samples)
    tools.counter_to_csv(counter, os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome)))

    # coverage_matrix = df.progress_apply(lambda col: compare_matrix(col, df), axis=0)
    # pairwise_coverage = coverage_matrix.where(
    #     np.triu(np.ones(coverage_matrix.shape), k=0).astype(bool)).stack().reset_index()
    # create_histogram(pairwise_coverage.loc[:, 0], patient, chromosome, df.shape[0], output)


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
        t1 = datetime.datetime.now()
        create_pairwise_coverage(file, output)
        t2 = datetime.datetime.now()
        print('time: ', t2 - t1)


if __name__ == '__main__':
    main()

"""/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin/python /cs/usr/liorf/PycharmProjects/proj_scwgbs/format_files/coverage_between_pairs.py --cpg_format_folder /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/cpg_format/threshold/CRC09/ --output_folder /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/stats/coverage_histograms/CRC09/"""
