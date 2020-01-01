#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python

import argparse
import glob
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import datetime
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
from commons import tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"  # TODO remove the mini
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
HISTOGRAM_FORMAT = "histogram_%s_%s.csv"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


# should move to commons
def create_histogram(series, patient, chromosome, num_of_bins, output):
    series.plot.hist(bins=num_of_bins)
    plt.xlabel('number of cells covering both location pairs')
    plt.xticks(range(num_of_bins))
    plt.title(HISTOGRAM_FORMAT % (patient, chromosome))
    plt.savefig(os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome)))
    # plt.show()
    plt.close()


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
    # tqdm.pandas()
    df = pd.read_pickle(cpg_format_file)
    converted_matrix = np.where(~np.isnan(df), 1, 0)
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
    counter = {}
    amount_of_samples = converted_matrix.shape[0]
    for i in range(amount_of_samples + 1):
        counter[i] = 0

    # for col in trange(converted_matrix.shape[1], desc='location progress'):
    for col in range(converted_matrix.shape[1]):
        col_coverage = compare_matrix(converted_matrix[:, col], converted_matrix)
        counter = work_counter(col_coverage, counter, amount_of_samples)
    tools.counter_to_csv(counter, os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome)))


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    chr = args.chr

    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        if chr:
            cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % chr)

        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)
    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    # for file in tqdm(all_cpg_format_file_paths, desc='files'):
    for file in all_cpg_format_file_paths:
        create_pairwise_coverage(file, output)


if __name__ == '__main__':
    t1 = datetime.datetime.now()
    main()
    t2 = datetime.datetime.now()
    print('start time: ', t1, ' end time: ', t2)
