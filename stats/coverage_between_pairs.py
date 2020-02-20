#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python
import argparse
import glob
import os
import pickle
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import commons.slurm_tools

sys.path.append(os.path.dirname(os.getcwd()))
from commons import files_tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
HISTOGRAM_FORMAT = "histogram_%s_%s.csv"
PICKLE_FORMAT = "dict_%s_%s.pickle"
WINDOWS_SIZE = 500


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


def compare_matrix(matrix, index):
    """
    Counts per CpG for how many of the samples it had a read.
    :param matrix:
    :param index:
    :return:
    """
    series = matrix[:, index]
    series_val_ind = np.where(series == 1)[0]
    masked_matrix = matrix[series_val_ind, :]
    return np.sum(masked_matrix, axis=0)


def work_counter(col_coverage, i):
    return np.where(col_coverage == i)[0].size


def create_pairwise_coverage(cpg_format_file):
    """
    Creates a matrix where each cell holds the number of cells both locations covered
    note: The correlation function used for the calculations returns 1 on th diagonal, but the diagonal isn't used so can be ignored.
    :param cpg_format_file:
    :return:
    """
    # tqdm.pandas()
    df = pd.read_pickle(cpg_format_file)
    converted_matrix = np.where(~np.isnan(df), 1, 0)

    counter = {}
    amount_of_samples = converted_matrix.shape[0]
    for i in range(amount_of_samples + 1):
        counter[i] = 0

    # for col in trange(converted_matrix.shape[1], desc='location progress'):
    for col in range(converted_matrix.shape[1]):
        col_coverage = compare_matrix(converted_matrix, col)
        for i in range(amount_of_samples + 1):
            c = work_counter(col_coverage, i)
            counter[i] += c

    return counter


def nc_pairwise_coverage(cpg_format_file):
    df = pd.read_pickle(cpg_format_file)
    normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    normal_df = df.loc[normal_cell_ids, :]
    converted_matrix = np.where(~np.isnan(normal_df), 1, 0)


    amount_of_samples = converted_matrix.shape[0]
    num_of_cpg = 3000
    # num_of_cpg = converted_matrix.shape[1]
    counter_matrix = np.empty((num_of_cpg, amount_of_samples + 1))
    counter = {}
    for i in range(amount_of_samples + 1):
        counter[i] = 0
    for window in tqdm(range(0, num_of_cpg, WINDOWS_SIZE)):
        window_matrix = converted_matrix[:, window:min(window + WINDOWS_SIZE, num_of_cpg)]
        for col in range(window_matrix.shape[1]):
            col_coverage = compare_matrix(window_matrix, col)

            for i in range(amount_of_samples + 1):
                c = work_counter(col_coverage, i)
                counter_matrix[window+col, i] = c
                # todo add another column of number of reads in the nc - just sum the rows in the converted matrix and add as a column.
                counter[i] += c

    return counter


def write_output(chromosome, counter, output, patient):
    output_path = os.path.join(output, patient, PICKLE_FORMAT % (patient, chromosome))
    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    with open(output_path, "wb") as of:
        pickle.dump(counter, of)


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
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file)[0]
        counter = nc_pairwise_coverage(file)
        # counter = create_pairwise_coverage(file)
        write_output(chromosome, counter, output, patient)


if __name__ == '__main__':
    commons.slurm_tools.init_slurm(main)
