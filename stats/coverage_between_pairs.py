# TODO:lior
# !/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python
import argparse
import glob
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

# Needed for imports
sys.path.append(os.path.dirname(os.getcwd()))
from commons import files_tools, consts, data_tools, slurm_tools

HISTOGRAM_FORMAT = "histogram_%s_%s.csv"
PICKLE_FORMAT = "dict_%s_%s.pickle"
NC_PICKLE_FORMAT = "NC_pairwise_coverage_%s_%s_numNC_%d_window_%d.dummy.pkl.zip"
WINDOWS_SIZE = 500


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


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
    df = pd.read_pickle(cpg_format_file)
    converted_matrix = np.where(~np.isnan(df), 1, 0)

    counter = {}
    amount_of_samples = converted_matrix.shape[0]
    for i in range(amount_of_samples + 1):
        counter[i] = 0

    for col in range(converted_matrix.shape[1]):
        col_coverage = compare_matrix(converted_matrix, col)
        for i in range(amount_of_samples + 1):
            c = work_counter(col_coverage, i)
            counter[i] += c

    return counter


def write_output(chromosome, counter, output, patient):
    output_path = os.path.join(output, patient, PICKLE_FORMAT % (patient, chromosome))
    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    with open(output_path, "wb") as of:
        pickle.dump(counter, of)


def nc_pairwise_coverage(cpg_format_file):
    df = pd.read_pickle(cpg_format_file)
    normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    normal_df = df.loc[normal_cell_ids, :]
    converted_matrix = np.where(~np.isnan(normal_df), 1, 0)

    amount_of_samples = converted_matrix.shape[0]
    num_of_cpg = converted_matrix.shape[1]
    counter_matrix = np.empty((num_of_cpg, amount_of_samples + 2))
    for window in tqdm(range(0, num_of_cpg, WINDOWS_SIZE)):
        window_matrix = converted_matrix[:, window:min(window + WINDOWS_SIZE, num_of_cpg)]
        for col in range(window_matrix.shape[1]):
            col_coverage = compare_matrix(window_matrix, col)

            for i in range(amount_of_samples + 1):
                c = work_counter(col_coverage, i)
                counter_matrix[window + col, i] = c
    counter_matrix[:, -1] = converted_matrix.sum(axis=0)
    counter_df = pd.DataFrame(counter_matrix, index=df.columns)
    return counter_df


def write_nc_output(chromosome, counter, output, patient):
    output_path = os.path.join(output, patient,
                               NC_PICKLE_FORMAT % (patient, chromosome, counter.shape[1] - 2, WINDOWS_SIZE))
    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    column_names = list(counter.columns)
    column_names[-1] = '#NC'
    counter.columns = column_names
    counter.to_pickle(output_path, compression='zip')


def main():
    args = parse_input()

    input_folder, output, chromosome = args.cpg_format_files, args.output_folder, args.chr

    if chromosome:
        formatted_pattern = consts.SCWGBS_FILE_FORMAT % ('*', chromosome)
    else:
        formatted_pattern = consts.SCWGBS_FILE_FORMAT % ('*', '*')

    all_cpg_format_file_paths = files_tools.get_files_to_work(input_folder, formatted_pattern)

    for file in all_cpg_format_file_paths:
        patient, chromosome = consts.DATA_FILE_SCWGBS_RE.findall(file)[0]
        # counter = create_pairwise_coverage(file)
        # write_output(chromosome, counter, output, patient)
        counter = nc_pairwise_coverage(file)
        write_nc_output(chromosome, counter, output, patient)


if __name__ == '__main__':
    slurm_tools.init_slurm(main)
