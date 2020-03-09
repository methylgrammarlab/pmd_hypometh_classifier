#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python
import argparse
import collections
import copy
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
# from commons import data_tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_chr%d.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CSV_FILE = "average_methylation_of_nc_%s.csv"
PATIENTS = ['CRC01', 'CRC11']  # 'CRC04']#, 'CRC10', 'CRC13' ]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def average_nc_methylation(cpg_format_file):
    df = pd.read_pickle(cpg_format_file)
    normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    normal_df = df.loc[normal_cell_ids, :]
    average = np.mean(normal_df, axis=0)
    return average.dropna().values


def avg_main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    cpg_format_file_path = os.path.join(args.cpg_format_folder, CPG_FORMAT_FILE_FORMAT)
    all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    counter = collections.Counter()

    patient = os.path.split(args.cpg_format_folder)[-1]

    # for file in tqdm(all_cpg_format_file_paths, desc='files'):
    for file in tqdm(all_cpg_format_file_paths):
        average = average_nc_methylation(file)
        counter.update(average)
    data_tools.counter_to_csv(counter, os.path.join(output, CSV_FILE % patient))


def methylation_diff(patients_list, f):
    all_met_ind = []

    not_nan = None
    for patient in patients_list:
        df = pd.read_pickle(patient)
        normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        normal_df = df.loc[normal_cell_ids, :]
        average = np.mean(normal_df, axis=0)
        met_ind = average.index[np.where(average >= 0.6)[0]]
        if not_nan is None:
            not_nan = average.index[np.where(average.notnull())[0]]
        else:
            not_nan &= average.index[np.where(average.notnull())[0]]
        all_met_ind.append(met_ind)

    if len(all_met_ind) == 0:
        return None

    met_and = copy.deepcopy(all_met_ind[0])
    met_or = copy.deepcopy(all_met_ind[0])
    for i in range(1, len(all_met_ind)):
        met_ind = copy.deepcopy(all_met_ind[i])
        met_and &= met_ind
        met_or |= met_ind

    return len(met_and & not_nan) / len(met_or & not_nan) * 100


def diff_main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    f = open(os.path.join(output, "avg_methylation_60.out"), "w")

    all_cpg_format_file_paths = dict()
    for chr in range(1, 23):
        all_cpg_format_file_paths[chr] = []

    for patient in PATIENTS:
        for chr in range(1, 23):
            cpg_format_file_path = os.path.join(args.cpg_format_folder, patient, CPG_FORMAT_FILE_FORMAT % chr)
            all_cpg_format_file_paths[chr] += glob.glob(cpg_format_file_path)

    for chr in all_cpg_format_file_paths:
        v = methylation_diff(all_cpg_format_file_paths[chr], f)
        # f.write("chr%s:%s\n" % (chr, v))
        print(v)

    f.close()


if __name__ == '__main__':
    # avg_main()
    diff_main()
