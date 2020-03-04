#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python
import argparse
import glob
import os
import pickle
import re
import sys
import collections

import numpy as np
import pandas as pd
from tqdm import tqdm

import commons.slurm_tools

sys.path.append(os.path.dirname(os.getcwd()))
from commons import files_tools, data_tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_chr%d.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CSV_FILE = "average_methylation_of_nc_%s.csv"
PATIENTS = ['CRC01', 'CRC04', 'CRC10', 'CRC11', 'CRC13']


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


def methylation_diff(patients_list):
    all_met_ind =  []

    for patient in patients_list:
        df = pd.read_pickle(patient)
        normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        normal_df = df.loc[normal_cell_ids, :]
        average = np.mean(normal_df, axis=0)
        met_ind = set(average.index[np.where(average >= 0.75)[0]])
        all_met_ind.append(met_ind)

    met_and = all_met_ind[0]
    met_or = all_met_ind[0]
    for i in range(1, len(all_met_ind)):
        met_and &= all_met_ind[i]
        met_or |= all_met_ind[i]
    pass


def diff_main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    all_cpg_format_file_paths = dict()
    for chr in range(1, 23):
        all_cpg_format_file_paths[chr] = []

    for patient in PATIENTS:
        for chr in range(1, 23):
            cpg_format_file_path = os.path.join(args.cpg_format_folder, patient, CPG_FORMAT_FILE_FORMAT % chr)
            all_cpg_format_file_paths[chr] += glob.glob(cpg_format_file_path)

    for chr in all_cpg_format_file_paths:
        methylation_diff(all_cpg_format_file_paths[chr])


if __name__ == '__main__':
    # avg_main()
    diff_main()
