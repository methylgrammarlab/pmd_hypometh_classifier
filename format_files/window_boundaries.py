import argparse
import copy
import glob
import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools

SCWGBS_FILE_FORMAT = "all_cpg_ratios_*_chr%d.dummy.pkl.zip"
PATIENTS = ['CRC01', 'CRC13', 'CRC11']
BEDGRAPH_FILE_FORMAT = "window_boundries.bedgraph"
PICKLE_FILE_FORMAT = "window_boundries.dummy.pkl.zip"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_folder', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def choose_cpgs(files):
    all_nan_cpgs = []
    columns = []
    for file in files:
        df = pd.read_pickle(file)
        columns.append(df.columns)
        tumour_cell_ids = [cell_id for cell_id in df.index if not cell_id.startswith('NC')]
        tumour_df = df.loc[tumour_cell_ids, :]
        masked_tumour = np.where(np.isnan(tumour_df), 0, 1)
        nan_cpgs = np.sum(masked_tumour, axis=0)
        all_nan_cpgs.append(tumour_df.columns[nan_cpgs == 0])

    if len(all_nan_cpgs) == 0:
        return

    met_or = copy.deepcopy(all_nan_cpgs[0])
    for i in range(1, len(all_nan_cpgs)):
        met_ind = copy.deepcopy(all_nan_cpgs[i])
        met_or |= met_ind
    pass


def window_boundries(files):
    if len(files) == 0:
        return []
    df = pd.read_pickle(files[0])
    num_of_cpg = df.shape[1]
    chr_windows = []
    window_size = 5000
    index_advance_value = 50000

    for i in tqdm(range(0, num_of_cpg, index_advance_value)):
        window_indices = (df.columns[i], df.columns[min(i + window_size, num_of_cpg) - 1])  # Get the indexes
        chr_windows.append(window_indices)
    return chr_windows


def save_windows(all_window_boundries, output_folder):
    pickle_output_path = os.path.join(output_folder, PICKLE_FILE_FORMAT)
    bedgraph_output_path = os.path.join(output_folder, BEDGRAPH_FILE_FORMAT)
    files_tools.save_as_compressed_pickle(pickle_output_path, all_window_boundries)
    index = []
    for key in all_window_boundries:
        index += [key for i in range(len(all_window_boundries[key]))]

    columns = ['start', 'end', 'val']
    df = pd.DataFrame(columns=columns, index=index)
    ind = 0
    for key in all_window_boundries:
        for tup in all_window_boundries[key]:
            start = tup[0]
            df.iloc[ind, :] = [start, start + 1, 1]
            ind += 1
    df.to_csv(bedgraph_output_path, header=False, sep='\t')


def main():
    args = parse_input()

    all_cpg_format_file_paths = dict()
    for chr in range(1, 23):
        all_cpg_format_file_paths[chr] = []

    for patient in PATIENTS:
        for chr in range(1, 23):
            cpg_format_file_path = os.path.join(args.cpg_format_folder, patient, SCWGBS_FILE_FORMAT % chr)
            all_cpg_format_file_paths[chr] += glob.glob(cpg_format_file_path)

    all_window_boundries = {}

    for chr in tqdm(all_cpg_format_file_paths):
        # choose_cpgs(all_cpg_format_file_paths[chr])
        all_window_boundries[chr] = window_boundries(all_cpg_format_file_paths[chr])
    save_windows(all_window_boundries, args.output_folder)


if __name__ == '__main__':
    main()
