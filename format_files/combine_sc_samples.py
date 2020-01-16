import argparse
import glob
import os
import re
import sys

from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

import numpy as np
import pandas as pd

import commons.files_tools as tools

PATIENT_FILE_FORMAT = "CRC*.pickle.zlib"
FULL_CPG_FILE_FORMAT = "full_cpg_seq_chr*.pickle.zlib"
OUTPUT_FILE_FORMAT = "all_cpg_ratios_%s_%s.dummy.pkl.zip"
PATIENT_FILE_DETAILS_RE = re.compile(".+(CRC\d+)_.+_(chr\d+|chrX|chrY).pickle.zlib")
PATIENT_CELL_NAME_RE = re.compile(".+(CRC\d+)_(\w+_\d+)_(chr\d+)|chrX|chrY.pickle.zlib")
FULL_CPG_FILE_DETAILS_RE = re.compile(".+full_cpg_seq_(chr\d+|chrX|chrY).pickle.zlib")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sc_folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--genomic_folder', help='Path to folder of genomic data', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--blacklisted_files', help='Path of file with list of files to filter out', required=False)
    parser.add_argument('--coverage_threshold', help='Number of reads to include (inclusive)', type=int, required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_chr_df(chr_file_list, chr_cpg_pos, threshold):
    chr_full_cpg = chr_cpg_pos[:, 0]

    chr_all_cells = np.full((len(chr_file_list) + 1, chr_full_cpg.size), None, dtype=np.float)
    chr_all_cells[0, :] = chr_full_cpg
    cell_names = []

    i = 1
    for cell_path in chr_file_list:
        cell_names.append(PATIENT_CELL_NAME_RE.findall(cell_path)[0][1])
        cell = tools.load_compressed_pickle(cell_path)
        if threshold:
            indices = np.where(cell[:,1] <= threshold)
            cell = cell[indices]
        match = np.isin(chr_full_cpg, cell[:, 0])
        ratio = cell[:, 2] / cell[:, 1]
        try:
            chr_all_cells[i, match] = ratio
        except Exception:
            cell_names[-1] = "missing_data"
            with open("missing_data_%s.txt" % PATIENT_CELL_NAME_RE.findall(cell_path)[0][0], "a") as missing:
                missing.write("Patient: %s, Cell: %s, Chromosome: %s\n" % PATIENT_CELL_NAME_RE.findall(cell_path)[0])
        i += 1

    df = pd.DataFrame(data=chr_all_cells[1:, :], columns=chr_all_cells[0, :].astype(np.int), index=cell_names,
                      dtype=np.float16)
    return df


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    # get dict chr#: array[cpg]
    full_cpg_files_path = os.path.join(args.genomic_folder, FULL_CPG_FILE_FORMAT)
    all_full_cpg_file_paths = glob.glob(full_cpg_files_path)
    chr_pos_dict = {}
    for file_path in all_full_cpg_file_paths:
        chr = FULL_CPG_FILE_DETAILS_RE.findall(file_path)[0]
        if chr not in chr_pos_dict:
            chr_pos_dict[chr] = []

        full_cpg_seq = tools.load_compressed_pickle(file_path)
        chr_pos_dict[chr].append(full_cpg_seq)

    # patient files path
    patient_files_path = os.path.join(args.sc_folder, PATIENT_FILE_FORMAT)
    all_patient_file_paths = glob.glob(patient_files_path)
    patient_chr_dict = {}

    # files to filter
    if args.blacklisted_files:
        with open(args.blacklisted_files) as f:
            blacklisted_files = set(os.path.join(args.sc_folder, line.rstrip()) for line in f)
        all_patient_file_paths = set(all_patient_file_paths) - blacklisted_files


    # get dict patient#chr#: array[path]
    for file_path in all_patient_file_paths:
        name = PATIENT_FILE_DETAILS_RE.findall(file_path)[0]
        if name[1].endswith('X') or name[1].endswith('Y'):
            continue
        if name not in patient_chr_dict:
            patient_chr_dict[name] = []

        patient_chr_dict[name].append(file_path)

    # get pandas df
    for patient_and_chr in tqdm(patient_chr_dict, desc="chr"):
        if args.chr and args.chr != patient_and_chr[1]:
            continue
        chr_all_cells = create_chr_df(patient_chr_dict[patient_and_chr], chr_pos_dict[patient_and_chr[1]][0],
                                      args.coverage_threshold)
        output_path = os.path.join(output, OUTPUT_FILE_FORMAT % (patient_and_chr[0], patient_and_chr[1]))
        chr_all_cells.to_pickle(output_path, compression='zip')


if __name__ == '__main__':
    main()
