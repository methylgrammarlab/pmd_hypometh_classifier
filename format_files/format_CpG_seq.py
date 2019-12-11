import argparse
import sys
import os
import glob
import re
import numpy as np

import commons.tools as tools

PATIENT_FILE_FORMAT = "CRC*.pickle.zlib"
FULL_CPG_FILE_FORMAT = "full_cpg_seq_chr*.pickle.zlib"
OUTPUT_FILE_FORMAT = "full_cpg_ratios_patient%chr%s.pickle.zlib"
PATIENT_FILE_DETAILS_RE = re.compile(".+(CRC\d+)_.+_(chr\d+).pickle.zlib")
FULL_CPG_FILE_DETAILS_RE = re.compile(".+full_cpg_seq_(chr\d+).pickle.zlib")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--patient', help='Name of the patient, all if not provided', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided', required=False)
    args = parser.parse_args()
    return args


def create_chr_df(chr_file_list, chr_cpg_pos):
    chr_full_cpg = chr_cpg_pos[:, 0]

    chr_all_cells = np.full((len(chr_file_list) + 1, chr_full_cpg.size), None, dtype=np.float)
    chr_all_cells[0, :] = chr_full_cpg

    i = 1
    for chr_path in chr_file_list:
        cell = tools.load_compressed_pickle(chr_path)
        cell = cell[cell[:,0].argsort()]
        # cell_indexes = np.where(
        #     np.logical_and(cell[:, POSITION_INDEX] > start, cell[:, POSITION_INDEX] < end))
        # pmd = cell[cell_indexes]
        match = np.isin(chr_full_cpg, cell[:, 0])
        ratio = cell[:, 2] / cell[:, 1]
        chr_all_cells[i, match] = ratio
        i += 1

    return chr_all_cells


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    # get dict chr#: array[cpg]
    full_cpg_files_path = os.path.join(args.folder, FULL_CPG_FILE_FORMAT)
    all_full_cpg_file_paths = glob.glob(full_cpg_files_path)
    chr_pos_dict = {}
    for file_path in all_full_cpg_file_paths:
        chr = FULL_CPG_FILE_DETAILS_RE.findall(file_path)[0]
        if chr not in chr_pos_dict:
            chr_pos_dict[chr] = []

        full_cpg_seq = tools.load_compressed_pickle(file_path)
        chr_pos_dict[chr].append(full_cpg_seq)

    # patient files path
    patient_files_path = os.path.join(args.folder, PATIENT_FILE_FORMAT)
    if args.patient:
        patient_files_path = os.path.join(args.folder, "%s" % args.patient + PATIENT_FILE_FORMAT)
    all_patient_file_paths = glob.glob(patient_files_path)
    patient_chr_dict = {}

    # get dict patient#: array[cell]
    for file_path in all_patient_file_paths:
        name = PATIENT_FILE_DETAILS_RE.findall(file_path)[0]
        if name not in patient_chr_dict:
            patient_chr_dict[name] = []

        patient_chr_dict[name].append(file_path)

    for patient_chr in patient_chr_dict:
        chr_all_cells = create_chr_df(patient_chr_dict[patient_chr], chr_pos_dict[patient_chr[1]][0])
        tools.save_as_compressed_pickle()


if __name__ == '__main__':
    main()