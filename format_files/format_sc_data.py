import argparse
import csv
import gc
import glob
import os
import pickle
import re
import sys
import zlib

import numpy as np

FILE_SUFFIX = "*.singleC.cpg.txt"
FILE_DETAILS_RE = re.compile(".+(CRC\d+)_(\w+)_(\d+).singleC.cpg")
CHR_INDEX = 0
POS_INDEX = 1
STRAND_INDEX = 3
TOTAL_INDEX = 4
MET_INDEX = 5


def input_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--patient', help='Name of the patient, all if not provided', required=False)
    args = parser.parse_args()
    return args


def format_file(file_path):
    chr_dict = extract_cols(file_path)
    chr_dict = combine_strands(chr_dict)
    return chr_dict


def extract_cols(file_path):
    chr_dict = {}
    with open(file_path) as f:
        csv_file = csv.DictReader(f, delimiter="\t")
        for _ in csv_file.reader:
            break

        for line in csv_file.reader:
            chr, pos, strand, total, met = line[CHR_INDEX], int(line[POS_INDEX]), line[STRAND_INDEX], \
                                           int(line[TOTAL_INDEX]), int(line[MET_INDEX])

            if strand == "-":
                pos -= 1

            if chr not in chr_dict:
                chr_dict[chr] = []

            chr_dict[chr].append([pos, total, met])

    return chr_dict


def combine_strands(chr_dict):
    # Combine same pos to one
    chr_dict_new = {}
    for chr in chr_dict:
        temp_array = np.array(chr_dict[chr], dtype=np.uint32)
        positions, positions_count = np.unique(temp_array[:, 0], return_counts=True)
        only_once = positions[np.where(positions_count == 1)]
        twice = positions[np.where(positions_count != 1)]
        chr_array = np.empty((positions.size, temp_array.shape[1]), dtype=np.uint32)
        i = 0
        for position in twice:
            arr = np.sum(temp_array[np.where(temp_array[:, 0] == position)], 0)
            arr[0] = position
            chr_array[i] = arr
            i += 1

        once_pos = np.in1d(temp_array[:, 0], only_once)
        chr_array[i:] = temp_array[once_pos]
        chr_dict_new[chr] = chr_array

    return chr_dict_new


def save_data(file_path, chr, chr_data, output_dir):
    patient, cell, num = FILE_DETAILS_RE.findall(file_path)[0]
    file_name = "%s_%s_%s_%s.pickle.zlib" % (patient, cell, num, chr)
    with open(os.path.join(output_dir, file_name), "wb") as patient_file:
        patient_file.write(zlib.compress(pickle.dumps(chr_data, pickle.HIGHEST_PROTOCOL), 9))


def main():
    args = input_parser()
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    files_path = os.path.join(args.raw, FILE_SUFFIX)
    if args.patient:
        files_path = os.path.join(args.raw, "*%s" % args.patient + FILE_SUFFIX)

    all_file_paths = glob.glob(files_path)
    patient_dict = {}
    for file_path in all_file_paths:
        name, cell, num = FILE_DETAILS_RE.findall(file_path)[0]
        if name not in patient_dict:
            patient_dict[name] = []

        patient_dict[name].append(file_path)

    for patient in patient_dict:
        for file_path in patient_dict[patient]:
            chr_dict = format_file(file_path)

            for chr in chr_dict:
                save_data(file_path, chr, chr_dict[chr], output)

        gc.collect()


if __name__ == '__main__':
    main()
