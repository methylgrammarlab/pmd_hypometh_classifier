import argparse
import gc
import glob
import os
import re
import sys
import pandas

import numpy as np

import format_files.tools as tools

OUTPUT_FILE_FORMAT = "%s_%s_%s_%s.pickle.zlib"

# The format of the input files and the way to extract data from them
FILE_SUFFIX = "*.singleC.cpg.txt"
FILE_DETAILS_RE = re.compile(".+(CRC\d+)_(\w+)_(\d+).singleC.cpg")

# Indexes match the single cell columns in the file
CHR_INDEX = 0
POS_INDEX = 1
STRAND_INDEX = 3
TOTAL_INDEX = 4
MET_INDEX = 5


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--patient', help='Name of the patient, all if not provided', required=False)
    parser.add_argument('--validate', help='Just validate the information', required=False,
                        action='store_true', default=False)
    args = parser.parse_args()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    files_path = os.path.join(args.raw, FILE_SUFFIX)
    if args.patient:
        files_path = os.path.join(args.raw, "*%s" % args.patient + FILE_SUFFIX)

    return output, files_path, args.validate


def format_scwgbs_file(file_path):
    """
    Format a scwgbs file to a more usable manner
    :param file_path: The path of the file to format
    :type file_path: str
    :return: A dict where each key is a chr and the value is an array with all the scwgbs reads
    :rtype: dict
    """
    chr_dict = extract_cols(file_path)
    chr_dict = combine_strands(chr_dict)
    return chr_dict


def extract_cols(file_path):
    """
    Extract the relevant columns from the file and save them by chromosome
    :param file_path: The path of the file to read
    :type file_path: str
    :return: A dict where each key is chr and the value is a list where each value is a line in the file
    :rtype: dict[str:list]
    """
    chr_dict = {}
    with open(file_path) as f:
        pandas_file = pandas.read_csv(f, sep="\t")
        data = pandas_file._values
        minus_index = np.where(data[:, 3] == "-")
        data[minus_index, 1] = data[minus_index, 1] - 1
        chrs = np.unique(data[:, 0])
        for chr in chrs:
            data_to_dict = data[np.where(data[:, 0] == chr), :][0][:, [1, 4, 5]]
            chr_dict[chr] = data_to_dict.astype(dtype=np.uint32)

        return chr_dict


def combine_strands(chr_dict):
    """
    Combine the information from the + strand with the information from the - strand
    :param chr_dict: A dict where each key is chr and the value is a list where each value is a line in the
    file
    :type chr_dict: dict
    :return: A dict where each key is chr and the value is a np array where each value is a line in the file
    """
    chr_dict_new = {}  # Can't change dict while running on it

    # Go over the chr, create new array with the size of the unique position and sum the positions we saw
    # twice
    for chr in chr_dict:
        temp_array = chr_dict[chr]
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


def get_patients_files_dict(files_path):
    """
    Get all the files and save them by the patient
    :param files_path: The path of the files to read
    :type files_path: str
    :return: A dict where each key is a patient and the values is a list of pathes
    :rtype: dict{str:list[str]}
    """
    all_file_paths = glob.glob(files_path)
    patient_dict = {}
    for file_path in all_file_paths:
        name, cell, num = FILE_DETAILS_RE.findall(file_path)[0]
        if name not in patient_dict:
            patient_dict[name] = []

        patient_dict[name].append(file_path)

    return patient_dict


def validate_only(output, patient_dict):
    log_file = open("log_file.txt", "w")
    for patient in patient_dict:
        for file_path in patient_dict[patient]:
            patient, cell, num = FILE_DETAILS_RE.findall(file_path)[0]
            chr_dict = extract_cols(file_path)

            for chr in chr_dict:
                file_name = OUTPUT_FILE_FORMAT % (patient, cell, num, chr)
                output_path = os.path.join(output, patient, file_name)
                if not os.path.exists(output_path):
                    log_file.write("%s,%s\n" %(patient, file_path))

    log_file.close()


def main():
    output, files_path, validate = format_args()
    patient_dict = get_patients_files_dict(files_path)

    if validate:
        validate_only(output, patient_dict)

    for patient in patient_dict:
        for file_path in patient_dict[patient]:
            patient, cell, num = FILE_DETAILS_RE.findall(file_path)[0]
            expected_file = OUTPUT_FILE_FORMAT % (patient, cell, num, "chr16")
            if os.path.exists(os.path.join(output, expected_file)):
                continue

            print("here")
            chr_dict = format_scwgbs_file(file_path)

            for chr in chr_dict:
                file_name = OUTPUT_FILE_FORMAT % (patient, cell, num, chr)
                output_path = os.path.join(output, file_name)
                tools.save_as_compressed_pickle(output_path, chr_dict[chr])

        gc.collect()


if __name__ == '__main__':
    main()
