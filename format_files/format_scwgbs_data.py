"""
Format the scwgbs data (txt files) to a pandas data frames where each file contains one chromosome for one
patient
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas
from tqdm import tqdm

sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.getcwd()))
from commons import files_tools as tools

OUTPUT_FILE_FORMAT = "{patient}_{cell}_{num}_{chromosome}_hg19.pickle.zlib"

# The format of the input files and the way to extract data from them
SCWGBS_FILE_RE = "*.singleC.cpg.txt"
SCWGBS_FILENAME_DETAILS_RE = re.compile(".+(CRC\d+)_(\w+)_(\d+).singleC.cpg")

# Indexes match the single cell columns in the file
CHR_INDEX = 0
START_INDEX = 1
STRAND_INDEX = 3
TOTAL_INDEX = 4
MET_INDEX = 5


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--patient', help='Name of the patient, all if not provided', required=False)
    parser.add_argument('--validate', help='Just validate the information', required=False,
                        action='store_true', default=False)
    parser.add_argument('--sort', action='store_true', default=False)
    args = parser.parse_args()

    output = args.output_folder
    files_path = os.path.join(args.raw, SCWGBS_FILE_RE)
    if args.patient:
        files_path = os.path.join(args.raw, "*%s" % args.patient + SCWGBS_FILE_RE)

    return output, files_path, args.validate, args.sort


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
    Note: this is really slow
    :param chr_dict: A dict where each key is chr and the value is a list where each value is a line in the
    file
    :type chr_dict: dict
    :return: A dict where each key is chr and the value is a np array where each value is a line in the file
    """
    chr_dict_new = {}  # Can't change dict while running on it

    # Go over the chr, create new array with the size of the unique position and sum the positions we saw
    # twice
    for chr in tqdm(chr_dict.keys(), desc="chr"):
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

        once_pos = np.isin(temp_array[:, 0], only_once)
        chr_array[i:] = temp_array[once_pos]
        chr_dict_new[chr] = chr_array[chr_array[:, 0].argsort()]

    return chr_dict_new


def get_patients_files_dict(files_path):
    """
    Get all the files and save them in a dictionary by the patient
    :param files_path: The path of the files to read
    :type files_path: str
    :return: A dict where each key is a patient and the values is a list of pathes
    :rtype: dict{str:list[str]}
    """
    all_file_paths = glob.glob(files_path)
    patient_dict = {}
    for file_path in all_file_paths:
        name, cell, num = SCWGBS_FILENAME_DETAILS_RE.findall(file_path)[0]
        if name not in patient_dict:
            patient_dict[name] = []

        patient_dict[name].append(file_path)

    return patient_dict


def validate_only(files_to_validate_folder, patient_dict):
    """
    Validate that we have all the files
    :param files_to_validate_folder: path to the folder which should contain all the files
    :param patient_dict: A dictionary with all the original files files
    """
    with open("log_file.txt", "w") as log_file:
        for patient in patient_dict:
            for file_path in patient_dict[patient]:
                patient, cell, num = SCWGBS_FILENAME_DETAILS_RE.findall(file_path)[0]
                chr_dict = extract_cols(file_path)

                for chromosome in chr_dict:
                    file_name = OUTPUT_FILE_FORMAT.format(**{"patient": patient, "cell": cell, "num": num,
                                                             "chromosome": chromosome})
                    output_path = os.path.join(files_to_validate_folder, patient, file_name)
                    if not os.path.exists(output_path):
                        log_file.write("%s,%s\n" % (patient, file_path))
                        break


def re_sort_files(output_files_folder):
    """
    Needed this to re-sort the data, first version wasn't sort
    :param output_files_folder: The folder of all the files - expected to have the following tree; output\patient\all_files
    """
    all_files_path = os.path.join(output_files_folder, "*", "*.pickle.zlib")
    all_files = glob.glob(all_files_path)
    for file_path in all_files:
        data = tools.load_compressed_pickle(file_path)
        data = data[data[:, 0].argsort()]
        tools.save_as_compressed_pickle(file_path, data)


def main():
    output, files_path, validate, sort = format_args()
    patient_dict = get_patients_files_dict(files_path)

    if validate:
        validate_only(output, patient_dict)

    elif sort:
        re_sort_files(output)

    else:
        pbar = tqdm(patient_dict.keys())
        for patient in patient_dict:
            pbar.set_description("Patients")
            fbar = tqdm(patient_dict[patient])
            for file_path in patient_dict[patient]:
                patient, cell, num = SCWGBS_FILENAME_DETAILS_RE.findall(file_path)[0]
                fbar.set_description("Files")

                # Skip files already parsed
                if os.path.exists(os.path.join(output, patient,
                                               OUTPUT_FILE_FORMAT % (patient, cell, num, "chr16"))):
                    continue

                chr_dict = format_scwgbs_file(file_path)

                # Save a file for each chromosome
                for chr in chr_dict:
                    file_name = OUTPUT_FILE_FORMAT % (patient, cell, num, chr)

                    output_folder = os.path.join(output, patient)
                    if not os.path.exists(output_folder):
                        os.mkdir(output_folder)

                    output_path = os.path.join(output_folder, file_name)
                    tools.save_as_compressed_pickle(output_path, chr_dict[chr])


if __name__ == '__main__':
    main()
