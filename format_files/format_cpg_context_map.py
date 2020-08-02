"""
Format a chromosome files to numpy object for future handling
"""

import argparse
import csv
import glob
import itertools
import os
import sys

import numpy as np
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import files_tools as tools
from commons import consts

INPUT_FILES_FORMAT = "*seq_context*"

# Indexes match the seq context columns in the file
POS_INDEX = 2
ORPH_INDEX_START = 3
ORPH_INDEX_END = 16
CONTEXT_INDEX = -1

NUMBER_OF_ORPH_PER_INDEX = [1, 2, 3, 4, 5, 10, 20, 35, 50, 75, 100, 150, 200]

CONTEXT_INT_TO_CHR_DICT = {}
all_possibilities = itertools.product(consts.NUCLEOTIDE_TO_NUMBER.keys(), repeat=8)


def get_orph_order():
    """
    Get the orph order in case someone needs it
    :return: A list where each value is the distance another CpG was searched
    :rtype: list[str]
    """
    return NUMBER_OF_ORPH_PER_INDEX


def is_weak(context):
    return False if consts.WEAK_FORMAT_RE.match(context) is None else True


def is_strong(context):
    return False if consts.STRONG_FORMAT_RE.match(context) is None else True


def convert_context_to_int(context):
    return "".join(consts.NUCLEOTIDE_TO_NUMBER[letter] for letter in context)


def convert_context_int_to_str(i):
    return CONTEXT_INT_TO_CHR_DICT[str(i)]


def get_context_as_int_for_chr(chr_info):
    return chr_info[:, -3]


def get_context_as_str_for_chr(chr_info):
    return [convert_context_int_to_str(i) for i in chr_info[:, -3]]


def get_context_as_str_for_chr_2(chr_info):
    return [convert_context_int_to_str(i)[1:-1] for i in chr_info[:, -3]]


def get_context_as_str_for_chr_1(chr_info):
    return [convert_context_int_to_str(i)[2:-2] for i in chr_info[:, -3]]


def get_weak_column(chr_info):
    return chr_info[:, -2]


def get_strong_column(chr_info):
    return chr_info[:, -1]


def get_orph_35_column(chr_info):
    return chr_info[:, -9]


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of chr_cpg', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()

    output_folder = args.output_folder
    files_paths = os.path.join(args.raw, INPUT_FILES_FORMAT)

    return output_folder, files_paths


def format_cpg_seq_file(cpg_seq_path):
    """
    Format a CpG seq file
    :param cpg_seq_path: The path of the cpg file
    :type cpg_seq_path: str
    :return: An array of the relevant data from a cpg seq file
    :rtype: np.ndarray
    """
    cpg_lines = []
    with open(cpg_seq_path, "r") as cpg_file:
        csv_file = csv.DictReader(cpg_file, delimiter="\t")
        for _ in csv_file.reader:
            break

        for line in csv_file.reader:
            pos, orph, context = line[POS_INDEX], line[ORPH_INDEX_START:ORPH_INDEX_END], line[CONTEXT_INDEX]
            orph.append(context)
            orph.insert(0, pos)
            cpg_lines.append(orph)

    cpg_array = np.array(cpg_lines, dtype=np.str)
    return cpg_array


def main():
    output_folder, files_path = format_args()
    all_file_paths = glob.glob(files_path)

    convert_context_to_int_vc = np.vectorize(convert_context_to_int)
    is_strong_vc = np.vectorize(is_strong)
    is_weak_vc = np.vectorize(is_weak)

    for cpg_file in tqdm(all_file_paths):
        cpg_matrix = format_cpg_seq_file(cpg_file)

        # Convert the contextto int and add columns of is_weak and is_strong
        context_col_int = convert_context_to_int_vc(cpg_matrix[:, -1])
        is_weak_col = is_weak_vc(cpg_matrix[:, -1]).reshape(cpg_matrix.shape[0], 1)
        is_strong_col = is_strong_vc(cpg_matrix[:, -1]).reshape(cpg_matrix.shape[0], 1)

        # Add new cols and convert matrix to u_int32
        cpg_matrix[:, -1] = context_col_int
        cpg_matrix = cpg_matrix.astype(np.uint32)
        cpg_matrix = np.hstack((cpg_matrix, is_weak_col))
        cpg_matrix = np.hstack((cpg_matrix, is_strong_col))

        # Save the file
        chr_number = consts.CHR_FULL_NAME_RE.findall(cpg_file)[0]
        output_path = os.path.join(output_folder, consts.FULL_CPG_CONTEXT_FILE_FORMAT % chr_number)
        tools.save_as_compressed_pickle(output_path, cpg_matrix)


if __name__ == '__main__':
    main()

for combination in all_possibilities:
    value = "".join(combination)
    CONTEXT_INT_TO_CHR_DICT[convert_context_to_int(value)] = value
