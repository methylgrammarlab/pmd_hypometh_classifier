import argparse
import csv
import glob
import os
import re
import sys

from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

import numpy as np

import commons.tools as tools

OUTPUT_FILE_FORMAT = "full_cpg_seq_chr%s.pickle.zlib"

GLOB_FORMAT = "*seq_context*"
FILE_FORMAT_RE = re.compile("chr(\d+).*")

# Regex for weak and strong context cpg
WEAK_FORMAT_RE = re.compile("\w\w[AT]CG[AT]\w\w")
STRONG_FORMAT_RE = re.compile("\w\w[CG]CG[CG]\w\w")

# Mapping between nucleotide and number
NUCLEOTIDE_TO_NUMBER = {"A": "1",
                        "C": "2",
                        "G": "3",
                        "T": "4",
                        "N": "5"}  # N means that this is unknown

# Indexes match the seq context columns in the file
POS_INDEX = 2
ORPH_INDEX_START = 3
ORPH_INDEX_END = 16
CONTEXT_INDEX = -1

NUMBER_OF_ORPH_PER_INDEX = [1, 2, 3, 4, 5, 10, 20, 35, 50, 75, 100, 150, 200]


def is_weak(context):
    return False if WEAK_FORMAT_RE.match(context) is None else True


def is_strong(context):
    return False if STRONG_FORMAT_RE.match(context) is None else True


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of chr_cpg', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()

    output_folder = args.output_folder
    if not output_folder:
        output_folder = os.path.dirname(sys.argv[0])

    files_paths = os.path.join(args.raw, GLOB_FORMAT)
    return output_folder, files_paths


def convert_context_to_int(context):
    return "".join(NUCLEOTIDE_TO_NUMBER[letter] for letter in context)


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


def get_orph_order():
    """
    Get the orph order in case someone needs it
    :return: A list where each value is the distance another CpG was searched
    :rtype: list[str]
    """
    return NUMBER_OF_ORPH_PER_INDEX


def main():
    convert_context_to_int_vc = np.vectorize(convert_context_to_int)
    is_strong_vc = np.vectorize(is_strong)
    is_weak_vc = np.vectorize(is_weak)

    output_folder, files_path = format_args()
    all_file_paths = glob.glob(files_path)

    for cpg_file in tqdm(all_file_paths):
        cpg_matrix = format_cpg_seq_file(cpg_file)
        context_col_int = convert_context_to_int_vc(cpg_matrix[:, -1])
        is_weak_col = is_weak_vc(cpg_matrix[:, -1]).reshape(cpg_matrix.shape[0], 1)
        is_strong_col = is_strong_vc(cpg_matrix[:, -1]).reshape(cpg_matrix.shape[0], 1)

        # Add new cols and convert matrix to u_int32
        cpg_matrix[:, -1] = context_col_int
        cpg_matrix = cpg_matrix.astype(np.uint32)
        cpg_matrix = np.hstack((cpg_matrix, is_weak_col))
        cpg_matrix = np.hstack((cpg_matrix, is_strong_col))

        chr_number = FILE_FORMAT_RE.findall(cpg_file)[0]
        output_path = os.path.join(output_folder, OUTPUT_FILE_FORMAT % chr_number)
        tools.save_as_compressed_pickle(output_path, cpg_matrix)


if __name__ == '__main__':
    main()
