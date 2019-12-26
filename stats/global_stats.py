import argparse
import collections
import glob
import os
import re
import sys

import numpy as np
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

from commons.consts import ALL_SEQ_PATH
import commons.tools as tools


def get_files_format(folder):
    """
    This is unique for us, format the files path by how the user saved them
    :param folder: The path of the folder
    :type folder: str
    :return: The formatted files locations ready to glob
    :rtype: str
    """
    return os.path.join(folder, "*", "*.pickle.zlib")


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()

    if not args.output_folder:
        output = os.path.join(os.path.dirname(__file__), "output")
    else:
        output = args.output_folder

    return args.files_folder, output


def update_total_read_count(data, reads_counter):
    unique, counts = np.unique(data[:, 1], return_counts=True)
    reads_counter.update(dict(zip(unique, counts)))
    return reads_counter


def update_ratio_count(data, ratio_counter, read_limit=-1):
    if read_limit == -1:
        unique, counts = np.unique(data[:, 2] / data[:, 1], return_counts=True)
        ratio_counter.update(dict(zip(unique, counts)))

    else:
        indexes = np.where(data[:, 1] <= read_limit)
        unique, counts = np.unique(data[indexes, 2] / data[indexes, 1], return_counts=True)
        ratio_counter.update(dict(zip(unique, counts)))

    return ratio_counter


def reads_stats(all_files, output_folder):
    total_read_counter = collections.Counter()
    total_ratio_counter = collections.Counter()
    total_ratio_for_limit_counter = collections.Counter()

    for file_path in tqdm(all_files):
        data = tools.load_compressed_pickle(file_path)
        total_read_counter = update_total_read_count(data, total_read_counter)
        total_ratio_counter = update_ratio_count(data, total_ratio_counter)
        total_ratio_for_limit_counter = update_ratio_count(data, total_ratio_for_limit_counter, 10)

    tools.counter_to_csv(total_read_counter, os.path.join(output_folder, "read_counter.csv "))
    tools.counter_to_csv(total_ratio_counter, os.path.join(output_folder, "ratio_counter.csv "))
    tools.counter_to_csv(total_ratio_for_limit_counter, os.path.join(output_folder,
                                                                     "ratio_counter_limit.csv "))


def percentage_of_positions(all_files, output_folder):
    # Get the number of positions per chr
    seq_files = glob.glob(os.path.join(ALL_SEQ_PATH, "*"))
    positions_in_chr = {}
    files_dict = {}
    all_file = open(os.path.join(output_folder, "chr_all_stats.csv"), "a")
    all_file.write("file name, chr cpg, read cpg, ratio\n")
    for seq_file in tqdm(seq_files, desc="full seq files"):
        chr_number = re.findall("\d+", seq_file)[0]
        seq_data = tools.load_compressed_pickle(seq_file)
        positions_in_chr[chr_number] = len(seq_data[:, 0])
        files_dict[chr_number] = open(os.path.join(output_folder, "chr%s_stats.csv" % chr_number), "a")
        files_dict[chr_number].write("file name, chr cpg, read cpg, ratio\n")

    # Get the number of positions per patient cell
    for read_file in tqdm(all_files, desc="patient files"):
        chr_number = re.findall("chr(\d+)|chrX|chrY", read_file)[0]
        # Skip X,Y which we don't have seq
        if chr_number not in positions_in_chr:
            continue

        num_of_expected_positions = positions_in_chr[chr_number]
        data = tools.load_compressed_pickle(read_file)
        num_of_read_positions = len(data[:, 0])
        ratio = num_of_read_positions / num_of_expected_positions
        files_dict[chr_number].write(
            "%s,%s,%s,%s\n" % (
                os.path.basename(read_file), num_of_expected_positions, num_of_read_positions, ratio))
        all_file.write(
            "%s,%s,%s,%s\n" % (
                os.path.basename(read_file), num_of_expected_positions, num_of_read_positions, ratio))

    for file_obj in files_dict.values():
        file_obj.close()


def main():
    files_folder, output_folder = format_args()
    all_files_path = get_files_format(files_folder)
    all_files = glob.glob(all_files_path)
    reads_stats(all_files, output_folder)
    percentage_of_positions(all_files, output_folder)


if __name__ == '__main__':
    tools.init_slurm(main)
