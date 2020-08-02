"""
Print several csv with information about the ratio and count number in the different files, this was used
for simple statistics
"""

import argparse
import collections
import os
import re
import sys

import numpy as np
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

from commons import files_tools, data_tools, consts


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()

    return args.files_folder, args.output_folder


def update_total_read_count(array, reads_counter):
    """
    Update a counter based on numpy data
    :param array: The np array
    :param reads_counter: The counter
    :return: The updated counter
    """
    unique, counts = np.unique(array[:, 1], return_counts=True)
    reads_counter.update(dict(zip(unique, counts)))
    return reads_counter


def update_ratio_count(data, ratio_counter, ):
    """
    Update the counter ratio
    :param data: The np array
    :param ratio_counter: The counter
    :return: The updated counter
    """
    unique, counts = np.unique(data[:, 2] / data[:, 1], return_counts=True)
    ratio_counter.update(dict(zip(unique, counts)))

    return ratio_counter


def reads_stats(all_files, output_folder):
    """
    Create csv with statistics about the read count, the ratio count
    :param all_files: All the files to read
    :param output_folder: The output folder
    """
    total_read_counter = collections.Counter()
    total_ratio_counter = collections.Counter()

    for file_path in tqdm(all_files):
        data = files_tools.load_compressed_pickle(file_path)
        total_read_counter = update_total_read_count(data, total_read_counter)
        total_ratio_counter = update_ratio_count(data, total_ratio_counter)

    data_tools.counter_to_csv(total_read_counter, os.path.join(output_folder, "read_counter.csv "))
    data_tools.counter_to_csv(total_ratio_counter, os.path.join(output_folder, "ratio_counter.csv "))


def percentage_of_positions(all_files, output_folder):
    """
    Get the number of positions per chr
    :param all_files: All the files
    :param output_folder: The output path
    """
    positions_in_chr = {}
    files_dict = {}
    all_file = open(os.path.join(output_folder, "chr_all_stats.csv"), "a")
    all_file.write("file name, chr cpg, read cpg, ratio\n")
    context_files = files_tools.get_cpg_context_map(only_locations=True)
    for chromosome_nam in context_files:
        chromosome = re.findall("\d+", chromosome_nam)[0]
        positions_in_chr[chromosome] = len(context_files[chromosome_nam][:, 0])
        files_dict[chromosome] = open(os.path.join(output_folder, "chr%s_stats.csv" % chromosome), "a")
        files_dict[chromosome].write("file name, chr cpg, read cpg, ratio\n")

    # Get the number of positions per patient cell
    for read_file in tqdm(all_files, desc="patient files"):
        chr_name = consts.CHR_FULL_NAME_RE.findall(read_file)[0]

        # Skip X,Y which we don't have seq
        if chr_name not in positions_in_chr:
            continue

        num_of_expected_positions = positions_in_chr[chr_name]
        data = files_tools.load_compressed_pickle(read_file)
        num_of_read_positions = len(data[:, 0])
        ratio = num_of_read_positions / num_of_expected_positions
        files_dict[chr_name].write(
            "%s,%s,%s,%s\n" % (
                os.path.basename(read_file), num_of_expected_positions, num_of_read_positions, ratio))
        all_file.write(
            "%s,%s,%s,%s\n" % (
                os.path.basename(read_file), num_of_expected_positions, num_of_read_positions, ratio))

    for file_obj in files_dict.values():
        file_obj.close()


def main():
    files_folder, output_folder = format_args()
    all_files = files_tools.get_files_to_work(files_folder, os.path.join("*", "*.pickle.zlib"))
    reads_stats(all_files, output_folder)
    percentage_of_positions(all_files, output_folder)


if __name__ == '__main__':
    main()
