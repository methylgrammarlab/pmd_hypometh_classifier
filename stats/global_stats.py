import argparse
import collections
import glob
import os

import numpy as np
import pandas as pd

import commons.tools as tools


def get_files_format(folder):
    """
    This is unique for us, format the files path by how the user saved them
    :param folder: The path of the folder
    :type folder: str
    :return: The formatted files locations ready to glob
    :rtype: str
    """
    return os.path.join(folder, "CRC09", "*.pickle.zlib")


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()

    if not args.output_folder:
        output = os.path.dirname(__file__)
    else:
        output = args.output_folder

    return args.files_folder, output


def update_total_read_count(data, reads_counter):
    unique, counts = np.unique(data[:, 1], return_counts=True)
    reads_counter.update(dict(zip(unique, counts)))
    return reads_counter


def update_ratio_count(data, ratio_counter):
    unique, counts = np.unique(data[:, 2] / data[:, 1], return_counts=True)
    ratio_counter.update(dict(zip(unique, counts)))
    return ratio_counter


def main():
    files_folder, output_folder = format_args()
    all_files_path = get_files_format(files_folder)
    all_files = glob.glob(all_files_path)
    total_read_counter = collections.Counter()
    total_ratio_counter = collections.Counter()
    for file_path in all_files:
        data = tools.load_compressed_pickle(file_path)
        total_read_counter = update_total_read_count(data, total_read_counter)
        total_ratio_counter = update_ratio_count(data, total_ratio_counter)

    read_counter_df = pd.DataFrame.from_dict(total_read_counter, orient='index').reset_index()
    read_counter_df.columns = ["read_number", "counter"]
    ratio_counter_df = pd.DataFrame.from_dict(total_ratio_counter, orient='index').reset_index()
    ratio_counter_df.columns = ["ratio", "counter"]
    read_counter_df.to_csv(os.path.join(output_folder, "read_counter.csv "))
    ratio_counter_df.to_csv(os.path.join(output_folder, "ratio_counter.csv "))


if __name__ == '__main__':
    main()
