import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr_(\d+).*")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedgraph_files', help='Path to bedgraph files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_files_to_work(files):
    if os.path.isdir(files):
        file_path = os.path.join(files, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [files]

    return all_file_paths


def get_files_in_dict(all_file_paths):
    d = {}
    for file_path in all_file_paths:
        try:
            patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
        except:
            continue

        if chromosome not in d:
            d[chromosome] = []

        d[chromosome].append(file_path)

    return d


def main():
    args = parse_input()

    all_file_paths = get_files_to_work(args.bedgraph_files)
    all_files_dict = get_files_in_dict(all_file_paths)

    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    temp_dict = {}

    for chromosome in all_files_dict:
        temp_dict[chromosome] = []
        windows_data = global_windows_data[chromosome]
        for file_path in all_files_dict[chromosome]:
            patient, _ = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
            covariance_pmd_df = handle_pmds.get_covariance_pmd_df(file_path, chromosome)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.index >= start) & (covariance_pmd_df.index <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            temp_dict[chromosome].append((patient, covariance_pmd_df.loc[prev_mask, :]))

    # Until here we have a list of tuples of pd with data after filtering windows with high cov and pmds


if __name__ == '__main__':
    main()
