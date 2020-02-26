import argparse
import collections
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())
from commons import data_tools

CSV_FILE = "common_cpg_in_cov_matrix_%s_chr_%s_region_%s.csv"

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_chr(\d+).dummy.pkl.zip")

ALL = 'ALL'
ALL_PT = "PT"

# The min number of pairs needed in common between cells to count the covariance
MIN_PERIODS = 10
SUBLINEAGE_MIN_PERIODS = 5


def get_files_to_work(cpg_format_files):
    if os.path.isdir(cpg_format_files):
        cpg_format_file_path = os.path.join(cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    else:
        all_cpg_format_file_paths = [cpg_format_files]

    return all_cpg_format_file_paths


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))

    args = parser.parse_args()
    return args


def main():
    args = parse_input()

    all_cpg_format_file_paths = get_files_to_work(args.cpg_format_files)

    nc_counter = collections.Counter()
    pt_counter = collections.Counter()
    for file_path in all_cpg_format_file_paths:
        patient, chro = CPG_FORMAT_FILE_RE.findall(file_path)[0]
        df = pd.read_pickle(file_path)

        nc_index = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        pt_index = [cell_id for cell_id in df.index if cell_id.startswith("PT")]

        nc_df = df.loc[nc_index, :]
        pt_df = df.loc[pt_index, :]

        nc_zero_one = np.where(~np.isnan(nc_df), 1, 0)
        pt_zero_one = np.where(~np.isnan(pt_df), 1, 0)

        nc_sum = np.sum(nc_zero_one, axis=0)
        pt_sum = np.sum(pt_zero_one, axis=0)

        nc_counter.update(nc_sum)
        pt_counter.update(pt_sum)

    data_tools.counter_to_csv(nc_counter, os.path.join(args.output_folder, "%s_nc.csv" % patient))
    data_tools.counter_to_csv(pt_counter, os.path.join(args.output_folder, "%s_pt.csv" % patient))


if __name__ == '__main__':
    main()
