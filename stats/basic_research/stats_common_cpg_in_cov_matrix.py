"""
Was used to find how many pairs in the covariance matrix make sense for minimum number to make valid cpg
"""

import argparse
import collections
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
from format_files import format_sublineage_info
from commons import data_tools

CSV_FILE = "common_cpg_in_cov_matrix_%s_chr_%s_region_%s.csv"

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_chr(\d+).dummy.pkl.zip")

ALL = 'ALL'
ALL_PT = "PT"

# The min number of pairs needed in common between cells to count the covariance
MIN_PERIODS = 10
SUBLINEAGE_MIN_PERIODS = 5

PERIOD_MIN = 5
PERIOD_MAX = 15


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
    parser.add_argument('--run_on_sublineage', required=False, default=False, action='store_true',
                        help='Should the code create one bedgraph or one per sublineage')
    parser.add_argument('--window_size', help='The window size, default is 500', default=5000,
                        required=False, type=int)

    args = parser.parse_args()
    return args


def create_histogram(df_path, sublineage_cells, sublineage_name, output_path, window_size, chromosome):
    """
    Create histogram for the number of CpG each CpG has enough common things to calculate covariance
    :param df_path: The filename of the file with the parsed scWGBS data
    :param sublineage_cells: The cells in the sublineage
    :param sublineage_name: The name of the sublineage_name
    :param window_size: The window_size
    :param output_path: The output folder
    :param chromosome: The name of the chromosome we are working on
    """
    df = pd.read_pickle(df_path)
    min_periods = MIN_PERIODS

    # All cells
    if sublineage_name == ALL:
        region_df = df

    # NC + PT
    elif sublineage_name == ALL_PT:
        region_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        region_cell_ids.extend(cell_id for cell_id in df.index if cell_id.startswith("PT"))
        region_df = df.loc[region_cell_ids, :]

    # NC + sublineage
    else:
        region_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        for sublineage_cell in sublineage_cells:
            region_cell_ids.extend(cell_id for cell_id in df.index if cell_id.startswith(sublineage_cell))

        region_df = df.loc[region_cell_ids, :]
        min_periods = SUBLINEAGE_MIN_PERIODS

    num_of_cpg = region_df.shape[1]

    counter = collections.Counter()
    for i in tqdm(range(0, num_of_cpg, window_size)):
        window_indexes = region_df.columns[i:min(i + window_size, num_of_cpg)]  # Get the indexes
        covariance_matrix = region_df.loc[:, window_indexes].cov(min_periods=min_periods)  # Create cov
        np.fill_diagonal(covariance_matrix.values, np.nan)  # Remove the diagonal
        counter.update(np.sum(~covariance_matrix.isnull()))

    data_tools.counter_to_csv(counter, output_path)


def create_cancer_histogram_multiple_min_periods(df_path, output_path, window_size, patient):
    df = pd.read_pickle(df_path)
    counters = [collections.Counter() for i in range(PERIOD_MIN, PERIOD_MAX)]
    col_list = ["min_pairs_of_%s" %i for i in range(PERIOD_MIN, PERIOD_MAX)]

    region_cell_ids = [cell_id for cell_id in df.index if not cell_id.startswith('NC')]

    region_df = df.loc[region_cell_ids, :]
    num_of_cpg = region_df.shape[1]

    for i in tqdm(range(0, num_of_cpg, window_size)):
        window_indexes = region_df.columns[i:min(i + window_size, num_of_cpg)]  # Get the indexes

        for min_periods in range(PERIOD_MIN, PERIOD_MAX):
            covariance_matrix = region_df.loc[:, window_indexes].cov(min_periods=min_periods)  # Create cov
            np.fill_diagonal(covariance_matrix.values, np.nan)  # Remove the diagonal
            counters[min_periods - PERIOD_MIN].update(np.sum(~covariance_matrix.isnull()))

    lines = []
    for i in range(window_size+1):
        line = [i] + [c.get(i) for c in counters]
        lines.append(line)

    new_df = pd.DataFrame.from_records(lines, columns=["number_of_pairs"] + col_list)
    new_df.to_csv(os.path.join(output_path, "%s_compare_min_pairs.csv" %patient))


def create_histograms(filename, output_path, window_size=500, run_on_sublineage=False):
    """
    Goes over all the regions and calls the function to create a bedgraph for each one.
    :param filename: The filename of the file with the parsed scWGBS data
    :param window_size: The size of the window to work with
    :param output_path: The output folder
    :param run_on_sublineage: Should we run the information on all the sublinage
    """
    sublineage_info = format_sublineage_info.get_sublineage_info()
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(filename)[0]

    if run_on_sublineage:
        patient_info = sublineage_info[patient]
        for sublineage_name in patient_info:
            output_filename = os.path.join(output_path, CSV_FILE % (patient, chromosome, 'NCand%s' %
                                                                    sublineage_name))
            create_histogram(df_path=filename, sublineage_cells=patient_info[sublineage_name],
                             sublineage_name=sublineage_name, output_path=output_filename,
                             window_size=window_size, chromosome=chromosome)
    else:
        output_filename = os.path.join(output_path, CSV_FILE % (patient, chromosome, 'NCandPT'))

        create_histogram(df_path=filename, sublineage_cells=[], sublineage_name=ALL_PT,
                         output_path=output_filename, window_size=window_size, chromosome=chromosome)


def main():
    args = parse_input()

    all_cpg_format_file_paths = get_files_to_work(args.cpg_format_files)

    for file_path in all_cpg_format_file_paths:
        # create_histograms(filename=file_path, output_path=args.output_folder, window_size=args.window_size,
        #                   run_on_sublineage=args.run_on_sublineage)
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]
        create_cancer_histogram_multiple_min_periods(file_path, output_path=args.output_folder, window_size=args.window_size, patient=patient)


if __name__ == '__main__':
    main()
