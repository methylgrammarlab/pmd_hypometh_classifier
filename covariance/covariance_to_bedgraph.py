import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.getcwd())
from format_files import format_sublineage_info
from commons import consts

BEDGRAPH_LINE_FORMAT = "{chr_name}\t{start}\t{end}\t{number}\n"
BEDGRAPH_OUTPUT_FILE_FORMAT = "average_covariance_between_cpg_%s_chr_%s_region_%s.bedgraph"

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
    parser.add_argument('--run_on_sublineage', help='Should the code create one bedgraph or one per '
                                                    'sublineage', required=False, default=False, type=bool)
    parser.add_argument('--window_size', help='The window size, default is 500', default=500,
                        required=False, type=int)

    args = parser.parse_args()
    return args


def create_region_bedgraph(df_path, sublineage_cells, sublineage_name, output_path, window_size, chromosome):
    """
    Creates a bedgraph of the mean covariance of CpG's, using data only from one region at a time,
    adding the NC. The covariance is calculated in windows.
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
    elif sublineage_cells == ALL_PT:
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

    with open(output_path, "w") as output_file:
        for i in tqdm(range(0, num_of_cpg, window_size)):
            window_indexes = region_df.columns[i:min(i + window_size, num_of_cpg)]  # Get the indexes
            covariance_matrix = region_df.loc[:, window_indexes].cov(min_periods=min_periods)  # Create cov
            np.fill_diagonal(covariance_matrix.values, np.nan)  # Remove the diagonal
            average_covariance = covariance_matrix.mean()  # Create the avg

            # Write the data
            for cpg in window_indexes:
                number = average_covariance[cpg]
                if not np.isnan(number):
                    line = BEDGRAPH_LINE_FORMAT.format(chr_name=chromosome, start=cpg, end=cpg + 1,
                                                       number=number)
                    output_file.write(line)


def create_bedgraphs(filename, output_path, window_size=500, run_on_sublineage=False):
    """
    Goes over all the regions and calls the function to create a bedgraph for each one.
    :param filename: The filename of the file with the parsed scWGBS data
    :param window_size: The size of the window to work with
    :param output_path: The output folder
    :param run_on_sublineage: Should we run the information on all the sublinage
    """
    sublineage_info = format_sublineage_info.get_sublineage_info(consts.SUBLINEAGE_FILE_LOCAL_DROR)
    patient, chromosome = CPG_FORMAT_FILE_RE.findall(filename)[0]

    if run_on_sublineage:
        patient_info = sublineage_info[patient]
        for sublineage_name in patient_info:
            output_filename = os.path.join(output_path,
                                           BEDGRAPH_OUTPUT_FILE_FORMAT % (patient, chromosome, 'NCand%s' %
                                                                          sublineage_name))
            create_region_bedgraph(df_path=filename, sublineage_cells=patient_info[sublineage_name],
                                   sublineage_name=sublineage_name, output_path=output_filename,
                                   window_size=window_size, chromosome=chromosome)
    else:
        output_filename = os.path.join(output_path, BEDGRAPH_OUTPUT_FILE_FORMAT %
                                       (patient, chromosome, 'NCandPT'))

        create_region_bedgraph(df_path=filename, sublineage_cells=[], sublineage_name=ALL_PT,
                               output_path=output_filename, window_size=window_size, chromosome=chromosome)


def main():
    args = parse_input()

    all_cpg_format_file_paths = get_files_to_work(args.cpg_format_files)

    for file_path in all_cpg_format_file_paths:
        create_bedgraphs(filename=file_path, output_path=args.output_folder, window_size=args.window_size,
                         run_on_sublineage=args.run_on_sublineage)


if __name__ == '__main__':
    main()