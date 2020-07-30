"""
Take scWGBS dataframes with all the cells and convert to a bedgraph with the covaraince is the value of the CpG
This code support taking several cells based on sublineage or region name(e.g. PT, LN, NC, PT&LN, etc...)
This code also allow the user to define of the covariance should be caluclated:
- Min numbers of samples for a valid cpg
- Min numbers of pairs between CpG in a window to be a valid CpG
- Window size
- Advancement index which allow for overlapping windows
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
from format_files import format_sublineage_info
from commons import consts, files_tools, utils

BEDGRAPH_LINE_FORMAT = "{chr_name}\t{start}\t{end}\t{number}\n"
BEDGRAPH_OUTPUT_FILE_FORMAT = "average_covariance_between_cpg_%s_%s_region_%s.bedgraph"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--extract_per_sublineage', required=False, default=False, type=bool,
                        help='Should the code create one bedgraph for all sublineage or one per sublineage')
    parser.add_argument('--region_name', default=utils.ALL, required=False, type=str,
                        help='Region name for bedgraph cells, is invalid if using extract_per_sublineage. '
                             'Default is ALL')
    parser.add_argument('--window_size', help='The window size, default is 5000', default=5000,
                        required=False, type=int)
    parser.add_argument('--index_advance_value', default=5000, required=False, type=int,
                        help='How many indexes to advance per run, if we want to do windows which contains '
                             'same CpG we need to change it, default is 5000')
    parser.add_argument('--min_pairs', help='Minimun number of pairs to be a valid CpG, default is 750',
                        default=750, required=False, type=int)
    parser.add_argument('--min_samples', help='Minimun number of cells per cpg to be a valid CpG, default is '
                                              '10', default=10, required=False, type=int)

    args = parser.parse_args()
    return args


def create_region_bedgraph(df, output_path, chromosome, window_size=5000, index_advance_value=5000,
                           min_numbers_of_pairs=450, min_samples_for_cpg=10):
    """
    Creates a bedgraph of the mean covariance of CpG's, using data only from one region at a time,
    adding the NC. The covariance is calculated in windows.
    :param df: The df to work on
    :param output_path: The output path for the bedgrph
    :param chromosome: The name of the chromosome we are working on
    :param window_size: The size of each window in the bedgraph for calculation
    :param index_advance_value: Amount of CpG index to advace each run, this is used if we want to create
    overlapping windows
    :param min_samples_for_cpg: Minimum samples needed for a cpg to be valid
    :param min_numbers_of_pairs: Min pairs for CpG in a window to be valid
    """
    num_of_cpg = df.shape[1]

    cpg_wrote = 0
    nans_removed = 0
    not_enough_pairs_removed = 0

    output_file = open(output_path, "w")
    for i in range(0, num_of_cpg, index_advance_value):
        window_indexes = df.columns[i:min(i + window_size, num_of_cpg)]  # Get the indexes
        covariance_matrix = df.loc[:, window_indexes].cov(min_periods=min_samples_for_cpg)  # Create cov
        np.fill_diagonal(covariance_matrix.values, np.nan)  # Remove the diagonal
        average_covariance = covariance_matrix.mean()  # Create the avg
        pairs_matrix = covariance_matrix.notnull().sum()

        # Write the data
        for cpg in window_indexes:
            number = average_covariance[cpg]
            if not np.isnan(number) and pairs_matrix[cpg] > min_numbers_of_pairs:
                line = BEDGRAPH_LINE_FORMAT.format(chr_name=chromosome, start=cpg, end=cpg + 1, number=number)
                output_file.write(line)
                cpg_wrote += 1

            else:
                if np.isnan(number):
                    nans_removed += 1
                else:
                    not_enough_pairs_removed += 1

    output_file.close()
    total = cpg_wrote + nans_removed + not_enough_pairs_removed
    print("Total cpg in chr: %s: %s, written: %s (%s), removed because of nan: %s (%s), removed because of "
          "little pairs: %s (%s)" % (chromosome, total, cpg_wrote, cpg_wrote / total * 100,
                                     nans_removed, nans_removed / total * 100,
                                     not_enough_pairs_removed, not_enough_pairs_removed / total * 100))


def create_bedgraphs(scwgbs_file_path, output_folder, extract_per_sublineage=False,
                     region_name=utils.ALL, window_size=5000, index_advance_value=5000,
                     min_numbers_of_pairs=450, min_samples_for_cpg=10):
    """
    Goes over all the sublineages and calls the function to create a bedgraph for each one.
    :param scwgbs_file_path: The filename of the file with the parsed scWGBS data
    :param window_size: The size of the window to work with
    :param output_folder: The output folder
    :param extract_per_sublineage: Should we run the information on all the sublinage
    :param region_name: Name of the region to export cells of
    :param window_size: The size of each window in the bedgraph for calculation
    :param index_advance_value: Amount of CpG index to advace each run, this is used if we want to create
    overlapping windows
    :param min_samples_for_cpg: Minimum samples needed for a cpg to be valid
    :param min_numbers_of_pairs: Min pairs for CpG in a window to be valid
    """
    patient, chromosome = consts.PATIENT_CHR_NAME_RE.findall(scwgbs_file_path)[0]
    sublineage_data = format_sublineage_info.get_sublineage_info() if extract_per_sublineage else {}

    df = pd.read_pickle(scwgbs_file_path)

    if extract_per_sublineage:
        patient_info = sublineage_data[patient]
        for sublineage_name in patient_info:
            output_filename = os.path.join(output_folder, BEDGRAPH_OUTPUT_FILE_FORMAT % (patient, chromosome,
                                                                                         sublineage_name))

            region_df = utils.filter_df_based_on_cells(df, patient_info[sublineage_name])
            create_region_bedgraph(df=region_df, chromosome=chromosome, output_path=output_filename,
                                   window_size=window_size, index_advance_value=index_advance_value,
                                   min_numbers_of_pairs=min_numbers_of_pairs,
                                   min_samples_for_cpg=min_samples_for_cpg)

    else:
        output_filename = os.path.join(output_folder, BEDGRAPH_OUTPUT_FILE_FORMAT % (patient, chromosome,
                                                                                     region_name))
        region_df = utils.filter_df_based_on_region_name(df, region_name)
        create_region_bedgraph(df=region_df, chromosome=chromosome, output_path=output_filename,
                               window_size=window_size, index_advance_value=index_advance_value,
                               min_numbers_of_pairs=min_numbers_of_pairs,
                               min_samples_for_cpg=min_samples_for_cpg)


def main():
    args = parse_input()

    all_cpg_format_file_paths = files_tools.get_files_to_work(args.cpg_format_files,
                                                              pattern=consts.SCWGBS_FILE_FORMAT % ('*', '*'))

    for file_path in all_cpg_format_file_paths:
        create_bedgraphs(scwgbs_file_path=file_path, output_folder=args.output_folder,
                         extract_per_sublineage=args.extract_per_sublineage, region_name=args.region_name,
                         window_size=args.window_size, index_advance_value=args.index_advance_value,
                         min_samples_for_cpg=args.min_samples, min_numbers_of_pairs=args.min_pairs
                         )


if __name__ == '__main__':
    main()
