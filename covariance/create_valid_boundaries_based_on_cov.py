"""
First:
Create a dataframe with the average covariance of patient per window
start        end     CRC11     CRC01     CRC13
 0     20304036   20958376  0.053869  0.056528  0.086328

The output is a dict with the chromosome name and a df with this data

Second:
Create a dict of chromosome: list of boundaries valid based on covarinace
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools, consts

BEDGRPH_FILE_FORMAT = os.path.join("*", "norm", "*.bedgraph")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Path to bedgraph files or folder', required=False)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--window_boundaries', help='File with the window boundries', required=False)
    parser.add_argument('--chromosome_df_file', help='File with chromosome_df_file saved', required=False)
    parser.add_argument('--min_cov', help='minimum covariance to be a valid window', required=False,
                        default=0.15)
    parser.add_argument('--max_cov', help='max covariance to be a valid window', required=False,
                        default=np.inf)

    args = parser.parse_args()
    return args


def get_chr_avg_cov_df(bedgraphs_paths, window_boundaries):
    """
    Create a df with the mean covariance value across different patient using bedgraph
    :param bedgraphs_paths: Path for the bedgraph files to combine
    :param window_boundaries: The window boundaries list
    :type window_boundaries: list[tuple]
    :return: The df of the average covariance per window and patient
    """
    patients_dict = {}
    for file_path in bedgraphs_paths:
        patient, chromosome = consts.PATIENT_CHR_NAME_RE.findall(file_path)[0]
        input_file = files_tools.load_bedgraph(file_path)
        patients_dict[patient] = input_file

    values_across_patients = []

    for boundary in window_boundaries:
        patient_windows_values = []
        start, end = boundary

        for patient in patients_dict:
            input_file = patients_dict[patient]

            try:
                window_values = input_file[np.logical_and(input_file.start > start, input_file.end < end)][
                    "coverage"]
                average_value = float(window_values.mean())
            except TypeError:  # Will happened if we only have nans in this window
                average_value = np.nan

            patient_windows_values.append([start, end, average_value])

        values_across_patients.append(patient_windows_values)

    df = pd.DataFrame(values_across_patients, columns=["start", "end"] + [i for i in patients_dict])
    return df


def create_valid_windows_dict(chromosomes_dict, min_cov=0.15, max_cov=np.inf):
    """
    Create a dictionary with chromosome: (start,end) for each window that has a covarinace value higher than
     min_cov
    :param chromosomes_dict: The dict of chromosome and the df matching it
    :param min_cov: The minimum covariance value to consider it a valid window
    :param max_cov: The max covariance value to consider it a valid window
    :return: A dict with chromosome and list of valid boundaries
    """
    windows_dict = {}

    for chromosome in chromosomes_dict:
        data = chromosomes_dict[chromosome]
        data = data.replace(np.nan, -1)
        only_patients_data = data.iloc[:, 2:]
        values = only_patients_data.min(axis=1)  # We only take places where the min is larger than min_cov

        valid_rows = data[np.logical_and(values >= min_cov, values <= max_cov)]
        windows_dict[chromosome] = list(valid_rows[["start", "end"]].to_records(index=False))

    return windows_dict


def main():
    args = parse_input()
    output_folder = args.output_folder
    chromosome_df_file = args.chromosome_df_file

    # If chromosome_df_file exists use it if not create it
    if chromosome_df_file:
        chromosome_df_dict = files_tools.load_compressed_pickle(chromosome_df_file)

    else:
        all_file_paths = files_tools.get_files_to_work(args.input, BEDGRPH_FILE_FORMAT)
        all_files_dict = files_tools.convert_paths_list_to_chromosome_based_dict(all_file_paths)
        window_boundaries = files_tools.load_compressed_pickle(args.window_boundaries)

        chromosome_df_dict = {}

        for chromosome in tqdm(all_files_dict):
            df = get_chr_avg_cov_df(all_files_dict[chromosome], window_boundaries[int(chromosome)])
            chromosome_df_dict[chromosome] = df

        files_tools.save_as_compressed_pickle(os.path.join(output_folder, "chromosome_df_file_hg19.pkl"),
                                              chromosome_df_dict)

    windows_dict = create_valid_windows_dict(chromosome_df_dict, args.min_cov)
    output_file = os.path.join(output_folder, "windows_with_cov_over_%s.pkl" % args.min_cov)
    files_tools.save_as_compressed_pickle(output_file, windows_dict)


if __name__ == '__main__':
    main()
