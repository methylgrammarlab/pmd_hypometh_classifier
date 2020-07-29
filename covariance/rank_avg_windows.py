"""
Take different bedgraph files which represent windows which were normalized ( so one value per window) and
rank them across patient to see if they match
"""

import argparse
import glob
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
    parser.add_argument('--input', help='Path to bedgraph files or folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--window_boundaries', help='File with the window boundries', required=False)

    args = parser.parse_args()
    return args


def get_bedgraph_files(input_path):
    """
    Get a list of all the bedgraph files to work on
    :param input_path: A dir path or file path of bedgraph
    :return: A list of all the paths to work on
    """
    if os.path.isdir(input_path):
        file_path = os.path.join(input_path, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [input_path]

    return all_file_paths


def rank_covariance_across_patients(files_paths, window_boundaries):
    """
    Get the covariance of a window and rank the different windows
    We base this on the fact that this is the normed bedgraph, so each window has only one value or nan
    :param files_paths: A list of paths to work on
    :param window_boundaries: The boundaries for the chromosome
    :return:
    """
    patients_dict = {}
    for file_path in files_paths:
        patient, chromosome = consts.PATIENT_CHR_NAME_RE.findall(file_path)[0]
        input_file = pd.read_csv(file_path, sep="\t", header=None, names=["chr", "start", "end", "cov"])
        values = []

        for i in window_boundaries:
            try:
                value = float(input_file[input_file.start == i[0]]["cov"])
            except TypeError:  # Will happened if we have all nans nans in this window
                value = -1

            values.append(value)

        patients_dict[patient] = pd.Series(values).rank(method="min")

    # Chose the first patient to be the baseline
    baseline_patient_name = list(patients_dict.keys())[0]
    baselines_p = np.copy(patients_dict[baseline_patient_name])

    for patient in patients_dict:
        patients_dict[patient] = [x for _, x in sorted(zip(baselines_p, patients_dict[patient]))]

    return pd.DataFrame(patients_dict)



def main():
    args = parse_input()

    all_file_paths = get_bedgraph_files(args.input)
    all_files_dict = files_tools.convert_paths_list_to_chromosome_based_dict(all_file_paths)
    window_boundaries = files_tools.load_compressed_pickle(args.window_boundaries)

    for chromosome in tqdm(all_files_dict):
        df = rank_covariance_across_patients(all_files_dict[chromosome], window_boundaries[int(chromosome)])
        df.to_csv(os.path.join(args.output_folder, "covariance_rank_ch%s.csv" % chromosome))


if __name__ == '__main__':
    main()
