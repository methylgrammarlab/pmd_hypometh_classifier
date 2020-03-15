import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools

CSV_FILE = "common_cpg_in_cov_matrix_%s_chr_%s_region_%s.csv"

BEDGRPH_FILE_FORMAT = os.path.join("*", "norm", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")


def get_files_to_work(files):
    if os.path.isdir(files):
        file_path = os.path.join(files, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [files]

    return all_file_paths


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Path to bedgraph files or folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--window_boundaries', help='File with the window boundries', required=False)

    args = parser.parse_args()
    return args


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


def rank_covariance_across_patients(files_paths, window_boundries):
    patients_dict = {}
    for file_path in files_paths:
        patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
        input_file = pd.read_csv(file_path, sep="\t", header=None, names=["chr", "start", "end", "cov"])
        values = []
        for i in window_boundries:
            try:
                value = float(input_file[input_file.start == i[0]]["cov"])
            except TypeError:  # Will happened if we had nans in this window
                value = -1
            values.append(value)

        patients_dict[patient] = pd.Series(values).rank(method="min")

    baseline_patient_name = list(patients_dict.keys())[0]
    baselines_p = np.copy(patients_dict[baseline_patient_name])

    for patient in patients_dict:
        patients_dict[patient] = [x for _, x in sorted(zip(baselines_p, patients_dict[patient]))]

    df = pd.DataFrame(patients_dict)
    df.to_csv("ch%s.csv" % chromosome)


def main():
    args = parse_input()

    all_file_paths = get_files_to_work(args.input)
    all_files_dict = get_files_in_dict(all_file_paths)
    window_boundries = files_tools.load_compressed_pickle(args.window_boundaries)

    for ch in tqdm(all_files_dict):
        rank_covariance_across_patients(all_files_dict[ch], window_boundries[int(ch)])


if __name__ == '__main__':
    main()
