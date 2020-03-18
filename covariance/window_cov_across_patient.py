import argparse
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
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
    parser.add_argument('--input', help='Path to bedgraph files or folder', required=False)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--window_boundaries', help='File with the window boundries', required=False)
    parser.add_argument('--data_path', help='File with data saved', required=False)

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


def window_cov_across_patients(files_paths, window_boundries):
    patients_dict = {}
    for file_path in files_paths:
        patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
        input_file = pd.read_csv(file_path, sep="\t", header=None, names=["chr", "start", "end", "cov"])
        patients_dict[patient] = input_file

    values_across_patients = []

    for i in window_boundries:
        values = []
        start, end = i

        values.extend([start, end])
        for patient in patients_dict:
            input_file = patients_dict[patient]

            try:
                value = float(input_file[input_file.start == start]["cov"])
            except TypeError:  # Will happened if we had nans in this window
                value = np.nan

            values.append(value)

        values_across_patients.append(values)

    df = pd.DataFrame(values_across_patients, columns=["start", "end"] + [i for i in patients_dict])
    return df


def create_min_bedgraph(chromosomes_dict):
    bedgraph_between_005_to_015 = open("coverage_btw_005_to_015_windows_5000_hg19.bed", "w")
    bedgraph_over_015 = open("coverage_over_015_windows_5000_hg19.bed", "w")

    for chromosome in chromosomes_dict:
        data = chromosomes_dict[chromosome]
        only_patients_data = data.iloc[:, 2:]
        values = only_patients_data.min(axis=1)

        between_005_and_015 = data[np.logical_and(values < 0.15, values >= 0.05)]
        over_015 = data[values >= 0.15]

        for index in between_005_and_015.index:
            try:
                row = between_005_and_015.loc[index]
                bedgraph_between_005_to_015.write("%s\t%s\t%s\n" % (chromosome, int(row.start), int(row.end)))
            except Exception:
                pass

        for index in over_015.index:
            row = over_015.loc[index]
            bedgraph_over_015.write("%s\t%s\t%s\n" % (chromosome, int(row.start), int(row.end)))

    bedgraph_over_015.close()
    bedgraph_between_005_to_015.close()


def create_min_histogram(chromosomes_dict):
    valid_windows = []
    nans = 0
    total = 0
    for chromosome in chromosomes_dict:
        data = chromosomes_dict[chromosome]
        only_patients_data = data.iloc[:, 2:]
        values = only_patients_data.min(axis=1)
        nans += np.sum(values.isnull())
        total += values.size
        valid_windows.extend(values[~values.isnull()].tolist())

    print("Number of windows without value is: %s/%s which is %s percent" % (nans, total, nans / total * 100))

    valid_df = pd.DataFrame(valid_windows)
    df_without_o = valid_df.where(valid_df <= valid_df.median() + 2 * valid_df.std())[0]
    _ = plt.hist(df_without_o.to_list(), bins='auto')
    plt.style.use('ggplot')

    plt.xlabel("Covariance value")
    plt.title("Histogram of min covariance value across the patients")
    plt.savefig("min_covariance_value_across_patients.png")


def main():
    args = parse_input()

    data_file = args.data_path

    if data_file:
        chromosomes_dict = files_tools.load_compressed_pickle(data_file)
        # new_dict = {}
        # for ch in chromosomes_dict:
        #     new_dict[ch] = chromosomes_dict[ch].replace(-1,np.nan)
        # files_tools.save_as_compressed_pickle("5000_window_cov_across_patients_hg19.pkl", new_dict)

    else:
        all_file_paths = get_files_to_work(args.input)
        all_files_dict = get_files_in_dict(all_file_paths)
        window_boundries = files_tools.load_compressed_pickle(args.window_boundaries)

        chromosomes_dict = {}

        for ch in tqdm(all_files_dict):
            df = window_cov_across_patients(all_files_dict[ch], window_boundries[int(ch)])
            chromosomes_dict[ch] = df

        files_tools.save_as_compressed_pickle("5000_window_cov_across_patients_hg19.pkl", chromosomes_dict)

    create_min_histogram(chromosomes_dict)
    create_min_bedgraph(chromosomes_dict)


if __name__ == '__main__':
    main()
