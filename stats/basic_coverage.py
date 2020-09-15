import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm
import time

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools
import commons.consts as consts


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--data_type', help='sc or bulk', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def sc_get_patient_df_dict(all_file_paths):
    """
    Loads all the methylation DF and returns them as a dictionary
    :param all_file_paths: List of all the methylation files
    :return: A dictionary were the keys are the patient and the values a dictionary of chromosome and the methylation df
    """
    patient_df_dict = {}
    for file_path in tqdm(all_file_paths):
        patient, chromosome = consts.DATA_FILE_SCWGBS_RE.findall(file_path)[0]
        if patient not in patient_df_dict:
            patient_df_dict[patient] = {}

        df = pd.read_pickle(file_path)

        patient_df_dict[patient][chromosome] = df

    return patient_df_dict


def bulk_get_df_dict(all_file_paths):
    """
    Loads all the methylation DF and returns them as a dictionary
    :param all_file_paths: List of all the methylation files
    :return: A dictionary were the keys are the patient and the values a dictionary of chromosome and the methylation df
    """
    df_dict = {}
    for file_path in tqdm(all_file_paths):
        chromosome = consts.DATA_FILE_BULK_RE.findall(file_path)[0]
        df = pd.read_pickle(file_path)

        df_dict[chromosome] = df

    return df_dict


def calculate_patient_chr(df, patient, chromosome):
    cell_coverage = np.sum(np.isfinite(df), axis=1) / df.shape[1]
    loc_coverage = np.sum(np.isfinite(df), axis=0) / df.shape[0]
    return patient, chromosome, \
           loc_coverage.shape[0], loc_coverage.mean(), loc_coverage.median(), \
           cell_coverage.shape[0], cell_coverage.mean(), cell_coverage.median()


def sc_calculate_all_meth_mean_df(meth_df_dict):
    index = [patient + "chr%d" % chromosome for patient in meth_df_dict for chromosome in range(1, 23)]
    columns = ["patient", "chromosome", "loc num", "loc mean", "loc median", "cell num", "cell mean", "cell median"]
    all_patients = pd.DataFrame(index=index, columns=columns)
    for patient in tqdm(meth_df_dict):
        for chromosome in meth_df_dict[patient]:
            df = meth_df_dict[patient][chromosome]
            row = patient + chromosome
            all_patients.loc[row, :] = calculate_patient_chr(df, patient, chromosome)
    return all_patients


def bulk_calculate_all_meth_mean_df(meth_df_dict):
    index = ["chr%d" % chromosome for chromosome in range(1, 23)]
    columns = ["patient", "chromosome", "loc num", "loc mean", "loc median", "cell num", "cell mean", "cell median"]
    all_patients = pd.DataFrame(index=index, columns=columns)
    for chromosome in meth_df_dict:
        df = meth_df_dict[chromosome]
        all_patients.loc[chromosome, :] = calculate_patient_chr(df, 'bulk', chromosome)
    return all_patients


def save_to_file(output_folder, all_patients_mean):
    """
    Saves both the mean DF and the median df to a pickle and a csv
    """
    print("saving...")
    path = os.path.join(output_folder, "basic_coverage_stats.csv")
    all_patients_mean.to_csv(path)


def sc_main(methylation_folder):
    all_files = files_tools.get_files_to_work(methylation_folder, pattern=os.path.join("*", "*.pkl.zip"))
    patients_dict = sc_get_patient_df_dict(all_files)
    all_patients_mean = sc_calculate_all_meth_mean_df(patients_dict)
    return all_patients_mean


def bulk_main(methylation_folder):
    all_files = files_tools.get_files_to_work(methylation_folder, pattern="*.pkl.zip")
    patients_dict = bulk_get_df_dict(all_files)
    all_patients_mean = bulk_calculate_all_meth_mean_df(patients_dict)
    return all_patients_mean


def main():
    args = parse_input()

    if args.data_type == "sc":
        all_patients_mean = sc_main(args.methylation_folder)
    elif args.data_type == "bulk":
        all_patients_mean = bulk_main(args.methylation_folder)

    save_to_file(args.output_folder, all_patients_mean)


if __name__ == '__main__':
    t0 = time.time()
    main()
    t1 = time.time()
    print("total time:", t1 - t0)
