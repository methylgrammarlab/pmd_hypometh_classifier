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
    parser.add_argument('--prediction_file', help='Path to prediction file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def get_patient_df_dict(all_file_paths):
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


def get_sublineage(cell, sublineage_info, patient):
    """
    :param cell: The cell name
    :param sublineage_info: Dictionary of cell names -> sublineage
    :param patient: Patient name
    :return: The cells sublineage if there is, undefined otherwise
    """
    try:
        return sublineage_info[patient + '_' + cell]
    except:
        if cell.startswith("NC"):
            return "NC"
        return "undefined"


def add_meta(mean_series, coverage_series, sublineage, patient, hypo):
    mean_df = pd.DataFrame(index=list(mean_series.index))
    mean_df.loc[:, "mean"] = mean_series
    mean_df.loc[:, "coverage"] = coverage_series
    mean_df.loc[:, "hypo"] = hypo
    mean_df.loc[:, "patient"] = patient
    mean_df.loc[:, "sublineage"] = mean_df.index
    mean_df.loc[:, "sublineage"] = mean_df.loc[:, 'sublineage'].apply(get_sublineage, args=(sublineage, patient))
    return mean_df


def calculate_all_meth_mean_df(meth_df_dict, loc_df, sublineage):
    all_patients = []
    for patient in tqdm(meth_df_dict):
        patient_cell_names = get_patient_cells(meth_df_dict[patient])
        prone_mean, resistant_mean, prone_coverage, resistant_coverage = calculate_patient_meth_mean(
            meth_df_dict[patient], loc_df, patient_cell_names)
        prone_mean = add_meta(prone_mean, prone_coverage, sublineage, patient, "prone")
        resistant_mean = add_meta(resistant_mean, resistant_coverage, sublineage, patient, "resistant")
        all_patients.append(prone_mean)
        all_patients.append(resistant_mean)
    return pd.concat(all_patients)


def get_patient_cells(patient_meth_df_dict):
    """
    Find cell names
    :return: Set of all the patients cell names
    """
    patient_cell_names = set()
    for chromosome in patient_meth_df_dict:
        patient_chromosome_df = patient_meth_df_dict[chromosome]
        patient_cell_names |= set(patient_chromosome_df.index)
    return patient_cell_names


def calculate_patient_meth_mean(patient_meth_df_dict, loc_df, patient_cell_names):
    patient_prone = []
    patient_resistant = []
    for chromosome in patient_meth_df_dict:
        chromosome_prone = subset_df(patient_meth_df_dict, loc_df, chromosome, "prone")
        chromosome_resistant = subset_df(patient_meth_df_dict, loc_df, chromosome, "resistant")
        patient_prone.append(chromosome_prone)
        patient_resistant.append(chromosome_resistant)
    prone_mean, prone_coverage = calculate_df_list_mean(patient_prone, patient_cell_names)
    resistant_mean, resistant_coverage = calculate_df_list_mean(patient_resistant, patient_cell_names)
    return prone_mean, resistant_mean, prone_coverage, resistant_coverage


def calculate_df_list_mean(patient_df_list, patient_cell_names):
    sum_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    num_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    tot_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    for i in range(len(patient_df_list)):
        chromosome_df = patient_df_list[i]
        sum_methylation_df.loc[chromosome_df.index, i] = np.sum(chromosome_df, axis=1)
        num_methylation_df.loc[chromosome_df.index, i] = np.sum(np.isfinite(chromosome_df), axis=1)
        tot_methylation_df.loc[chromosome_df.index, i] = chromosome_df.shape[1]
    mean_df = np.sum(sum_methylation_df, axis=1) / np.sum(num_methylation_df, axis=1)
    coverage_df = np.sum(num_methylation_df, axis=1) / np.sum(tot_methylation_df, axis=1)
    return mean_df, coverage_df


def subset_df(patient_meth_df_dict, loc_df, chromosome, hypo):
    chromosome_meth_df = patient_meth_df_dict[chromosome]
    chromosome_loc = loc_df.loc[loc_df.loc[:, "chromosome"] == chromosome, :]
    if hypo == "prone":
        loc_list = chromosome_loc.loc[chromosome_loc.loc[:, "pred"] > 0.5, "location"]
    else:
        loc_list = chromosome_loc.loc[chromosome_loc.loc[:, "pred"] < 0.5, "location"]
    return chromosome_meth_df.loc[:, loc_list]


def save_to_file(output_folder, all_patients_mean):
    """
    Saves both the mean DF and the median df to a pickle and a csv
    """
    print("saving...")
    path = os.path.join(output_folder, "prone_resistant_avg_meth_coverage.csv")
    all_patients_mean.to_csv(path)


def main():
    args = parse_input()

    sublineage = files_tools.load_compressed_pickle(consts.CONVERT_SUBLINEAGE_LIOR_AQUARIUM)
    predictions = pd.read_csv(args.prediction_file, index_col=0)

    all_files = files_tools.get_files_to_work(args.methylation_folder, pattern=os.path.join("*", "*.pkl.zip"))
    patients_dict = get_patient_df_dict(all_files)
    all_patients_mean = calculate_all_meth_mean_df(patients_dict, predictions, sublineage)
    save_to_file(args.output_folder, all_patients_mean)


if __name__ == '__main__':
    t0 = time.time()
    main()
    t1 = time.time()
    print("total time:", t1 - t0)
