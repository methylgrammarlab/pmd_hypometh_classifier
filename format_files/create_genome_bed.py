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
    parser.add_argument('--QC_path', help='Path to QC file', required=True)
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


def get_QC(QC_path):
    qc_df = pd.read_csv(QC_path)
    qc_dict = {patient: qc_df.loc[qc_df.patient == patient, "X"] for patient in qc_df.patient.unique()}
    return qc_dict


def combine_df(df_list):
    merged = pd.concat(df_list)
    return merged


def save_to_file(patient, all_patient_chr, out_path):
    filename = "filtered_cpg_ratios_%s_hg19.bedgraph.gz"
    path = os.path.join(out_path, filename % patient)
    all_patient_chr.to_csv(path, compression='gzip', sep='\t', index=False)


def create_bed(patients_dict, passed_QC, out_path):
    for patient in tqdm(patients_dict):
        df_list = create_df_list(patients_dict[patient], passed_QC[patient])
        all_patient_chr = combine_df(df_list)
        save_to_file(patient, all_patient_chr, out_path)


def create_df_list(patient_dict, passed_QC):
    patient_list = []
    for chr in patient_dict:
        chr_df = patient_dict[chr]
        # filter
        filtered_df = chr_df.loc[chr_df.index & passed_QC, :]
        # transpose
        transposed_df = filtered_df.T
        # add bed columns
        transposed_df.insert(0, 'chr', chr)
        transposed_df.insert(1, 'start', transposed_df.index.astype(int))
        transposed_df.insert(2, 'end', transposed_df.start + 1)
        patient_list.append(transposed_df)
    return patient_list


def main():
    args = parse_input()

    passed_QC = get_QC(args.QC_path)

    all_files = files_tools.get_files_to_work(args.methylation_folder, pattern=os.path.join("*", "*.pkl.zip"))
    patients_dict = get_patient_df_dict(all_files)
    create_bed(patients_dict, passed_QC, args.output_folder)


if __name__ == '__main__':
    t0 = time.time()
    main()
    t1 = time.time()
    print("total time:", t1 - t0)
