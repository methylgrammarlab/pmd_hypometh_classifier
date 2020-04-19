import argparse
import glob
import os
import re
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from covariance import covariance_to_bedgraph

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
SEQ_SIZE = 150
TOP_LOW_PERCENTAGE_TO_REMOVE = 5


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_bedgraph_files(files):
    if os.path.isdir(files):
        file_path = os.path.join(files, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [files]

    return all_file_paths


def get_patient_dict(all_file_paths):
    d = {}
    for file_path in all_file_paths:
        try:
            patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
        except:
            continue

        if patient not in d:
            d[patient] = []

        d[patient].append((chromosome, file_path))

    return d


def get_cancer_methylation_of_patient(patient, chromosome, indexes, methylation_folder):
    methylation_file_path = os.path.join(methylation_folder, patient,
                                         METHYLATION_FILE_FORMAT % (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    _, df = covariance_to_bedgraph.get_region_df(df, sublineage_cells=[],
                                                 sublineage_name=covariance_to_bedgraph.ONLY_PT)
    mdf = df.mean()
    return mdf.loc[indexes].values


def main():
    args = parse_input()
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    methylation_folder = args.methylation_folder
    # all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "all*.pkl.zip")))
    all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)

    # get df per patient only for valid PMD
    for patient in all_files_dict:
        if patient not in ["CRC01", "CRC04", "CRC10", "CRC13", "CRC14"]:
            continue
        patients_dict = {}
        for t in all_files_dict[patient]:
            chromosome, file_path = t
            try:
                chromosome = str(chromosome)
                windows_data = global_windows_data[chromosome]
            except Exception:
                chromosome = int(chromosome)
                windows_data = global_windows_data[chromosome]

            patient_chr_df = pd.read_pickle(file_path)
            try:
                chromosome = str(chromosome)
                covariance_pmd_df = handle_pmds.get_pmd_df(patient_chr_df, chromosome)
            except:
                chromosome = int(chromosome)
                covariance_pmd_df = handle_pmds.get_pmd_df(patient_chr_df, chromosome)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.columns >= start) & (covariance_pmd_df.columns <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            _, df = covariance_to_bedgraph.get_region_df(covariance_pmd_df.loc[:, prev_mask],
                                                         sublineage_cells=[],
                                                         sublineage_name=covariance_to_bedgraph.ALL_CANCER)

            cell_mean = df.mean(axis=1).sort_values()
            bottom = min(len(cell_mean) - 1, int((len(cell_mean) / 100) * 15))
            top = int((len(cell_mean) / 100) * 15)
            low_15 = list(cell_mean.iloc[:bottom].index)
            high_15 = list(cell_mean.iloc[-top - 1:-1].index)
            patients_dict[chromosome] = {"low": low_15, "high": high_15}
        output_path = os.path.join(output, "%s_15perc_low_high_avg_methylation_cells.pickle.zlib" % patient)
        files_tools.save_as_compressed_pickle(output_path, patients_dict)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
