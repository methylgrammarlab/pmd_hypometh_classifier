import argparse
import glob
import os
import pickle
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

    methylation_folder = args.methylation_folder
    all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    patients_dict = {}
    patients_dict_f = {}

    # get df per patient only for valid PMD
    for patient in all_files_dict:
        patients_dict[patient] = []
        patients_dict_f[patient] = []
        for t in all_files_dict[patient]:
            chromosome, file_path = t
            try:
                chromosome = str(chromosome)
                windows_data = global_windows_data[chromosome]
            except Exception:
                chromosome = int(chromosome)
                windows_data = global_windows_data[chromosome]

            dff = pd.read_pickle(file_path)
            try:
                chromosome = str(chromosome)
                covariance_pmd_df = handle_pmds.get_pmd_df(dff, chromosome)
            except:
                chromosome = int(chromosome)
                covariance_pmd_df = handle_pmds.get_pmd_df(dff, chromosome)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.columns >= start) & (covariance_pmd_df.columns <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            _, df = covariance_to_bedgraph.get_region_df(covariance_pmd_df.loc[:, prev_mask],
                                                         sublineage_cells=[],
                                                         sublineage_name=covariance_to_bedgraph.ALL_CANCER)

            coverage = np.sum(~pd.isnull(df), axis=1) / df.shape[1]
            meth = np.mean(df, axis=1)
            patients_dict[patient].append((chromosome, meth, coverage, set(coverage.index)))

            # n_cells = coverage.shape[0]
            # cells_to_ignore = int(n_cells * 2.5 / 100)

    pa_cells = set([])
    for pa in patients_dict:
        for t in patients_dict[pa]:
            cells = t[3]
            pa_cells |= cells

        pa_df_meth = pd.DataFrame(index=list(pa_cells))
        pa_df_cov = pd.DataFrame(index=list(pa_cells))
        for t in patients_dict[pa]:
            chromosome, df, coverage, cells = t
            pa_df_cov.loc[coverage.index, chromosome] = coverage
            pa_df_meth.loc[df.index, chromosome] = df

        cov_sorted = np.mean(pa_df_cov, axis=1).sort_values()
        n_cells = cov_sorted.shape[0]
        cells_to_ignore = int(n_cells * TOP_LOW_PERCENTAGE_TO_REMOVE / 100)
        cells_to_use = cov_sorted.index[cells_to_ignore:-cells_to_ignore]

        patients_dict_f[pa] = (pa_df_cov, pa_df_meth, cells_to_use)

    for pa in patients_dict_f:
        pa_df_meth = patients_dict_f[pa][1]
        cells_to_use = patients_dict_f[pa][2]
        values = pa_df_meth.loc[cells_to_use].median(axis=1).sort_values()
        x = [i for i in range(values.shape[0])]
        plt.bar(x, values)
        plt.xlabel("tumor cells (self sorted)")
        plt.ylabel("avg meth")
        plt.title("cell methylation means for %s" % pa)
        plt.savefig("cell_methylation_means_for_%s.png" % pa)
        plt.close()

    # This save dict with patient name as key and value of:
    # (coverage per chromosome, methylation per chromosome, cells to use after removing 5% of cov)
    with open("methylation_avg_per_cells.pickle", "wb") as output:
        pickle.dump(patients_dict_f, output)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
