import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools, slurm_tools
from classifier import create_data
from covariance import covariance_to_bedgraph

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr_(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
SEQ_SIZE = 150


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedgraph_files', help='Path to bedgraph files', required=True)
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


def get_bedgraph_in_dict(all_file_paths):
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


def get_seq_info(ind, chromosome):
    seq = []
    for i in ind:
        seq.append(create_data.get_seq_for_cpg(chromosome, i, SEQ_SIZE))

    return seq


def get_cancer_methylation_of_patient(patient, chromosome, indexes, methylation_folder):
    methylation_file_path = os.path.join(methylation_folder, patient,
                                         METHYLATION_FILE_FORMAT % (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    _, df = covariance_to_bedgraph.get_region_df(df, sublineage_cells=[],
                                                 sublineage_name=covariance_to_bedgraph.ALL_CANCER)
    mdf = df.mean()
    return mdf.loc[indexes].values


def main():
    args = parse_input()

    all_file_paths = get_bedgraph_files(args.bedgraph_files)
    all_files_dict = get_bedgraph_in_dict(all_file_paths)

    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    cov_dict = {}

    for chromosome in all_files_dict:
        cov_dict[chromosome] = []
        windows_data = global_windows_data[chromosome]
        for file_path in all_files_dict[chromosome]:
            patient, _ = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
            covariance_pmd_df = handle_pmds.get_covariance_pmd_df(file_path, chromosome, True)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.index >= start) & (covariance_pmd_df.index <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            cov_dict[chromosome].append((patient, covariance_pmd_df.loc[prev_mask, :]))

    sum_list = []
    for chromosome in cov_dict:
        sum_df = pd.DataFrame(columns=["chromosome", "location", "pmd_index"])
        ind = None
        for patient_info in cov_dict[chromosome]:
            patient, cov_df = patient_info
            if ind is not None:
                ind = set(cov_df.index.values) | ind
            else:
                ind = set(cov_df.index.values)

        sum_df["location"] = list(ind)
        sum_df["chromosome"] = chromosome
        sum_df["sequence"] = get_seq_info(ind, chromosome)

        sum_df = sum_df.set_index("location")
        for patient_info in cov_dict[chromosome]:
            patient, cov_df = patient_info
            sum_df.loc[cov_df.index, "cov%s" % patient[-2:]] = cov_df["coverage"]
            sum_df.loc[cov_df.index, "pmd_index"] = cov_df["pmd_index"]
            sum_df.loc[cov_df.index, "meth%s" % patient[-2:]] = \
                get_cancer_methylation_of_patient(patient, chromosome, cov_df.index, args.methylation_folder)

        sum_list.append(sum_df)

    sum_df = pd.concat(sum_list)
    try:
        sum_df.to_pickle("/vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror"
                         "/covariance/cancer/5000_windows_2500_overlap_750_min_pairs/valid_cpg.pkl")
        sum_df.to_csv("/vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror"
                      "/covariance/cancer/5000_windows_2500_overlap_750_min_pairs/valid_cpg.csv")
    except:
        sum_df.to_pickle("valid_cpg.pkl")

if __name__ == '__main__':
    slurm_tools.init_slurm(main)
