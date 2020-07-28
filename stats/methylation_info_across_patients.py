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
from commons import files_tools
from classifier import create_data
from covariance import covariance_to_bedgraph

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
SEQ_SIZE = 150


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
                                                 sublineage_name=covariance_to_bedgraph.ONLY_PT)
    mdf = df.mean()
    return mdf.loc[indexes].values


def main():
    args = parse_input()

    methylation_folder = args.methylation_folder
    all_files_dict = get_bedgraph_in_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    cov_dict = {}

    for chromosome in all_files_dict:
        cov_dict[chromosome] = []
        try:
            chromosome = str(chromosome)
            windows_data = global_windows_data[chromosome]
        except Exception:
            chromosome = int(chromosome)
            windows_data = global_windows_data[chromosome]

        for file_path in all_files_dict[chromosome]:
            patient, _ = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
            dff = pd.read_pickle(file_path)
            try:
                chromosome = str(chromosome)
                covariance_pmd_df = handle_pmds.filtered_out_non_pmd(dff, chromosome)
            except:
                chromosome = int(chromosome)
                covariance_pmd_df = handle_pmds.filtered_out_non_pmd(dff, chromosome)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.columns >= start) & (covariance_pmd_df.columns <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            _, df = covariance_to_bedgraph.get_region_df(covariance_pmd_df.loc[:, prev_mask],
                                                         sublineage_cells=[],
                                                         sublineage_name=covariance_to_bedgraph.ONLY_PT)
            mdf = df.mean()
            cov_dict[chromosome].append((patient, mdf))

    sum_list = []
    for chromosome in cov_dict:
        sum_df = pd.DataFrame(columns=["chromosome", "location"])
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
            sum_df.loc[cov_df.index, "meth%s" % patient[-2:]] = \
                get_cancer_methylation_of_patient(patient, chromosome, cov_df.index, args.methylation_folder)

        sum_list.append(sum_df)

    sum_df = pd.concat(sum_list)
    try:
        sum_df.to_pickle(os.path.join(args.output_folder, "all_patients_pmd_meth.pkl"))
        sum_df.to_csv(os.path.join(args.output_folder, "all_patients_pmd_meth.pkl"))
    except:
        sum_df.to_pickle("all_patients_pmd_meth.pkl")


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
