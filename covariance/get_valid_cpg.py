import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import pyfaidx

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from covariance import covariance_to_bedgraph
from commons import consts

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr_(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
SEQ_SIZE = 150
MEAN = 1
VAR = 2

genome = pyfaidx.Fasta(consts.GENOME_FILE_LOCAL_DROR, sequence_always_upper=True, as_raw=True)


def get_seq_for_cpg(chr_num, i, seq_size):
    chr_info = genome[chr_num]
    # assert chr_info[i] == "G"
    seq_size = int((seq_size - 2) / 2)
    return chr_info[i - seq_size - 1:i + seq_size + 1]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedgraph_files', help='Path to bedgraph files', required=False)
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--nc_files', help='Path to nc files', required=True)
    parser.add_argument('--cells_to_use', help='Path to cells to use file', required=False)
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
            patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]
        except:
            continue

        if patient not in d:
            d[patient] = []

        d[patient].append((chromosome[3:], file_path))

    return d


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
        seq.append(get_seq_for_cpg(chromosome, i, SEQ_SIZE))

    return seq


def get_methylation_of_patient(patient, chromosome, indexes, methylation_folder, sublineage_name, oper=MEAN):
    methylation_file_path = os.path.join(methylation_folder, patient,
                                         METHYLATION_FILE_FORMAT % (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    _, df = covariance_to_bedgraph.get_region_df(df, sublineage_cells=[],
                                                 sublineage_name=sublineage_name)

    if oper == MEAN:
        mdf = df.mean()

    elif oper == VAR:
        mdf = df.var()

    return mdf.loc[indexes]


def get_nc_avg(chromosome, indexes, nc_files):
    f = glob.glob(os.path.join(nc_files, "*chr%s.dummy*" % chromosome))[0]

    df = pd.read_pickle(f)
    mdf = df.mean(axis=1)
    return mdf.loc[indexes]


def main2():
    # trying to create file after file
    args = parse_input()

    all_file_paths = get_bedgraph_files(args.bedgraph_files)
    all_files_dict = get_bedgraph_in_dict(all_file_paths)

    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    cov_dict = {}

    for chromosome in all_files_dict:
        cov_dict[chromosome] = []
        windows_data = global_windows_data[int(chromosome)]
        for file_path in all_files_dict[chromosome]:
            patient, _ = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
            covariance_pmd_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(file_path, chromosome,
                                                                                   True)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.index >= start) & (covariance_pmd_df.index <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            cov_dict[chromosome].append((patient, covariance_pmd_df.loc[prev_mask, :]))

    sum_list = []
    for chromosome in cov_dict:
        ind = None
        for patient_info in cov_dict[chromosome]:
            patient, cov_df = patient_info
            if ind is not None:
                ind = set(cov_df.index.values) | ind
            else:
                ind = set(cov_df.index.values)

        indexes = list(ind)
        indexes.sort()
        sum_df = pd.DataFrame(columns=["chromosome", "location"])
        sum_df["location"] = indexes
        sum_df["chromosome"] = chromosome
        sum_df["sequence"] = get_seq_info(indexes, chromosome)

        sum_df = sum_df.set_index("location")
        sums = []
        for patient_info in cov_dict[chromosome]:
            patient, cov_df = patient_info
            sum_df.loc[cov_df.index, "cov%s" % patient[-2:]] = cov_df["coverage"]
            sum_df.loc[cov_df.index, "pmd_index"] = cov_df["pmd_index"]
            methylation = get_methylation_of_patient(patient, chromosome, cov_df.index,
                                                     args.methylation_folder,
                                                     sublineage_name=covariance_to_bedgraph.ALL_CANCER)
            sum_df.loc[methylation.index, "meth%s" % patient[-2:]] = methylation

            nc_methylation = get_methylation_of_patient(patient, chromosome, cov_df.index,
                                                        args.methylation_folder,
                                                        sublineage_name=covariance_to_bedgraph.ONLY_NC)
            sum_df.loc[nc_methylation.index, "nc_meth%s" % patient[-2:]] = nc_methylation

            cancer_var = get_methylation_of_patient(patient, chromosome, cov_df.index,
                                                    args.methylation_folder,
                                                    sublineage_name=covariance_to_bedgraph.ALL_CANCER,
                                                    oper=VAR)
            sum_df.loc[cancer_var.index, "cancer_var%s" % patient[-2:]] = cancer_var

            nc_meth_avg = get_nc_avg(chromosome, cov_df.index, args.nc_files)
            sum_df.loc[nc_meth_avg.index, "nc_avg"] = nc_meth_avg

        sum_list.append(sum_df.reset_index())

    final_df = pd.concat(sum_list)
    try:
        final_df.to_pickle(os.path.join(args.output_folder, "valid_cpg.pkl"))
        # final_df.to_csv(os.path.join(args.output_folder, "valid_cpg.csv" ))
    except:
        final_df.to_pickle("valid_cpg.pkl")



if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    # slurm_tools.init_slurm(main3)
    main2()
