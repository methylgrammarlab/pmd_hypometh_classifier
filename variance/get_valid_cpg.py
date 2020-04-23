import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import pyfaidx

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import files_tools, consts, slurm_tools
from covariance import covariance_to_bedgraph
from format_files import handle_pmds

MEAN = 1
VAR = 2
METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

SEQ_SIZE = 150

genome = pyfaidx.Fasta(consts.GENOME_FILE_LOCAL_DROR, sequence_always_upper=True, as_raw=True)


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--nc_files', help='Path to nc files', required=True)
    parser.add_argument('--cells_to_use', help='Path to cells to use file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


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


def get_seq_for_cpg(chr_num, i, seq_size):
    chr_info = genome[chr_num]
    # assert chr_info[i] == "G"
    seq_size = int((seq_size - 2) / 2)
    return chr_info[i - seq_size - 1:i + seq_size + 1]


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


def get_seq_info(ind, chromosome):
    seq = []
    for i in ind:
        seq.append(get_seq_for_cpg(chromosome, i, SEQ_SIZE))

    return seq


def main():
    args = parse_input()
    cov_dict = {}

    methylation_folder = args.methylation_folder
    cells_to_use_path = args.cells_to_use
    if not cells_to_use_path:
        sys.exit("missing the following arg: --cells_to_use")

    cells_to_use = files_tools.load_compressed_pickle(cells_to_use_path)

    all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    all_files_dict_2 = {}

    for k in all_files_dict:
        if k in cells_to_use:
            all_files_dict_2[k] = [all_files_dict[k][0]]

    all_files_dict = all_files_dict_2
    patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                global_windows_data,
                                                                                )

    data_file = open("summary_info.csv", "w")
    data_file.write("patient, chromosome, total_cpg, low_filter, high_filter, both_filter\n")

    for patient in patients_dict:
        patient_cells = cells_to_use[patient]
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            total_cells = []
            low_cells = patient_cells["low"]
            high_cells = patient_cells["high"]
            total_cells.extend(low_cells)
            total_cells.extend(high_cells)

            filtered_df = df.loc[total_cells]
            cpg_coverage_low = np.sum(~pd.isnull(df.loc[low_cells]), axis=0)
            cpg_coverage_high = np.sum(~pd.isnull(df.loc[high_cells]), axis=0)
            low_filter = cpg_coverage_low >= 5
            high_filter = cpg_coverage_high >= 5
            both_filter = np.logical_and(cpg_coverage_high >= 5, cpg_coverage_low >= 5)

            data_file.write("%s,%s,%s,%s,%s,%s\n" %
                            (
                            patient, chromosome, low_filter.shape[0], np.sum(low_filter), np.sum(high_filter),
                            np.sum(both_filter)))

            filtered_df = filtered_df.loc[:, both_filter]  # we only want cpg with at least 5 points
            if chromosome not in cov_dict:
                cov_dict[chromosome] = []

            cov_dict[chromosome].append((patient, filtered_df))

    data_file.close()
    sum_list = []
    for chromosome in cov_dict:
        ind = set([])
        for patient_info in cov_dict[chromosome]:
            patient, meth_df = patient_info
            meth_df = meth_df.T
            ind = set(meth_df.index.values) | ind

        indexes = list(ind)
        indexes.sort()
        sum_df = pd.DataFrame(columns=["chromosome", "location"])
        sum_df["location"] = indexes
        sum_df["chromosome"] = chromosome
        sum_df["sequence"] = get_seq_info(indexes, chromosome)

        sum_df = sum_df.set_index("location")
        for patient_info in cov_dict[chromosome]:
            patient, meth_df = patient_info

            nc_meth_avg = get_nc_avg(chromosome, meth_df.columns, args.nc_files)
            sum_df.loc[nc_meth_avg.index, "nc_avg"] = nc_meth_avg

            pmd_index = handle_pmds.get_pmd_index(sum_df, chromosome, global_windows_data)
            sum_df.loc[pmd_index.index, "pmd_index"] = pmd_index

            sum_df.loc[meth_df.columns, "meth%s" % patient[-2:]] = meth_df.mean()
            sum_df.loc[meth_df.columns, "var%s" % patient[-2:]] = meth_df.var()
            nc_methylation = get_methylation_of_patient(patient, chromosome, meth_df.columns,
                                                        args.methylation_folder,
                                                        sublineage_name=covariance_to_bedgraph.ONLY_NC)
            sum_df.loc[nc_methylation.index, "nc_meth%s" % patient[-2:]] = nc_methylation

        sum_list.append(sum_df.reset_index())

    final_df = pd.concat(sum_list)
    try:
        final_df.to_pickle(os.path.join(args.output_folder, "valid_cpg.pkl"))
        # final_df.to_csv(os.path.join(args.output_folder, "valid_cpg.csv" ))
    except:
        final_df.to_pickle("valid_cpg.pkl")


if __name__ == '__main__':
    slurm_tools.init_slurm(main)
    # main()
