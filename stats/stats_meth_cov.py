import argparse
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

plt.style.use('seaborn')

sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
BEDGRAPH_FILE_NAME_RE = re.compile(".*_chr_(\d+)_.*")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--covariance_folder', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def format_args():
    """
    Format the args for this script
    :return: The path of the files and the output directory
    """
    args = parse_input()
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])
    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    covariance_dict = get_covariance_dict(args.covariance_folder)

    return all_cpg_format_file_paths, covariance_dict, output


def get_covariance_dict(covariance_path):
    covariance_path = os.path.join(covariance_path, "*.bedgraph")
    d = {}
    bedgraph_files = glob.glob(covariance_path)
    for f in bedgraph_files:
        chromosome = BEDGRAPH_FILE_NAME_RE.findall(f)[0]
        d["chr%s" % chromosome] = f

    return d


def plot_methylation_vs_covariance(input_files, covariance_dict, high_threshold=10.0, low_threshold=-10.0):
    mean_list = []
    covariance_list = []

    for file_path in tqdm(input_files):
        _, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]
        covariance_pmd_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(covariance_dict[chromosome],
                                                                               chromosome)

        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage < high_threshold]
        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage > low_threshold]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.filtered_out_non_pmd(df, chromosome)
        pmd_mean = pmd_df.mean(axis=0, skipna=True)
        mean_values = pmd_mean[covariance_pmd_df.index]._values
        covariance_values = covariance_pmd_df._values
        mean_list.append(mean_values)
        covariance_list.append(covariance_values)

    mean_c = np.concatenate(mean_list)
    covariance_c = np.concatenate(covariance_list)
    plt.plot(mean_c, covariance_c, linestyle='', marker='o', markersize=0.5)
    plt.title("Methylation level vs Covariance in PMD")
    plt.xlabel("Avg methylation level")
    plt.ylabel("Covariance in window")
    plt.savefig("methylation_vs_covariance")
    plt.close()


def plot_methylation_diff_vs_covariance(input_files, covariance_dict, high_threshold=10.0,
                                        low_threshold=-10.0):
    mean_list = []
    covariance_list = []

    for file_path in tqdm(input_files):
        _, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]
        covariance_pmd_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(covariance_dict[chromosome],
                                                                               chromosome)

        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage < high_threshold]
        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage > low_threshold]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.filtered_out_non_pmd(df, chromosome)

        nc_index = [cell_id for cell_id in pmd_df.index if cell_id.startswith('NC')]
        pt_index = [cell_id for cell_id in pmd_df.index if cell_id.startswith('PT')]

        nc_values = pmd_df.loc[nc_index, :].mean(axis=0, skipna=True)
        pt_values = pmd_df.loc[pt_index, :].mean(axis=0, skipna=True)
        diff = nc_values - pt_values

        mean_values = diff[covariance_pmd_df.index]._values
        covariance_values = covariance_pmd_df._values
        mean_list.append(mean_values)
        covariance_list.append(covariance_values)

    mean_c = np.concatenate(mean_list)
    covariance_c = np.concatenate(covariance_list)

    plt.plot(mean_c, covariance_c, linestyle='', marker='o', markersize=0.5)
    plt.title("Diff in Methylation level vs Covariance in PMD")
    plt.xlabel("Avg methylation diff level")
    plt.ylabel("Covariance in window")
    plt.savefig("diff_methylation_vs_covariance")
    plt.close()


def plot_covariance_hist(covariance_dict):
    covariance_list = []
    for chromosome in covariance_dict:
        covariance_pmd_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(covariance_dict[chromosome],
                                                                               chromosome)
        covariance_list.append(covariance_pmd_df._values)

    data = np.concatenate(covariance_list)
    plt.hist(data, bins=20)
    plt.savefig("covariance_hist")


def plot_cancer_covariance_vs_meth(input_files, covariance_dict, only_pmd=False):
    mean_list = []
    covariance_list = []
    title = "across all genome" if not only_pmd else "across pmd"

    for file_path in tqdm(input_files):
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]
        patient_df = pd.read_pickle(file_path)

        if only_pmd:
            covariance_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(covariance_dict[chromosome],
                                                                               chromosome)
            patient_df = handle_pmds.filtered_out_non_pmd(patient_df, chromosome)
        else:
            covariance_df = files_tools.load_bedgraph(covariance_dict[chromosome])

        cancer_samples_id = [cell_id for cell_id in patient_df.index if not cell_id.startswith('NC')]
        cancer_df = patient_df.loc[cancer_samples_id, :]

        z_cov = np.abs(stats.zscore(covariance_df._values))
        covariance_df_o = covariance_df[(z_cov < 3)]  # Remove outliers
        average = np.mean(cancer_df[covariance_df_o.index], axis=0)

        average = np.round(average * 500).astype(np.int) / 500
        # covariance_df_o = np.round(covariance_df_o * 100).astype(np.int) / 100
        mean_list.append(average)
        covariance_list.append(covariance_df_o)

    mean_c = np.concatenate(mean_list)
    covariance_c = np.concatenate(covariance_list)

    new_df = pd.DataFrame(
        data=np.hstack((mean_c.reshape(covariance_c.size, 1), covariance_c.reshape(covariance_c.size, 1))),
        columns=["methylation", "covariance"])

    plt.plot(new_df.groupby("methylation").mean(), linestyle='', marker='o', markersize=2)
    plt.title("Avg Methylation level vs Avg Covariance, %s, %s" % (patient, title))
    plt.xlabel("Avg methylation")
    plt.ylabel("Avg Covariance value")
    plt.savefig("methylation_vs_covariance_%s_%s" % (patient, title.replace(" ", "_")))
    plt.close()


def main():
    input_files, covariance_dict, output_dir = format_args()
    # plot_methylation_vs_covariance(input_files, covariance_dict, high_threshold=0.08, low_threshold=-0.03)
    # plot_methylation_diff_vs_covariance(input_files, covariance_dict, high_threshold=0.08,
    #                                     low_threshold=-0.03)
    # plot_covariance_hist(covariance_dict)
    plot_cancer_covariance_vs_meth(input_files, covariance_dict, only_pmd=False)
    plot_cancer_covariance_vs_meth(input_files, covariance_dict, only_pmd=True)


if __name__ == '__main__':
    main()
