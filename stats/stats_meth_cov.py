import argparse
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.getcwd())
from format_files import handle_pmds

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
        covariance_pmd_df = handle_pmds.get_covariance_pmd_df(covariance_dict[chromosome], chromosome)

        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage < high_threshold]
        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage > low_threshold]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.get_pmd_df(df, chromosome)
        pmd_mean = pmd_df.mean(axis=0, skipna=True)
        mean_values = pmd_mean[covariance_pmd_df.index]._values
        covariance_values = covariance_pmd_df._values
        mean_list.append(mean_values)
        covariance_list.append(covariance_values)

    plt.scatter(np.concatenate(mean_list), np.concatenate(covariance_list), alpha=0.2, s=0.8)
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
        covariance_pmd_df = handle_pmds.get_covariance_pmd_df(covariance_dict[chromosome], chromosome)

        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage < high_threshold]
        covariance_pmd_df = covariance_pmd_df[covariance_pmd_df.coverage > low_threshold]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.get_pmd_df(df, chromosome)

        nc_index = [cell_id for cell_id in pmd_df.index if cell_id.startswith('NC')]
        pt_index = [cell_id for cell_id in pmd_df.index if cell_id.startswith('PT')]

        nc_values = pmd_df.loc[nc_index, :].mean(axis=0, skipna=True)
        pt_values = pmd_df.loc[pt_index, :].mean(axis=0, skipna=True)
        diff = nc_values - pt_values

        mean_values = diff[covariance_pmd_df.index]._values
        covariance_values = covariance_pmd_df._values
        mean_list.append(mean_values)
        covariance_list.append(covariance_values)

    plt.scatter(np.concatenate(mean_list), np.concatenate(covariance_list), alpha=0.2, s=0.8)
    plt.title("Diff in Methylation level vs Covariance in PMD")
    plt.xlabel("Avg methylation diff level")
    plt.ylabel("Covariance in window")
    plt.savefig("diff_methylation_vs_covariance")
    plt.close()


def plot_covariance_hist(covariance_dict):
    covariance_list = []
    for chromosome in covariance_dict:
        covariance_pmd_df = handle_pmds.get_covariance_pmd_df(covariance_dict[chromosome], chromosome)
        covariance_list.append(covariance_pmd_df._values)

    data = np.concatenate(covariance_list)
    plt.hist(data, bins=20)
    plt.savefig("covariance_hist")


def main():
    input_files, covariance_dict, output_dir = format_args()
    plot_methylation_vs_covariance(input_files, covariance_dict, high_threshold=0.1, low_threshold=-0.05)
    plot_methylation_diff_vs_covariance(input_files, covariance_dict, high_threshold=0.1,
                                        low_threshold=-0.05)
    # plot_covariance_hist(covariance_dict)


if __name__ == '__main__':
    main()
