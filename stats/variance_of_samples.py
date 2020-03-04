import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
# from pandas_profiling import ProfileReport
from tqdm import tqdm

sys.path.append(os.getcwd())
from format_files import format_sublineage_info
from format_files.format_cpg_context_map import NUMBER_OF_ORPH_PER_INDEX
from format_files import format_cpg_context_map

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
SUB_LINEAGE_FILE = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\resource\sublineage\patient_sublineage_dict.pickle"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
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

    return all_cpg_format_file_paths, output


def plot_variance_histogram(df, patient, chromosome, group_name, group_samples, include_nc = False):
    df.reset_index(drop=True)
    samples = df.axes[0]
    pt_index = []
    nc_avg = None

    for sample in group_samples:
        pt_index.extend([i for i in range(len(samples)) if samples[i].startswith(sample)])

    if include_nc:
        nc_index = [i for i in range(len(samples)) if samples[i].startswith("NC")]
        nc_variance = df.iloc[nc_index, :].var(axis=0, skipna=True)
        nc_avg = nc_variance[~nc_variance.isnull()]._values.mean()
        plt.hist(nc_variance, 30, facecolor='b', label="nc", log=True)

    pt_variance = df.iloc[pt_index, :].var(axis=0, skipna=True)

    pt_avg = pt_variance[~pt_variance.isnull()]._values.mean()

    plt.hist(pt_variance, 30, facecolor='r', alpha=0.5, label="pt", log=True)

    plt.title("Histogram of variance nc(%s) vs %s(%s)" % (nc_avg, group_name, pt_avg))
    plt.xlabel("Variance value")
    plt.ylabel("Freq logged")
    plt.legend()
    plt.savefig("variance_of_%s_%s_nc_vs_%s.png" % (patient, chromosome, group_name))
    plt.close()

    return nc_avg, pt_avg


def plot_cancer_hist(df, patient, chromosome):
    df.reset_index(drop=True)
    samples = df.axes[0]
    pt_index = []

    pt_index.extend([i for i in range(len(samples)) if samples[i].startswith("LN") or samples[i].startswith("PT")])
    pt_variance = df.iloc[pt_index, :].var(axis=0, skipna=True)
    pt_avg = pt_variance[~pt_variance.isnull()]._values.mean()
    plt.hist(pt_variance, bins=30, facecolor='r', alpha=0.5, log=True)

    plt.title("Histogram of variance of cancer %s" %pt_avg)
    plt.xlabel("Variance value")
    plt.ylabel("Freq logged")
    plt.savefig("variance_of_%s_%s_cancer.png" % (patient, chromosome))
    plt.close()


def profile_large_variance(chromosome, cpg_dict, nc_variance, patient):
    strange_loc = cpg_dict[nc_variance == 0.5]
    orph_coloms = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]
    context_3 = format_cpg_context_map.get_context_as_str_for_chr(strange_loc)
    context_2 = format_cpg_context_map.get_context_as_str_for_chr_2(strange_loc)
    context_1 = format_cpg_context_map.get_context_as_str_for_chr_1(strange_loc)
    orph_info = strange_loc[:, 1:14]
    weak = format_cpg_context_map.get_weak_column(strange_loc)
    strong = format_cpg_context_map.get_strong_column(strange_loc)
    weak_or_strong = np.logical_xor(weak, strong)
    # Combine the different tables to one table and than convert to df
    final_table = np.hstack((orph_info, weak[:, None], strong[:, None], weak_or_strong[:, None],
                             np.array(context_1)[:, None], np.array(context_2)[:, None],
                             np.array(context_3)[:, None]))
    end_df = pd.DataFrame(final_table,
                          columns=orph_coloms + ["weak", "strong", "weak_or_strong", "context_1",
                                                 "context_2", "context_3"])
    profile = ProfileReport(end_df, title='Pandas Profiling Report', html={'style': {'full_width': True}})
    profile.to_file(output_file=f"variance_of_{patient}_{chromosome}.html")


def plot_variance_histogram_groups(df, patient, chromosome):
    sublineage_info = format_sublineage_info.get_sublineage_info(SUB_LINEAGE_FILE)
    patient_info = sublineage_info[patient]
    for sublineage in patient_info:
        nc_value, variance_value = plot_variance_histogram(df, patient, chromosome, sublineage,
                                                           patient_info[sublineage])
        print("%s:%s" % (sublineage, variance_value))


def main():
    input_files, output_dir = format_args()

    for file_path in tqdm(input_files):
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]

        df = pd.read_pickle(file_path)
        plot_cancer_hist(df, patient, chromosome)

if __name__ == '__main__':
    main()
