"""
Provide some basic statistics on the get_valid_cpgs_dataset file before converting it to a data for the nn
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')
sns.set_style('whitegrid')

# Change this when moving between different datasets
SC_PATIENTS = ["01", "10", "11", "13"]
BULK_PATIENTS = [""]

PATIENTS = SC_PATIENTS


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_file', help='Path to the cpg file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def print_basic_information(df):
    """
    Provide some basic statisic about the CpG in the dataframe
    :param df: The df to work on
    """
    # Give some statistics on the data

    num_of_pmd = df.groupby(["chromosome", "pmd_index"]).ngroups
    num_of_unique_seq = len(df["sequence"].unique())

    print("##########\nprinting information\n#######")
    print("%s CpG passed PMD filter (%s PMDs). %s are unique seq" % (df.shape[0], num_of_pmd,
                                                                     num_of_unique_seq))

    solo = df[df["sequence"].str.count("CG") == 1]
    solo_ss = solo[solo["seq4"].str.contains("[CG]CG[CG]", regex=True)]
    solo_ww = solo[solo["seq4"].str.contains("[AT]CG[AT]", regex=True)]

    print("We have %s(%.2f%%) solo CpG" % (solo.shape[0], solo.shape[0] / df.shape[0] * 100))
    print("solo WCGW: %s(%.2f%%) and solo SCGS: %s(%.2f%%)" %
          (solo_ww.shape[0], solo_ww.shape[0] / solo.shape[0] * 100,
           solo_ss.shape[0], solo_ss.shape[0] / solo.shape[0] * 100))

    for patient in PATIENTS:
        label = "var%s" % patient
        num_of_empty = np.sum(pd.isnull(df[label]))
        print("CRC%s has %s(%s) cpg without data" % (patient, num_of_empty, num_of_empty / df.shape[0] * 100))


def plot_variance_histogram(df, output_folder):
    """
    Plot histogram of variance
    :param df: The df
    :param output_folder: Output folder
    """
    for patient in PATIENTS:
        label = "var%s" % patient
        df[label].hist()

        plt.style.use('ggplot')
        plt.title("Hist of variance in PMD for crc%s" % patient)
        plt.xlabel("Variance value")
        plt.ylabel("Amount")
        plt.savefig(os.path.join(output_folder, "hist_variance_in_pmd_crc%s.png" % patient))
        plt.close()


def plot_variance_histogram_vs_ww_ss(df, output_folder):
    """
    Plot histogram of variance strong solo vs weak solo
    :param df: The df
    :param output_folder: Output folder
    """
    solo = df[df["sequence"].str.count("CG") == 1]
    solo_ss = solo[solo["seq4"].str.contains("[CG]CG[CG]", regex=True)]
    solo_ww = solo[solo["seq4"].str.contains("[AT]CG[AT]", regex=True)]

    for patient in PATIENTS:
        label = "var%s" % patient
        solo_ss[label].hist(label="SCGS")
        solo_ww[label].hist(label="WCGW")

        plt.style.use('ggplot')
        plt.title("Hist of variance in PMD for crc%s SCGS vs WCGW" % patient)
        plt.xlabel("Variance value")
        plt.ylabel("Amount")
        plt.legend()
        plt.savefig(os.path.join(output_folder, "hist_variance_in_pmd_crc%s_ss_ww.png" % patient))
        plt.close()


def plot_meth_density_violin(df, output_folder):
    """
    Plot the methylation value as a violin plot for different density
    :param df: The data frame
    :param output_folder: Output folder
    """
    for patient in PATIENTS:
        meth_label = "meth%s" % patient
        new_df = pd.DataFrame()

        for i in range(1, 6):
            if i != 5:
                df.loc[df["sequence"].str.count("CG") == i, "#CpG"] = "%s" % i
            else:
                df.loc[df["sequence"].str.count("CG") >= i, "#CpG"] = ">%s" % i

        new_df["Methylation level"] = df[meth_label]
        new_df["#CpG"] = df["#CpG"]
        sns.violinplot(y="Methylation level", x="#CpG", data=new_df, palette="muted", order=["1", "2", "3",
                                                                                             "4", ">5"])
        plt.title("Methylation density for for CRC%s" % patient)
        plt.savefig(os.path.join(output_folder, "density_meth_crc%s.png" % patient))
        plt.close()


def plot_meth_density(df, output_folder):
    """
    Plot the methylation density as a facotr of the cpg density
    :param df: The data frame
    :param output_folder: Output folder
    """
    for patient in PATIENTS:
        meth_label = "meth%s" % patient

        for i in range(1, 6):
            if i == 1:
                sns.distplot(df[df["sequence"].str.count("CG") == i][meth_label], hist=False, kde=True,
                             kde_kws={'linewidth': 3}, label="#CpG=%s (solo)" % (i - 1))
            elif i != 5:
                sns.distplot(df[df["sequence"].str.count("CG") == i][meth_label], hist=False, kde=True,
                             kde_kws={'linewidth': 3}, label="#CpG=%s" % (i - 1))

            else:
                sns.distplot(df[df["sequence"].str.count("CG") >= i][meth_label], hist=False, kde=True,
                             kde_kws={'linewidth': 3}, label="#CpG>=%s" % (i - 1))

        plt.title("Methylation density for CRC%s" % patient, fontsize=20)
        plt.xlabel("Methylation Level", fontsize=16)
        plt.grid(False)
        plt.ylabel("Distribution", fontsize=16)
        plt.legend()

        plt.savefig(os.path.join(output_folder, "density_meth_crc%s.png" % patient))
        plt.close()


def plot_meth_vs_var_jointplot(df, output_folder, i=""):
    """
    Plot methylation vs variance of all patients as a jointplot
    :param df: The df
    :param output_folder: path for output folder
    """
    for patient in PATIENTS:
        var_label = "var%s" % patient
        meth_label = "meth%s" % patient

        plt.subplots_adjust(top=0.9)
        sns_plot = sns.jointplot(x=df[meth_label], y=df[var_label], kind="kde")
        sns_plot.fig.suptitle("Methylation vs Variance in solo CpG for CRC%s" % patient, fontsize=20)
        sns_plot.set_axis_labels("Methylation level", "Variance", fontsize=16)

        sns_plot.savefig(os.path.join(output_folder, "dist_meth_vs_var_scatter_patient%s_%s.png" % (
            patient, i)))
        plt.close()


def remove_by_nc_methylation_info(df, nc_filter=0.6):
    """
    Check how many CpG will be removed with nc filter
    :param df: The df
    :param nc_filter: the nc filter value
    :return:
    """
    for patient in PATIENTS:
        nc_label = "nc_avg"
        var_label = "var%s" % patient
        patient_df = df[~pd.isnull(df[var_label])]
        amount_of_cpg = patient_df.shape[0]
        amount_of_cpg_with_nc = np.sum(patient_df[nc_label] >= nc_filter)

        print("CRC%s using methylation normal cell filter>=%s will give %s%% cpg which is %s" %
              (nc_filter, patient, amount_of_cpg_with_nc / amount_of_cpg * 100, amount_of_cpg))


def main():
    args = parse_input()
    df = pd.read_pickle(args.cpg_file)
    df["small_seq"] = df["sequence"].str[73:77]

    # Remove empty cpg
    methylation_columns = df[["meth%s" % i for i in PATIENTS]]
    df = df[~pd.isnull(methylation_columns.min(axis=1))]
    df = df[df["nc_avg"] > 0.6]

    # keep only solo
    solo_rows = df[df["sequence"].str.count("CG") == 1]

    # just solo
    # print_basic_information(df)
    # print_basic_information(solo_after_nc_rows)
    # plot_meth_vs_var_jointplot(solo_after_nc_rows, args.output_folder)
    # plot_meth_density(df, args.output_folder)


if __name__ == '__main__':
    main()
