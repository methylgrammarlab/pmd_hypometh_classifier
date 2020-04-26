import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')

PATIENTS = ["01", "10", "11", "13"]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_file', help='Path to the cpg file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def global_cpg_info(df):
    # Give some statistics on the data

    num_of_pmd = df.groupby(["chromosome", "pmd_index"]).ngroups
    num_of_unique_seq = len(df["sequence"].unique())

    print("##########\nprinting information\n#######")
    print("%s CpG passed PMD filter (%s PMDs). %s are unique seq" %
          (df.shape[0], num_of_pmd, num_of_unique_seq))

    solo = [seq for seq in df["sequence"] if seq.count("CG") == 1]
    weak_solo = [seq for seq in solo if seq[73] in ["A", "T"] and seq[76] in ["A", "T"]]
    strong_solo = [seq for seq in solo if seq[73] in ["C", "G"] and seq[76] in ["C", "G"]]

    print("We have %s(%.2f%%) solo CpG" % (len(solo), len(solo) / df.shape[0] * 100))
    print("solo WCGW: %s(%s) and solo SCGS: %s(%s)" %
          (len(weak_solo), len(weak_solo) / len(solo) * 100,
           len(strong_solo), len(strong_solo) / len(solo) * 100))

    for patient in PATIENTS:
        label = "var%s" % patient

        num_of_empty = np.sum(pd.isnull(df[label]))
        print("Crc%s has %s(%s) cpg without data" % (patient, num_of_empty, num_of_empty / df.shape[0] * 100))


def plot_variance_histogram(df, output):
    for patient in PATIENTS:
        label = "var%s" % patient
        meth_label = "meth%s" % patient
        nc_label = "nc_meth%s" % patient
        df = df[df[meth_label] < 0.8]
        df = df[df[nc_label] > 0.5]

        _ = plt.hist(df[label], color="b")
        plt.style.use('ggplot')
        plt.title("Hist of variance in PMD for crc%s. Median(%s)" % (
            patient, np.median(df[label][~pd.isnull(df[label])])))
        plt.xlabel("Variance value")
        plt.ylabel("Amount")
        plt.savefig(os.path.join(output, "hist_variance_in_pmd_crc%s.png" % patient))
        plt.close()


def plot_variance_histogram_vs_type(df, output):
    for patient in PATIENTS:
        label = "var%s" % patient
        meth_label = "meth%s" % patient
        nc_label = "nc_meth%s" % patient
        df = df[df[meth_label] < 0.8]
        df = df[df[nc_label] > 0.5]

        strong_rows = df[df["sequence"].str.contains("[CG]CG[CG]", regex=True)]
        weak_rows = df[df["sequence"].str.contains("[AT]CG[AT]", regex=True)]
        _ = plt.hist(weak_rows[label], label="WCGW", normed=True, alpha=0.5, color="b")
        _ = plt.hist(strong_rows[label], label="SCGS", normed=True, alpha=0.5, color="g")
        plt.legend()
        plt.style.use('ggplot')
        plt.title("Histogram of variance in strong vs weak for crc%s" % (patient))
        plt.xlabel("Variance")
        plt.ylabel("Amount")
        plt.savefig(os.path.join(output, "var_strong_vs_weak_crc%s.png" % patient))
        plt.close()


def plot_var_density(df, output):
    for patient in PATIENTS:
        label = "var%s" % patient
        meth_label = "meth%s" % patient
        nc_label = "nc_meth%s" % patient
        df = df[df[meth_label] < 0.8]
        df = df[df[nc_label] > 0.5]

        strong_rows = df[df["sequence"].str.contains("[CG]CG[CG]", regex=True)]
        weak_rows = df[df["sequence"].str.contains("[AT]CG[AT]", regex=True)]

        sns.distplot(weak_rows[label][~pd.isnull(weak_rows[label])],
                     hist=False, kde=True, kde_kws={'linewidth': 3}, label="WCGW")
        sns.distplot(strong_rows[label][~pd.isnull(strong_rows[label])],
                     hist=False, kde=True, kde_kws={'linewidth': 3}, label="SSCS")
        plt.title("Variance density for CRC%s" % patient)
        plt.savefig(os.path.join(output, "density_var_strong_vs_weak_crc%s.png" % patient))
        plt.close()


def plot_meth_vs_var_scatter(df, output):
    for patient in PATIENTS:
        var_label = "var%s" % patient
        meth_label = "meth%s" % patient
        # df = df[df[meth_label] < 0.8]
        # plt.plot(df[meth_label], df[var_label], linestyle='', marker='o', markersize=0.5)
        # plt.title("Methylation level vs Covariance in solo CpG for patient %s" % patient)
        # plt.xlabel("Methylation")
        # plt.ylabel("Variance")
        # plt.savefig(os.path.join(output, "solo_meth_vs_var_scatter_patient%s.png" % patient))
        # plt.close()

        sns_plot = sns.jointplot(x=meth_label, y=var_label, data=df, kind="kde")
        plt.title("Methylation level vs Covariance in solo CpG for patient %s" % patient)
        plt.xlabel("Methylation")
        plt.ylabel("Variance")
        sns_plot.savefig(os.path.join(output, "dist_solo_meth_vs_var_scatter_patient%s.png" % patient))
        plt.close()


def get_basic_info(df, name, output_folder):
    output_folder = os.path.join(output_folder, name)
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    #
    plot_variance_histogram(df, output_folder)
    plot_variance_histogram_vs_type(df, output_folder)
    plot_var_density(df, output_folder)


def remove_by_nc_methylation_info(df):
    for patient in PATIENTS:
        nc_label = "nc_avg"
        var_label = "var%s" % patient
        dfc = df.copy()
        dfc = dfc[~pd.isnull(dfc[var_label])]
        amount_of_cpg = dfc.shape[0]
        amount_of_cpg_with_nc = np.sum(dfc[nc_label] > 0.5)

        print("CRC%s using methylation in nc >0.5 will give %s%% cpg which is %s" %
              (patient, amount_of_cpg_with_nc / amount_of_cpg * 100, amount_of_cpg))


def main():
    args = parse_input()
    df = pd.read_pickle(args.cpg_file)
    df["small_seq"] = df["sequence"].str[73:77]

    # Remove empty cpg
    methylation_columns = df[["meth%s" % i for i in PATIENTS]]
    df = df[~pd.isnull(methylation_columns.min(axis=1))]

    # keep only solo
    solo_rows = df[df["sequence"].str.count("CG") == 1]

    remove_by_nc_methylation_info(solo_rows)

    solo_after_nc_rows = solo_rows[solo_rows["nc_avg"] > 0.5]

    # just solo
    # global_cpg_info(df)
    global_cpg_info(solo_after_nc_rows)
    plot_meth_vs_var_scatter(solo_after_nc_rows, args.output_folder)

    # get_basic_info(solo_rows, "solo", args.output_folder)


if __name__ == '__main__':
    main()
