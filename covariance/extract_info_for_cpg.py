import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('seaborn')


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_file', help='Path to the cpg file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def global_cpg_info(df):
    num_of_pmd = df.groupby(["chromosome", "pmd_index"]).ngroups
    num_of_unique_seq = len(df["sequence"].unique())

    print("We have a total of %s CpG" % df.shape[0])
    print("We have a total of %s PMDs" % num_of_pmd)
    print("We have %s unique sequences" % num_of_unique_seq)

    solo = [seq for seq in df["sequence"] if seq.count("CG") == 1]
    weak_solo = [seq for seq in solo if seq[73] in ["A", "T"] and seq[73] in ["A", "T"]]
    strong_solo = [seq for seq in solo if seq[73] in ["C", "G"] and seq[73] in ["C", "G"]]

    print("We have %s solo CpG which is %s" % (len(solo), len(solo) / df.shape[0] * 100))
    print("Include %s WCGW solo and %s SCGS solo" % (len(weak_solo), len(strong_solo)))


def histogram_of_coverage_across_patients(df):
    cov_columns = df[["cov01", "cov11", "cov13"]]
    values = pd.notnull(cov_columns).sum(axis=1).values
    _ = plt.hist(values, bins='auto')
    plt.style.use('ggplot')
    plt.title("Covariance coverage across samples")
    plt.xlabel("Covered by X patients")
    plt.ylabel("Amount of CpG")
    plt.savefig("hist_of_covariance_cov_across_samples.png")


def histogram_of_num_of_cg(df):
    num_of_cg = [seq.count("CG") for seq in df["sequence"]]
    _ = plt.hist(num_of_cg)
    plt.style.use('ggplot')
    plt.title("Number of CpG in seq include middle. Median (%s)" % np.median(num_of_cg))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_cg_in_seq_include_middle_all.png")

    num_of_cg = [i for i in num_of_cg if i < 5]
    _ = plt.hist(num_of_cg)
    plt.style.use('ggplot')
    plt.title("Number of CpG in seq include middle less than 5. Median (%s)" % np.median(num_of_cg))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_cg_in_seq_include_middle_lower_than_5.png")


def histogram_on_numb_of_cpg_per_pmd(df):
    values = df.groupby(["chromosome", "pmd_index"]).count()["sequence"].values
    values.sort()
    _ = plt.hist(values)
    plt.style.use('ggplot')
    plt.title("Hist of valid CpG per PMD. Median(%s)" % np.median(values))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_valid_cpg_per_pmd_all.png")

    lower_values = [i for i in values if i <= 5000]
    _ = plt.hist(lower_values)
    plt.style.use('ggplot')
    plt.title("Hist of valid CpG per PMD less than 5000 seq. Median (%s)" % np.median(lower_values))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_valid_cpg_per_pmd_less_than_5000_seq.png")


def plot_methylation_vs_covariance(df):
    cov_columns = df[["cov01", "cov11", "cov13"]]
    methylation_columns = df[["meth01", "meth11", "meth13"]]
    cov_mean = np.mean(cov_columns, axis=1)
    meth_mean = np.mean(methylation_columns, axis=1)

    # plt.plot(meth_mean, cov_mean, linestyle='', marker='o', markersize=0.5)
    # plt.title("Methylation level vs Covariance in PMD")
    # plt.xlabel("Avg methylation level")
    # plt.ylabel("Covariance in window")
    # # plt.savefig("methylation_vs_covariance")
    # plt.show()

    v = np.vstack((meth_mean.values, cov_mean.values))
    dfs = pd.DataFrame(v.T, columns=["meth", "cov"])
    # sns.jointplot(x="meth", y="cov", data=dfs, kind="kde")


def main():
    args = parse_input()
    df = pd.read_pickle(args.cpg_file)

    # global_cpg_info(df)
    # histogram_on_numb_of_cpg_per_pmd(df)
    # histogram_of_coverage_across_patients(df)
    # histogram_of_num_of_cg(df)
    plot_methylation_vs_covariance(df)


if __name__ == '__main__':
    main()
