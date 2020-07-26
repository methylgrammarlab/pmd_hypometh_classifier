import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')
sns.set_style('whitegrid')


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


def plot_variance_histogram_vs_type(df, output):
    label = "var"
    meth_label = "meth"

    strong_rows = df[df["sequence"].str.contains("[CG]CG[CG]", regex=True)]
    weak_rows = df[df["sequence"].str.contains("[AT]CG[AT]", regex=True)]
    _ = plt.hist(weak_rows[label], label="WCGW", normed=True, alpha=0.5, color="b")
    _ = plt.hist(strong_rows[label], label="SCGS", normed=True, alpha=0.5, color="g")
    plt.legend()
    plt.style.use('ggplot')
    plt.title("Histogram of variance in strong vs weak")
    plt.xlabel("Variance")
    plt.ylabel("Amount")
    plt.savefig(os.path.join(output, "var_strong_vs_weak.png"))
    plt.close()


def plot_meth_density_violin(df, output):
    meth_label = "meth"
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
    plt.title("Methylation density")
    plt.savefig(os.path.join(output, "density_meth.png"))
    plt.close()


def plot_meth_vs_var_scatter_dens(df, output):
    var_label = "var"
    meth_label = "meth"

    plt.subplots_adjust(top=0.9)
    sns_plot = sns.jointplot(x=df[meth_label], y=df[var_label], kind="kde")
    sns_plot.fig.suptitle("Methylation vs Variance in solo CpG", fontsize=20)
    sns_plot.set_axis_labels("Methylation level", "Variance", fontsize=16)

    sns_plot.savefig(os.path.join(output, "dist_solo_meth_vs_var_scatter.png"))
    plt.close()


def plot_meth_vs_var_scatter(df, output):
    var_label = "var"
    meth_label = "meth"

    _ = plt.scatter(x=df[meth_label], y=df[var_label], marker=".")
    plt.xlim(left=0)
    plt.xlabel("Methylation", fontsize=16)
    plt.ylabel("Variance", fontsize=16)
    plt.title("Methylation vs Variance in solo CpG")

    plt.savefig(os.path.join(output, "meth_vs_var_scatter.png"))
    plt.close()


def main():
    args = parse_input()
    df = pd.read_pickle(args.cpg_file)
    df["small_seq"] = df["sequence"].str[73:77]

    # keep only solo
    solo_rows = df[df["sequence"].str.count("CG") == 1]
    # global_cpg_info(df)
    # plot_variance_histogram_vs_type(solo_rows, args.output_folder)
    # plot_meth_density_violin(df, args.output_folder)
    # plot_meth_vs_var_scatter_dens(solo_rows, args.output_folder)
    plot_meth_vs_var_scatter(solo_rows, args.output_folder)


if __name__ == '__main__':
    main()
