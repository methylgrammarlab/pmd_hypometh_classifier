import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

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
    weak_solo = [seq for seq in solo if seq[73] in ["A", "T"] and seq[76] in ["A", "T"]]
    strong_solo = [seq for seq in solo if seq[73] in ["C", "G"] and seq[76] in ["C", "G"]]

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
    plt.close()


def histogram_of_num_of_cg(df):
    num_of_cg = [seq.count("CG") for seq in df["sequence"]]
    _ = plt.hist(num_of_cg)
    plt.style.use('ggplot')
    plt.title("Number of CpG in seq include middle. Median (%s)" % np.median(num_of_cg))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_cg_in_seq_include_middle_all.png")
    plt.close()

    num_of_cg = [i for i in num_of_cg if i < 5]
    _ = plt.hist(num_of_cg)
    plt.style.use('ggplot')
    plt.title("Number of CpG in seq include middle less than 5. Median (%s)" % np.median(num_of_cg))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_cg_in_seq_include_middle_lower_than_5.png")
    plt.close()


def histogram_on_numb_of_cpg_per_pmd(df):
    values = df.groupby(["chromosome", "pmd_index"]).count()["sequence"].values
    values.sort()
    _ = plt.hist(values)
    plt.style.use('ggplot')
    plt.title("Hist of valid CpG per PMD. Median(%s)" % np.median(values))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_valid_cpg_per_pmd_all.png")
    plt.close()

    lower_values = [i for i in values if i <= 5000]
    _ = plt.hist(lower_values)
    plt.style.use('ggplot')
    plt.title("Hist of valid CpG per PMD less than 5000 seq. Median (%s)" % np.median(lower_values))
    plt.xlabel("Number of CpG")
    plt.ylabel("Amount")
    plt.savefig("num_of_valid_cpg_per_pmd_less_than_5000_seq.png")
    plt.close()


def plot_methylation_vs_covariance(df):
    meth = df["meth_mean"]
    cov = df["cov_mean"]
    outliers_ind = np.abs(stats.zscore(cov)) < 3

    df_i = df[outliers_ind]
    meth_i = meth[outliers_ind]
    cov_i = cov[outliers_ind]

    # Draw 2d plot
    h = plt.hist2d(meth_i, cov_i)
    plt.colorbar(h[3])
    plt.savefig("2dhist_methylation_vs_covariance_global.png")
    plt.close()

    # Draw only PMD
    ind = df_i["pmd_index"] > 0
    meth_ip = meth_i[ind]
    cov_ip = cov_i[ind]

    plt.plot(meth_ip, cov_ip, linestyle='', marker='o', markersize=0.5)
    plt.title("Methylation level vs Covariance in PMD")
    plt.xlabel("Avg methylation level")
    plt.ylabel("Covariance in window")
    plt.savefig("methylation_vs_covariance_pmds.png")
    plt.close()

    v = np.vstack((meth_ip.values, cov_ip.values))
    dfs = pd.DataFrame(v.T, columns=["meth", "cov"])
    sns_plot = sns.jointplot(x="meth", y="cov", data=dfs, kind="kde")
    sns_plot.savefig("methylation_vs_covariance_pmds_cluster.png")
    plt.close()


def plot_methylation_vs_covariance_solo(df):
    solo_rows = df[df["sequence"].str.count("CG") == 1]
    solo_rows = solo_rows[solo_rows["pmd_index"] <= 5]
    meth_mean = solo_rows["meth_mean"]
    cov_mean = solo_rows["cov_mean"]

    # plt.plot(meth_mean, cov_mean, linestyle='', marker='o', markersize=0.5)
    # plt.title("Methylation level vs Covariance in solo CpG")
    # plt.xlabel("Avg methylation level")
    # plt.ylabel("Covariance in window")
    # plt.savefig("solo_meth_vs_cov_scatter.png")
    # plt.close()
    #
    # v = np.vstack((meth_mean.values, cov_mean.values))
    # dfs = pd.DataFrame(v.T, columns=["meth", "cov"])
    # sns_plot = sns.jointplot(x="meth", y="cov", data=dfs, kind="kde")
    # sns_plot.savefig("solo_meth_vs_cov.png")
    # plt.close()

    # meth_mean = np.round(meth_mean * 500).astype(np.int) / 500
    z_cov = np.abs(stats.zscore(cov_mean.values))
    outliers_ind = z_cov < 3
    plt.plot(meth_mean[outliers_ind], cov_mean[outliers_ind], linestyle='', marker='o', markersize=0.5)
    plt.title("Methylation level vs Covariance in solo CpG")
    plt.xlabel("Avg methylation level")
    plt.ylabel("Covariance in window")
    plt.savefig("solo_meth_vs_cov_scatter.png")
    plt.close()

    v = np.vstack((meth_mean[outliers_ind].values, cov_mean[outliers_ind].values))
    dfs = pd.DataFrame(v.T, columns=["meth", "cov"])
    sns_plot = sns.jointplot(x="meth", y="cov", data=dfs, kind="kde")
    sns_plot.savefig("solo_meth_vs_cov.png")
    plt.close()


def check_flip_seq(df):
    sequences = df["sequence"]
    seq_uniq = set(sequences)
    print("We have %s uniq seq to beging with" % len(seq_uniq))

    translation_table = {84: 65, 65: 84, 67: 71, 71: 67}
    flipped = []
    for s in seq_uniq:
        seq_translated = s.translate(translation_table)
        seq_flipped = seq_translated[::-1]
        flipped.append(seq_flipped)

    flipped_uniq = set(flipped)

    comb = flipped_uniq | seq_uniq
    print("Combined: %s" % len(comb))


def plot_methylation_vs_covariance_solo_vs_non_solo(df):
    pmd_df = df[df["pmd_index"] <= 10]
    solo_rows = pmd_df[pmd_df["sequence"].str.count("CG") == 1]
    non_solo_rows = pmd_df[pmd_df["sequence"].str.count("CG") > 1]

    solo_meth = solo_rows["meth_mean"]
    solo_cov = solo_rows["cov_mean"]
    nsolo_meth = non_solo_rows["meth_mean"]
    nsolo_cov = non_solo_rows["cov_mean"]

    solo_z_cov = np.abs(stats.zscore(solo_cov.values))
    solo_outliers_ind = solo_z_cov < 3

    nsolo_z_cov = np.abs(stats.zscore(nsolo_cov.values))
    nsolo_outliers_ind = nsolo_z_cov < 3

    solo = plt.scatter(solo_meth[solo_outliers_ind], solo_cov[solo_outliers_ind], marker='x', color='r')
    nsolo = plt.scatter(nsolo_meth[nsolo_outliers_ind], nsolo_cov[nsolo_outliers_ind], marker='o',
                        color='b', alpha=0.05)
    plt.title("Methylation level vs Covariance in solo(red) vs non-solo (blue)")
    plt.xlabel("Avg methylation level")
    plt.ylabel("Covariance in window")
    plt.savefig("methylation_vs_covariance_solo_vs_non_solo.png")
    plt.close()


def plot_methylation_vs_covariance_solo_weak_vs_strong(df):
    pmd_df = df[df["pmd_index"] <= 10]
    solo_rows = pmd_df[pmd_df["sequence"].str.count("CG") == 1]
    ou = np.abs(stats.zscore(solo_rows["cov_mean"].values)) < 3
    solo_rows = solo_rows[ou]

    strong_rows = solo_rows[solo_rows["sequence"].str.contains("[CG]CG[CG]", regex=True)]
    weak_rows = solo_rows[solo_rows["sequence"].str.contains("[AT]CG[AT]", regex=True)]

    weak_meth = weak_rows["meth_mean"]
    weak_cov = weak_rows["cov_mean"]
    weak_z_cov = np.abs(stats.zscore(weak_cov.values))
    weak_outliers_ind = weak_z_cov < 3

    weak_meth_i = weak_meth
    weak_cov_i = weak_cov

    strong_meth = strong_rows["meth_mean"]
    strong_cov = strong_rows["cov_mean"]
    strong_z_cov = np.abs(stats.zscore(strong_cov.values))
    strong_outliers_ind = strong_z_cov < 3

    strong_meth_i = strong_meth
    strong_cov_i = strong_cov

    solo = plt.scatter(weak_meth, weak_cov, marker='x', color='r')
    nsolo = plt.scatter(strong_meth, strong_cov, marker='o', color='b', alpha=0.05)
    plt.title("Methylation level vs Covariance in WSCW(red) vs SCGS(blue)")
    plt.xlabel("Avg methylation level")
    plt.ylabel("Covariance in window")
    plt.savefig("methylation_vs_covariance_strong_vs_weak.png")
    plt.close()


def plot_box_violin_all(df):
    pmd_df = df[np.logical_and(df["meth01"] >= 0.4, df["meth01"] <= 0.6)]  # not really, just middle lines
    # pmd_df = pmd_df [~pmd_df ["cov01"].isnull()]
    x = []

    solo_rows = pmd_df[pmd_df["sequence"].str.count("CG") == 1]
    cov_solo = solo_rows["cov01"]
    # x.extend(["solo"] * len(cov_solo))

    nsolo_rows = pmd_df[pmd_df["sequence"].str.count("CG") > 1]
    cov_nsolo = nsolo_rows["cov01"]
    # x.extend(["not-solo"] * len(cov_nsolo))

    weak_rows = solo_rows[solo_rows["sequence"].str.contains("[AT]CG[AT]", regex=True)]
    strong_rows = solo_rows[solo_rows["sequence"].str.contains("[CG]CG[CG]", regex=True)]

    cov_weak = weak_rows["cov01"]
    # x.extend(["solo-WCGW"] * len(cov_weak))
    cov_strong = strong_rows["cov01"]
    # x.extend(["solo-SCGS"] * len(cov_strong))
    #
    # data = np.append(cov_solo.values, cov_nsolo.values)
    # data = np.append(data, cov_weak.values)
    # data = np.append(data, cov_strong.values)
    #
    #
    # _ = plt.hist(cov_solo, bins=40)
    # plt.style.use('ggplot')
    # plt.title("Covariance of solo CpG")
    # # plt.xlabel("Number of CpG")
    # # plt.ylabel("Amount")
    # plt.savefig("solo_cpg_cov.png")
    # plt.close()
    #
    # _ = plt.hist(cov_nsolo,bins=40)
    # plt.style.use('ggplot')
    # plt.title("Covariance of non-solo cpv")
    # # plt.xlabel("Number of CpG")
    # # plt.ylabel("Amount")
    # plt.savefig("not_solo_cpg_cov.png")
    # plt.close()
    #
    # _ = plt.hist(cov_strong, bins=40)
    # plt.style.use('ggplot')
    # plt.title("Covariance of strong CpG")
    # # plt.xlabel("Number of CpG")
    # # plt.ylabel("Amount")
    # plt.savefig("scgs_cov.png")
    # plt.close()
    #
    # _ = plt.hist(cov_weak, bins=40)
    # plt.style.use('ggplot')
    # plt.title("Covariance of weak CpG")
    # # plt.xlabel("Number of CpG")
    # # plt.ylabel("Amount")
    # plt.savefig("wcgw_cov.png")
    # plt.close()

    # colors = ['#E69F00', '#56B4E9', '#F0E442', '#009E73']
    # names = ['not solo', 'solo' , 'strong', 'weak']
    # plt.hist([cov_nsolo, cov_solo, cov_strong, cov_weak], bins = int(180/15), normed=True,
    #          color = colors, label=names)
    # plt.legend()
    # plt.show()

    sns.distplot(cov_nsolo, hist=False, kde=True, kde_kws={'linewidth': 3}, label="not solo")
    sns.distplot(cov_solo, hist=False, kde=True, kde_kws={'linewidth': 3}, label="solo")
    sns.distplot(cov_strong, hist=False, kde=True, kde_kws={'linewidth': 3}, label="SCGS")
    sns.distplot(cov_weak, hist=False, kde=True, kde_kws={'linewidth': 3}, label="WSCW")
    plt.show()

    #
    # kwargs = dict(histtype='stepfilled', alpha=0.3, normed=True, bins=40)
    # _ = plt.hist(cov_solo, label="solo", **kwargs)
    # _ = plt.hist(cov_nsolo, label="not solo",**kwargs)
    # _ = plt.hist(cov_strong, label="strong",**kwargs)
    # _ = plt.hist(cov_weak, label="weak",**kwargs)
    # plt.style.use('ggplot')
    # plt.legend(loc="upper right")
    # plt.show()


def plot_box_violin_all_meth(df, meth_v):
    # pmd_df = df[np.logical_and(df["meth01"] >=0.4, df["meth01"] <= 0.6)]# not really, just middle lines
    pmd_df = df[~df[meth_v].isnull()]
    x = []

    solo_rows = pmd_df[pmd_df["sequence"].str.count("CG") == 1]
    cov_solo = solo_rows[meth_v]
    # x.extend(["solo"] * len(cov_solo))

    nsolo_rows = pmd_df[pmd_df["sequence"].str.count("CG") > 1]
    cov_nsolo = nsolo_rows[meth_v]

    weak_rows = solo_rows[solo_rows["sequence"].str.contains("[AT]CG[AT]", regex=True)]
    strong_rows = solo_rows[solo_rows["sequence"].str.contains("[CG]CG[CG]", regex=True)]

    cov_weak = weak_rows[meth_v]
    cov_strong = strong_rows[meth_v]

    sns.distplot(cov_nsolo, hist=False, kde=True, kde_kws={'linewidth': 3}, label="not solo")
    sns.distplot(cov_solo, hist=False, kde=True, kde_kws={'linewidth': 3}, label="solo")
    sns.distplot(cov_strong, hist=False, kde=True, kde_kws={'linewidth': 3}, label="SCGS")
    sns.distplot(cov_weak, hist=False, kde=True, kde_kws={'linewidth': 3}, label="WSCW")
    if meth_v != "meth_mean":
        plt.title("Methylation density for CRC%s" % meth_v[-2:])
    else:
        plt.title("Methylation density - avg")

    plt.savefig("methylation_dens_for_%s.png" % (meth_v))

    plt.close()


def plot_patients_meth(df):
    solo_rows = df[df["sequence"].str.count("CG") == 1]

    cov1 = solo_rows["meth13"]
    cov11 = solo_rows["meth11"]
    # plt.scatter(cov1.values, cov11.values)
    # plt.show()
    v = np.vstack((cov1.values, cov11.values))
    dfs = pd.DataFrame(v.T, columns=["meth13", "meth11"])
    # sns_plot = sns.jointplot(x="cov1", y="cov11", data=dfs, kind="kde")
    # sns_plot.savefig("cov1vs11.png")
    # plt.close()

    sample = dfs.sample(frac=0.005)
    plt.scatter(sample["meth13"], sample["meth11"])
    plt.xlabel("meth13")
    plt.ylabel("meth11")
    plt.title("sample of meth01 and meth11 methylation")
    plt.show()
    plt.close()
    #
    # avg = (cov1 + cov11) /2
    # diff = cov1 - cov11
    # v = np.vstack((avg.values, diff.values))
    # dfs = pd.DataFrame(v.T, columns=["avg", "diff"])
    # sns_plot = sns.jointplot(x="avg", y="diff", data=dfs, kind="kde")
    # sns_plot.savefig("avg_meth_vs_diff.png")
    # plt.close()


def get_num_of_seq_in_extreme_meth(df):
    # Only create solo 
    solo_rows = df[df["sequence"].str.count("CG") == 1]
    methylation_columns = solo_rows[["meth01", "meth11", "meth13"]]

    solo_rows["min_meth"] = np.min(methylation_columns, axis=1)
    solo_rows["max_meth"] = np.max(methylation_columns, axis=1)

    # v = pd.notnull(methylation_columns).sum(axis=1).values
    # full_info_rows = solo_rows[v == 3]

    extreme_index = np.logical_or(solo_rows["max_meth"] <= 0.2, solo_rows["min_meth"] >= 0.8)
    extreme_rows = solo_rows[extreme_index]
    print("We will be using %s/%s of seq" % (extreme_rows.shape[0], solo_rows.shape[0]))


def plot_densitiy_met_cov_for_cpg_numb(df):
    # df = df[df["chromosome"] =="1"]
    for density in range(1, 7):
        pmd_df = df[~df["meth01"].isnull()]

        rows = pmd_df[pmd_df["sequence"].str.count("CG") == density]
        meth_values = rows["meth01"]

        weak_rows = rows[rows["sequence"].str.contains("[AT]CG[AT]", regex=True)]
        strong_rows = rows[rows["sequence"].str.contains("[CG]CG[CG]", regex=True)]

        meth_weak = weak_rows["meth01"]
        meth_strong = strong_rows["meth01"]

        sns.distplot(meth_values, hist=False, kde=True, kde_kws={'linewidth': 3}, label="any")
        sns.distplot(meth_strong, hist=False, kde=True, kde_kws={'linewidth': 3}, label="SCGS")
        sns.distplot(meth_weak, hist=False, kde=True, kde_kws={'linewidth': 3}, label="WSCW")
        plt.title("Methylation density for CRC%swith cpg=%s" % ("meth01"[-2:], density))
        plt.savefig("methylation_dens_for_%s_cpg%s.png" % ("meth01", density))

        plt.close()

        pmd_df = pmd_df[~pmd_df["cov01"].isnull()]

        rows = pmd_df[pmd_df["sequence"].str.count("CG") == density]
        cov_values = rows["cov01"]

        weak_rows = rows[rows["sequence"].str.contains("[AT]CG[AT]", regex=True)]
        strong_rows = rows[rows["sequence"].str.contains("[CG]CG[CG]", regex=True)]

        cov_weak = weak_rows["cov01"]
        cov_strong = strong_rows["cov01"]
        sns.distplot(cov_values, hist=False, kde=True, kde_kws={'linewidth': 3}, label="any")
        sns.distplot(cov_strong, hist=False, kde=True, kde_kws={'linewidth': 3}, label="SCGS")
        sns.distplot(cov_weak, hist=False, kde=True, kde_kws={'linewidth': 3}, label="WSCW")
        plt.title("Cov density for CRC%swith cpg=%s" % ("cov01"[-2:], density))
        plt.savefig("cov_dens_for_%s_cpg%s.png" % ("cov01", density))

        plt.close()


def plot_densitiy_cov_between_02_to_06(df):
    meth_v = "meth01"
    pmd_df = df[np.logical_and(df["meth01"] >= 0.2, df["meth01"] <= 0.6)]  # not really, just middle lines

    solo_rows = pmd_df[pmd_df["sequence"].str.count("CG") == 1]
    cov_solo = solo_rows["cov01"]

    nsolo_rows = pmd_df[pmd_df["sequence"].str.count("CG") > 1]
    cov_nsolo = nsolo_rows["cov01"]

    weak_rows = solo_rows[solo_rows["sequence"].str.contains("[AT]CG[AT]", regex=True)]
    strong_rows = solo_rows[solo_rows["sequence"].str.contains("[CG]CG[CG]", regex=True)]

    cov_weak = weak_rows["cov01"]
    cov_strong = strong_rows["cov01"]

    sns.distplot(cov_nsolo[np.abs(stats.zscore(cov_nsolo.values)) < 3], hist=False, kde=True,
                 kde_kws={'linewidth': 3}, label="not solo")
    sns.distplot(cov_solo[np.abs(stats.zscore(cov_solo.values)) < 3], hist=False, kde=True,
                 kde_kws={'linewidth': 3}, label="solo")
    sns.distplot(cov_strong[np.abs(stats.zscore(cov_strong.values)) < 3], hist=False, kde=True,
                 kde_kws={'linewidth': 3}, label="SCGS")
    sns.distplot(cov_weak[np.abs(stats.zscore(cov_weak.values)) < 3], hist=False, kde=True,
                 kde_kws={'linewidth': 3}, label="WSCW")
    if meth_v != "meth_mean":
        plt.title("Methylation density for CRC%s" % meth_v[-2:])
    else:
        plt.title("Methylation density - avg")

    plt.savefig("methylation_dens_for_%s_02_to_06.png" % (meth_v))

    plt.close()

def main():
    args = parse_input()
    df = pd.read_pickle(args.cpg_file)
    cov_columns = df[["cov01"]]
    methylation_columns = df[["meth01"]]

    cov_mean = np.mean(cov_columns, axis=1)
    meth_mean = np.mean(methylation_columns, axis=1)

    df["cov_mean"] = cov_mean
    df["meth_mean"] = meth_mean

    # global_cpg_info(df)
    # histogram_on_numb_of_cpg_per_pmd(df)
    # histogram_of_coverage_across_patients(df)
    # histogram_of_num_of_cg(df)
    plot_methylation_vs_covariance(df)
    # plot_methylation_vs_covariance_solo(df)
    # check_flip_seq(df)
    # plot_methylation_vs_covariance_solo_vs_non_solo(df)
    # plot_methylation_vs_covariance_solo_weak_vs_strong(df)
    # plot_box_violin_all(df)
    # plot_box_violin_all_meth(df, "meth01")
    # plot_box_violin_all_meth(df, "meth11")
    # plot_box_violin_all_meth(df, "meth13")
    # plot_box_violin_all_meth(df, "meth_mean")
    #
    # for v in df.columns:
    #     if v.startswith("meth"):
    #         plot_box_violin_all_meth(df, v)

    # plot_patients_meth(df)
    # get_num_of_seq_in_extreme_meth(df)
    # plot_densitiy_met_cov_for_cpg_numb(df)
    # plot_densitiy_cov_between_02_to_06(df)

if __name__ == '__main__':
    main()
