import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')
sns.set_style('whitegrid')


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input_file', help='path for the input file', required=True)
    parser.add_argument('--add_perc', help='path for the input file', required=False, type=bool,
                        default=True)

    args = parser.parse_args()
    return args


def main():
    args = parse_input()
    input_path = args.input_file
    add_perc = args.add_perc
    full_df = pd.read_csv(input_path)

    label1 = "meth"
    step1 = 0.01

    label2 = "coveriance"
    step2 = 0.001

    df = full_df[full_df["ccpg"] < 4]
    df = df[df["orig_meth"] >= 0.7]
    solo_wcgw = np.logical_and(df["ccpg"] == 1, df["small_seq"].str.contains("[AT]CG[AT]", regex=True))
    df["solo-WCGW"] = solo_wcgw
    df["meth_diff"] = df["orig_meth"] - df["meth"]

    x = np.arange(df[label1].min(), df[label1].max(), step1)
    y = np.arange(float("%.2f" % df[label2].min()), float("%.2f" % df[label2].max()), step2)

    data_matrix = np.zeros(shape=(len(x), len(y)))
    size_matrix = np.zeros(shape=(len(x), len(y)))

    row = 0

    for i in x:
        column = len(y) - 1
        for j in y:
            indexes = np.logical_and(np.logical_and(df[label1] >= i, df[label1] < i + step1),
                                     np.logical_and(df[label2] >= j, df[label2] < j + step2))
            temp_df = df[indexes]
            num_of_cpg = temp_df.shape[0]
            perc_of_solo_wcgw = np.sum(temp_df["solo-WCGW"]) / num_of_cpg * 100 if num_of_cpg > 500 \
                else np.nan

            data_matrix[row, column] = perc_of_solo_wcgw
            size_matrix[row, column] = num_of_cpg
            column -= 1

            if add_perc:
                df.loc[indexes, "soloWCGW_perc_%s_%s" % (label1, label2)] = perc_of_solo_wcgw

        row += 1

    x = [float("%.2f" % i) for i in x]
    y_reverse = [float("%.2f" % i) for i in y]
    y_reverse.reverse()

    df_data_matrix = pd.DataFrame(data_matrix, index=x, columns=y_reverse)
    df_size_matrix = pd.DataFrame(size_matrix, index=x, columns=y_reverse)

    df_data_matrix.to_pickle(os.path.join(args.output_folder, "%s_%s_matrix.pkl" % (label1, label2)))
    df_size_matrix.to_pickle(os.path.join(args.output_folder, "%s_%s_size.pkl" % (label1, label2)))

    plt.subplots_adjust(top=0.9)
    sns.heatmap(df_data_matrix)
    plt.title("%s vs %s for %%CpG1" % (label1, label2), fontsize=20)
    plt.xlabel("%s" % label2, fontsize=16)
    plt.ylabel("%s" % label1, fontsize=16)

    plt.savefig(os.path.join(args.output_folder, "%s_vs_%s.png" % (label1, label2)))
    plt.close()

    plt.subplots_adjust(top=0.9)
    sns.heatmap(df_size_matrix)
    plt.title("%s vs %s total CpGs" % (label1, label2), fontsize=20)
    plt.xlabel("%s" % label2, fontsize=16)
    plt.ylabel("%s" % label1, fontsize=16)
    plt.savefig(os.path.join(args.output_folder, "%s_vs_%s_total_cpg.png" % (label1, label2)))
    plt.close()

    if add_perc:
        file_folder = os.path.dirname(input_path)
        file_name = os.path.basename(input_path).split(".")[0] + "with_perc" + ".pkl"
        df.to_pickle(os.path.join(file_folder, file_name))


if __name__ == '__main__':
    main()
