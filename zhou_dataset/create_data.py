import argparse
import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')
sns.set_style('whitegrid')

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

OUTPUT_FILE = "test.pkl"
TRANSLATION_TABLE = {84: 65, 65: 84, 67: 71, 71: 67}

COV = 0
VAR = 1


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input_file', help='path for the input file', required=True)

    args = parser.parse_args()
    return args


def print_statistics(train, test):
    train_size = train.shape[0]
    test_size = test.shape[0]
    total_size = train_size + test_size
    print("Total data: %s. train data: %s(%s%%). test data: %s(%s%%)"
          % (total_size, train_size, train_size / total_size * 100,
             test_size, test_size / total_size * 100))


def get_reverse_seq(df):
    sequences = df["sequence"]
    flipped = [s.translate(TRANSLATION_TABLE)[::-1] for s in sequences]
    return flipped


def split_df_by_pmd(df):
    total_cpg = df.shape[0] * 1.6
    max_pmd = df["pmd_index"].max()

    upper_limit = int(max_pmd) - 1
    lower_limit = 1
    num_of_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    while num_of_cpg / total_cpg <= 0.1:
        upper_limit -= 5
        lower_limit += 1
        num_of_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    while num_of_cpg / total_cpg > 0.1:
        lower_limit -= 1
        num_of_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    lower_limit += 1
    while num_of_cpg / total_cpg > 0.2:
        upper_limit += 1
        num_of_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    # Fixing the bound due to multiple patients and different in methylations
    lower_limit -= 2
    upper_limit += 2
    print("Going to use all pmd index: %s-%s" % (lower_limit, upper_limit))

    train_part = df[np.logical_and(df["pmd_index"] > lower_limit, df["pmd_index"] < upper_limit)]
    test_part = df[np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit)]

    return train_part, test_part


def add_reverse_strand(df):
    reversed_seq = get_reverse_seq(df)
    df_reverse = df.copy()
    df_reverse["sequence"] = reversed_seq
    df["data_type"] = "original"
    df_reverse["data_type"] = "reverse"

    combined_df = pd.concat([df, df_reverse])

    combined_df["seq4"] = combined_df["sequence"].str[73:77]
    combined_df["seq10"] = combined_df["sequence"].str[70:80]
    return combined_df


def flat_pd_v3(df, patients):
    l = []

    for patient in patients:
        train_df = pd.DataFrame()
        var_label = "var%s" % patient[-2:]
        meth_label = "meth%s" % patient[-2:]

        train_df["meth"] = df[meth_label]
        train_df["var"] = df[var_label]
        train_df["sequence"] = df["sequence"]
        train_df["ccpg"] = df["ccpg"]

        l.append(train_df)

    return pd.concat(l)


def label_based_on_extreme_v3(df):
    df.loc[np.logical_and(df["meth"] >= 0.75, df["var"] <= 0.06), "label"] = 0  # Group 0 - partial loss
    df.loc[np.logical_and(df["meth"] <= 0.55, df["var"] >= 0.1), "label"] = 1  # Group 1 - comp loss

    filtered_df = df[~pd.isnull(df["label"])]
    return filtered_df


def get_high_low_var(df):
    var = df["meth"].values
    var.sort()
    # We want around 800K at the end, so it's 400K seq which mean 200k to each type
    low_var = var[200000]
    high_var = var[len(var) - 200000]
    return low_var, high_var


def create_data():
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)
    df["start"] = df["location"]
    df["ccpg"] = df["sequence"].str.count("CG")

    # We start with solo which are methylated in NC
    df = df[df["ccpg"] == 1]

    # low_var, high_var = get_high_low_var(df)

    # Split the data to train and test based on pmd
    train, test = split_df_by_pmd(df)

    # Leave only the extreme
    train = label_based_on_extreme_v3(train)
    test = label_based_on_extreme_v3(test)

    # Add reverse strand
    train = add_reverse_strand(train)
    test = add_reverse_strand(test)

    print_statistics(train, test)

    output = {"train": train, "test": test}
    with open(os.path.join(args.output_folder, OUTPUT_FILE), "wb") as output_file:
        pickle.dump(output, output_file)

    return train


def plot_data(df):
    plt.subplots_adjust(top=0.9)
    sns.scatterplot(x=df["meth"], y=df["var"], hue=df["label"])

    plt.show()


def main():
    df = create_data()
    pass
    # plot_data(df)
    # df["label"].hist()
    # plt.show()


if __name__ == '__main__':
    main()
