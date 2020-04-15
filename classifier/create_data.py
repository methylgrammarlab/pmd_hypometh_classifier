import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

OUTPUT_FILE = "classifier_data_ccpg1.pkl"
TRANSLATION_TABLE = {84: 65, 65: 84, 67: 71, 71: 67}

TRAIN_PATIENT = ["CRC11", "CRC01"]
TEST_PATIENT = ["CRC13"]

TRAIN_EXTREME = 1 / 3
TEST_EXTREME = 1 / 4

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


def split_df_by_pmd(df, test_patients=TEST_PATIENT, train_patient=TRAIN_PATIENT):
    total_cpg = df.shape[0] * 2.5
    max_pmd = df["pmd_index"].max()

    upper_limit = int(max_pmd) - 1
    lower_limit = 1
    num_of_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    while num_of_cpg / total_cpg <= 0.2:
        upper_limit -= 5
        lower_limit += 1
        num_of_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    while num_of_cpg / total_cpg > 0.2:
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


def flat_pd_v1(df, patients):
    l = []

    for patient in patients:
        # Remove cpg which weren't methylated to begin with
        nc_meth_label = "nc_meth%s" % patient[-2:]
        patient_train_part = df[df[nc_meth_label] > 0.6]

        train_df = pd.DataFrame()
        meth_label = "meth%s" % patient[-2:]
        train_df["meth"] = patient_train_part[meth_label]
        train_df["sequence"] = patient_train_part["sequence"]
        train_df["ccpg"] = patient_train_part["ccpg"]
        train_df["nc_avg"] = patient_train_part["nc_avg"]

        l.append(train_df)

    return pd.concat(l)


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


def label_based_on_extreme_v1(df, extreme_value):
    sorted_values = df["meth"].sort_values()
    size_of_df = df.shape[0]

    low_index = int(np.round(size_of_df * extreme_value) + 1)
    high_index = int(np.round(size_of_df * (1 - extreme_value)) - 1)

    low_value, high_value = sorted_values.iloc[low_index], sorted_values.iloc[high_index]

    filtered_df = df[np.logical_or(df["meth"] <= low_value, df["meth"] >= high_value)]
    filtered_df.loc[filtered_df["meth"] >= high_value, "label"] = 0  # State didn't change
    filtered_df.loc[filtered_df["meth"] <= low_value, "label"] = 1  # State changed
    return filtered_df


def label_based_on_extreme_v2(df, hv=COV):
    # Group 1 methylation in nc > 0.5 (avg) and methylation in cancer < 0.1
    # group 2 methylation in cancer > 0.5 and cov\variance in top 10%

    if hv == COV:
        sorted_values = df["cov"].sort_values()
        sorted_values = sorted_values[~pd.isnull(sorted_values)]
        high_value = sorted_values.quantile(0.9)
        df.loc[np.logical_and(df["cov"] > high_value, df["meth"] > 0.5), "label"] = 0  # state kept - group 2



    elif hv == VAR:
        sorted_values = df["var"].sort_values()
        sorted_values = sorted_values[~pd.isnull(sorted_values)]
        high_value = sorted_values.quantile(0.9)
        df.loc[np.logical_and(df["var"] > high_value, df["meth"] > 0.5), "label"] = 0  # state kept - group 2

    df.loc[np.logical_and(df["meth"] < 0.1, df["nc_avg"] > 0.5), "label"] = 1  # State changed - group 1

    filtered_df = df[~pd.isnull(df["label"])]
    return filtered_df


def v1():
    # Only use methylation to create the labels
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)
    df["ccpg"] = df["sequence"].str.count("CG")

    # We start with solo
    df = df[df["ccpg"] == 1]

    # Split the data to train and test based on pmd
    train, test = split_df_by_pmd(df)
    train = flat_pd_v1(train, TRAIN_PATIENT)
    test = flat_pd_v1(test, TEST_PATIENT)

    # Leave only the extreme
    train = label_based_on_extreme_v1(train, TRAIN_EXTREME)
    test = label_based_on_extreme_v1(test, TEST_EXTREME)

    # Add reverse strand
    train = add_reverse_strand(train)
    test = add_reverse_strand(test)

    print_statistics(train, test)

    output = {"train": train, "test": test}
    with open(os.path.join(args.output_folder, OUTPUT_FILE), "wb") as output_file:
        pickle.dump(output, output_file)


def flat_pd_v2(df, TRAIN_PATIENT):
    l = []

    for patient in TRAIN_PATIENT:
        train_df = pd.DataFrame()
        meth_label = "meth%s" % patient[-2:]
        cov_label = "cov%s" % patient[-2:]
        var_label = "cancer_var%s" % patient[-2:]

        train_df["meth"] = df[meth_label]
        train_df["var"] = df[var_label]
        train_df["cov"] = df[cov_label]
        train_df["sequence"] = df["sequence"]
        train_df["ccpg"] = df["ccpg"]
        train_df["nc_avg"] = df["nc_avg"]

        l.append(train_df)

    return pd.concat(l)


def v2():
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)
    df["ccpg"] = df["sequence"].str.count("CG")

    # We start with solo
    df = df[df["ccpg"] == 1]

    # Split the data to train and test based on pmd
    train, test = split_df_by_pmd(df)
    train = flat_pd_v2(train, TRAIN_PATIENT)
    test = flat_pd_v2(test, TEST_PATIENT)

    # Leave only the extreme
    train = label_based_on_extreme_v2(train, hv=COV)
    test = label_based_on_extreme_v2(test, hv=COV)

    # Add reverse strand
    train = add_reverse_strand(train)
    test = add_reverse_strand(test)

    print_statistics(train, test)

    output = {"train": train, "test": test}
    with open(os.path.join(args.output_folder, OUTPUT_FILE), "wb") as output_file:
        pickle.dump(output, output_file)


def main():
    v2()

if __name__ == '__main__':
    main()
