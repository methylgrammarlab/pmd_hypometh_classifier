import argparse
import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tqdm

plt.style.use('seaborn')
sns.set_style('whitegrid')

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

OUTPUT_FILE = "test.pkl"
TRANSLATION_TABLE = {84: 65, 65: 84, 67: 71, 71: 67}

TRAIN_PATIENT = ["CRC01", "CRC11"]
TEST_PATIENT = ["CRC10", "CRC13"]
# TEST_PATIENT = ["CRC13"]

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
    total_cpg = df.shape[0] * 1.8
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


def flat_pd_v2(df, patients):
    l = []

    for patient in patients:
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


def label_based_on_extreme_v3(df, name):
    # Group 1 (partial loss)- variance > 0.2
    # Group 2 (completely loss)- variance < 0.1 with meth < 0.2

    df.loc[np.logical_and(df["meth"] <= 0.2, df["var"] < 0.05), "label"] = 1  # Group 1 - completely loss
    df.loc[df["var"] >= 0.2, "label"] = 0  # Group 0 - partial loss

    filtered_df = df[~pd.isnull(df["label"])]

    print("For %s we keep %s(%s%%) cpg after putting labels" % (name, filtered_df.shape[0],
                                                                filtered_df.shape[0] / df.shape[0] * 100))
    return filtered_df


def flat_and_label_train_based_on_match(df, patients):
    # This is a bad written code

    if len(patients) > 2:
        sys.exit("this code only work for 2 training patient")

    # var_labels = ["var%s" % patient[-2:] for patient in patients]
    #
    # not_null = np.mean(~pd.isnull(df[var_labels]), axis=1) == 1
    # print("We have %s(%s%%) shared cpg between %s" %
    #       (np.sum(not_null), np.sum(not_null) * 100 / df.shape[0], ",".join(patients)))

    # first_filter = df[not_null]
    first_filter = df

    for patient in patients:
        label = "label%s" % patient[-2:]
        var_label = "var%s" % patient[-2:]
        meth_label = "meth%s" % patient[-2:]
        first_filter[label] = 3  # misc
        first_filter.loc[first_filter[meth_label].isnull(), label] = 2  # empty

        # Group 1 - completely los
        first_filter.loc[np.logical_and(first_filter[meth_label] <= 0.2,
                                        first_filter[var_label] <= 0.05), label] = 1

        # Group 0 - partial loss
        first_filter.loc[first_filter[var_label] >= 0.2, label] = 0

    label1 = "label%s" % patients[0][-2:]
    label2 = "label%s" % patients[1][-2:]

    partial_partial = np.logical_and(df[label1] == 0, df[label2] == 0)
    total_total = np.logical_and(df[label1] == 1, df[label2] == 1)
    partial_empty = np.logical_and(df[label1] == 0, df[label2] == 2)
    total_empty = np.logical_and(df[label1] == 1, df[label2] == 2)
    empty_total = np.logical_and(df[label1] == 2, df[label2] == 0)
    empty_partial = np.logical_and(df[label1] == 2, df[label2] == 1)
    good_index = partial_partial | total_total | partial_empty | total_empty | empty_total | empty_partial

    # other_other = np.logical_and(df[label1] == 2, df[label2] == 2)
    # partial_total = np.logical_and(df[label1] == 0, df[label2] == 1)
    # total_partial = np.logical_and(df[label1] == 1, df[label2] == 0)
    # bad_indexes = np.logical_or(other_other, partial_total)
    # bad_indexes = np.logical_or(bad_indexes, total_partial)
    #
    # second_filter = first_filter[~bad_indexes]
    second_filter = first_filter[good_index]

    l = []

    for patient in patients:
        train_df = pd.DataFrame()
        label = "label%s" % patient[-2:]
        var_label = "var%s" % patient[-2:]
        meth_label = "meth%s" % patient[-2:]

        train_df["meth"] = second_filter[meth_label]
        train_df["var"] = second_filter[var_label]
        train_df["sequence"] = second_filter["sequence"]
        train_df["ccpg"] = second_filter["ccpg"]
        train_df["label"] = second_filter[label]
        train_df["chr"] = second_filter["chromosome"]
        train_df["start"] = second_filter["location"]

        l.append(train_df)

    conc = pd.concat(l)
    good_index = np.logical_and(conc["label"] != 2, conc["label"] != 3)
    return conc[good_index]


def v3_variance():
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)
    df["ccpg"] = df["sequence"].str.count("CG")

    # We start with solo which are methylated in NC
    df = df[df["ccpg"] == 1]
    df = df[df["nc_avg"] > 0.5]

    # Split the data to train and test based on pmd
    train, test = split_df_by_pmd(df)
    # train = flat_pd_v3(train, TRAIN_PATIENT)
    # test = flat_pd_v3(test, TEST_PATIENT)

    # Leave only the extreme
    # train = label_based_on_extreme_v3(train, "train")
    # test = label_based_on_extreme_v3(test, "test")

    train = flat_and_label_train_based_on_match(train, TRAIN_PATIENT)
    test = flat_and_label_train_based_on_match(test, TEST_PATIENT)

    # Add reverse strand
    train = add_reverse_strand(train)
    test = add_reverse_strand(test)

    print_statistics(train, test)

    output = {"train": train, "test": test}
    with open(os.path.join(args.output_folder, OUTPUT_FILE), "wb") as output_file:
        pickle.dump(output, output_file)

    return train


def calculate_diff():
    # This is very hard codded stuff
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)
    df["ccpg"] = df["sequence"].str.count("CG")

    # We start with solo which are methylated in NC
    df = df[df["ccpg"] == 1]
    df = df[df["nc_avg"] > 0.5]

    patients = ["11", "01", "13", "10"]

    meth_labels = ["meth%s" % i for i in patients]
    meth_col = df[meth_labels]
    not_null_meth = ~pd.isnull(meth_col)
    not_null_in_at_least_both = not_null_meth.sum(axis=1) > 1

    labels = []
    for patient in patients:
        label = "label_%s" % patient
        labels.append(label)
        df[label] = 2  # else
        df.loc[pd.isnull(df["meth%s" % patient]), label] = -1  # no information
        df.loc[np.logical_and(df["meth%s" % patient] <= 0.2, df["var%s" % patient] < 0.05), label] = 1  # cl
        df.loc[df["var%s" % patient] >= 0.2, label] = 0  # pl

    df = df[not_null_in_at_least_both]  # at least in 2

    num_of_agree = 0
    total = df.shape[0]
    # we can probably do it with numpy but it's hard to think about it
    for i in tqdm.tqdm(range(total)):
        p1 = df.iloc[i][labels[0]]
        p2 = df.iloc[i][labels[1]]
        p3 = df.iloc[i][labels[2]]
        p4 = df.iloc[i][labels[3]]

        agree = True
        use = [i for i in [p1, p2, p3, p4] if i != -1]

        last_value = use[0]
        for i in use[1:]:
            if last_value != i:
                agree = False

        if agree:
            num_of_agree += 1

    print(num_of_agree / total * 100)


def plot_data(df):
    plt.subplots_adjust(top=0.9)
    sns.scatterplot(x=df["meth"], y=df["var"], hue=df["label"])

    plt.show()


def main():
    df = v3_variance()
    plot_data(df)
    df["label"].hist()
    plt.show()
    # calculate_diff()


if __name__ == '__main__':
    main()
