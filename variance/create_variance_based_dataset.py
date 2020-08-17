"""
Create dataset for the nn, adding labels, splitting to test and train in a way that will consider the pmd
and add reverse strand.
This code work both for scwgbs and bulk
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools, consts, sequence_tools

# scWGBS
TRAIN_PATIENT = ["CRC01", "CRC11"]
TEST_PATIENT = ["CRC10", "CRC13"]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input_file', help='path for the input file', required=True)
    parser.add_argument('--parse_format', help='do we parse scwgb or bulk data', required=True, type=str,
                        choices=consts.PARSE_FORMAT_OPTIONS)
    parser.add_argument('--output_name', help='name of the output file', required=False,
                        default="nn_dataset.pkl")

    args = parser.parse_args()
    return args


def print_statistics(train, test):
    """
    Print some basic statistics about the train and test dataset
    :param train: The train df
    :param test: The test df
    """
    train_size = train.shape[0]
    test_size = test.shape[0]
    total_size = train_size + test_size
    print("Total data: %s. train data: %s(%s%%). test data: %s(%s%%)"
          % (total_size, train_size, train_size / total_size * 100,
             test_size, test_size / total_size * 100))


def split_df_by_pmd(df, size_wheel=3.0):
    """
    Split data frame based on pmd, we assume we have a column which is called pmd_index
    We don't want to have CpG from the same PMD in test and train - we want to make it dependent + we want
    to have PMD from all chromosome this mean we can't take only the first PMD (10% test data will cause
    that chromosome 16-22 only on test) and can't take only the last PMD(10% test data will cause that
    there is no chr16-22 in test) so we take some of the beginning and some from the end and try to make it
    best depend on the size_wheel which help make sure that we have enough data
    :param df: The dataframe with all the data
    :param size_wheel: This is used to change the splitting method and getting the right percentage of the
    data for test and train, due to the PMD and the fact that each PMD has a different size of CpG it's
    hard to put a number and get it splited exactly so we use this with trail and error until getting the
    right number
    :return: train_part, test_part
    """
    total_cpg = df.shape[0] * size_wheel
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


def double_with_reverse_strand(df):
    """
    Double the data by adding the reverse strand
    :param df: The data frame
    :return: A new data frame with the reverse strand
    """
    reversed_seq = sequence_tools.get_reverse_seq(df["sequence"])
    df_reverse = df.copy()
    df_reverse["sequence"] = reversed_seq
    df["data_type"] = "original"
    df_reverse["data_type"] = "reverse"

    combined_df = pd.concat([df, df_reverse])

    combined_df["seq4"] = combined_df["sequence"].str[73:77]
    combined_df["seq10"] = combined_df["sequence"].str[70:80]
    return combined_df


def label_sc_based_on_meth_var_flat(df, cl_max_meth=0.2, cl_max_var=0.05, pl_min_meth=0, pl_min_var=0.2):
    """
    Label the df based on meth and var values, this is only good for single cell where the cl has low
    methylation and low variance and the pl has high var and high meth
    :param df: The dataframe to label
    :param cl_max_meth: THe max meth for cl
    :param cl_max_var: The max variance for cl
    :param pl_min_meth: The min meth for pl
    :param pl_min_var: The minimum var for pl
    :return: The new dataframe with labels, removing CpG which didn't include
    """
    df.loc[np.logical_and(df["meth"] >= pl_min_meth, df["var"] >= pl_min_var), "label"] = \
        consts.LABEL_HYPO_RESISTANCE
    df.loc[np.logical_and(df["meth"] <= cl_max_meth, df["var"] <= cl_max_var), "label"] = \
        consts.LABEL_HYPO_PRONE

    return df[~pd.isnull(df["label"])]


def label_and_flat_sc_based_on_meth_var(df, patients, cl_max_meth=0.2, cl_max_var=0.05, pl_min_meth=0,
                                        pl_min_var=0.2):
    """
    Label the df based on meth and var values, this is only good for single cell where the cl has low
    methylation and low variance and the pl has high var and high meth
    Due to the nature of the sc dataset each patient will provide different number of sequences and maybe
    different label so we try to only take sequences which either match in the label or missing the label,
    this again is only for a single cell data
    :param patients: Patients to take from the df
    :param df: The dataframe to label
    :param cl_max_meth: THe max meth for cl
    :param cl_max_var: The max variance for cl
    :param pl_min_meth: The min meth for pl
    :param pl_min_var: The minimum var for pl
    :return: The new dataframe with labels, removing CpG which didn't include
    """
    # Notice: this code is tailored made for our dataset
    if len(patients) > 2:
        sys.exit("this code only work for 2 training patient")

    for patient in patients:
        patient_num = patient[-2:]
        label = "label%s" % patient_num
        var_label = "var%s" % patient_num
        meth_label = "meth%s" % patient_num
        df[label] = 3  # misc
        df.loc[df[meth_label].isnull(), label] = 2  # empty

        df.loc[np.logical_and(df[meth_label] <= cl_max_meth, df[var_label] <= cl_max_var), label] = \
            consts.LABEL_HYPO_PRONE

        df.loc[np.logical_and(df[meth_label] >= pl_min_meth, df[var_label] >= pl_min_var), label] = \
            consts.LABEL_HYPO_RESISTANCE

    label1 = "label%s" % patients[0][-2:]
    label2 = "label%s" % patients[1][-2:]

    # Get all indexes of different combinations of labels
    partial_partial = np.logical_and(df[label1] == 0, df[label2] == 0)
    total_total = np.logical_and(df[label1] == 1, df[label2] == 1)
    partial_empty = np.logical_and(df[label1] == 0, df[label2] == 2)
    total_empty = np.logical_and(df[label1] == 1, df[label2] == 2)
    empty_total = np.logical_and(df[label1] == 2, df[label2] == 0)
    empty_partial = np.logical_and(df[label1] == 2, df[label2] == 1)
    good_index = partial_partial | total_total | partial_empty | total_empty | empty_total | empty_partial

    matching_labels_cpg = df[good_index]
    dropped_index = df.shape[0] - matching_labels_cpg.shape[0]

    df_list = []

    # Combine data from different patients
    for patient in patients:
        patient_num = patient[-2:]
        temp_df = pd.DataFrame()
        label = "label%s" % patient_num
        var_label = "var%s" % patient_num
        meth_label = "meth%s" % patient_num

        temp_df["meth"] = matching_labels_cpg[meth_label]
        temp_df["var"] = matching_labels_cpg[var_label]
        temp_df["sequence"] = matching_labels_cpg["sequence"]
        temp_df["ccpg"] = matching_labels_cpg["ccpg"]
        temp_df["label"] = matching_labels_cpg[label]
        temp_df["chr"] = matching_labels_cpg["chromosome"]
        temp_df["location"] = matching_labels_cpg["location"]

        df_list.append(temp_df)

    combined_df = pd.concat(df_list)
    good_index = np.logical_and(combined_df["label"] != 2, combined_df["label"] != 3)
    final_df = combined_df[good_index]
    dropped_index += combined_df.shape[0] - final_df.shape[0]

    print("During filtering we dropped %s CpG due to mismatch or invalid labels" % dropped_index)
    return final_df


def label_bulk_based_on_meth_var_flat(df, cl_max_meth=0.55, cl_min_var=0.11, pl_min_meth=0.75,
                                      pl_max_var=0.055):
    df.loc[np.logical_and(df["meth"] >= pl_min_meth,
                          df["var"] <= pl_max_var), "label"] = consts.LABEL_HYPO_RESISTANCE
    df.loc[np.logical_and(df["meth"] <= cl_max_meth,
                          df["var"] >= cl_min_var), "label"] = consts.LABEL_HYPO_PRONE

    filtered_df = df[~pd.isnull(df["label"])]
    return filtered_df


def label_bulk_based_on_meth_covariance_flat(df):
    # this is hard coded to the line equation we choose
    hypo_prone_index = df.apply(lambda x: x['coveriance'] > 0.05 * x['meth'] - 0.0113, axis=1)
    hypo_resistance_index = df.apply(lambda x: x['coveriance'] < 0.05 * x['meth'] - 0.0263, axis=1)
    df.loc[hypo_prone_index, "label"] = consts.LABEL_HYPO_PRONE
    df.loc[hypo_resistance_index, "label"] = consts.LABEL_HYPO_RESISTANCE

    # df.loc[np.logical_and(df["meth"] >= resistance_min_meth,
    #                       df["coveriance"] <= resistance_max_cov), "label"] = consts.LABEL_HYPO_RESISTANCE
    # df.loc[np.logical_and(df["meth"] <= prone_max_meth,
    #                       df["coveriance"] >= prone_min_cov), "label"] = consts.LABEL_HYPO_PRONE

    filtered_df = df[~pd.isnull(df["label"])]
    return filtered_df


def create_nn_dataset():
    args = parse_input()
    parse_format = args.parse_format

    # Read and add features
    df = pd.read_pickle(args.input_file)

    # We start with solo which are methylated in NC
    df["ccpg"] = df["sequence"].str.count("CG")
    df = df[df["ccpg"] < 4]

    if parse_format == consts.SCWGBS:
        df = df[df["orig_meth_avg"] >= 0.7]
    else:
        df = df[df["orig_meth"] >= 0.7]

    # Split the data to train and test based on pmd
    train, test = split_df_by_pmd(df)

    if parse_format == consts.SCWGBS:
        # Flat based on patients and add labels
        train = label_and_flat_sc_based_on_meth_var(train, TRAIN_PATIENT)
        test = label_and_flat_sc_based_on_meth_var(test, TEST_PATIENT)
    else:
        # Label based on values
        # train = label_bulk_based_on_meth_var_flat(train)
        # test = label_bulk_based_on_meth_var_flat(test)
        train = label_bulk_based_on_meth_covariance_flat(train)
        test = label_bulk_based_on_meth_covariance_flat(test)

    # Add reverse strand
    train = double_with_reverse_strand(train)
    test = double_with_reverse_strand(test)

    nn_dataset = {"train": train, "test": test}
    files_tools.save_pickle(file_path=os.path.join(args.output_folder, args.output_name), data=nn_dataset)

    print_statistics(train, test)


def split_bulk_dataset(df, test_size=0.2):
    lower_limit, upper_limit = df["pmd_index"].min() + 1, df["pmd_index"].max() - 1
    total_cpg = df.shape[0]
    current_test_cpg = np.sum(np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    while current_test_cpg / total_cpg < test_size:
        lower_limit += 1
        upper_limit -= 5
        current_test_cpg = np.sum(
            np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    while current_test_cpg / total_cpg > test_size:
        lower_limit -= 1
        current_test_cpg = np.sum(
            np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    lower_limit += 1
    while current_test_cpg / total_cpg > test_size:
        upper_limit += 1
        current_test_cpg = np.sum(
            np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit))

    upper_limit -= 1

    print("Going to use all pmd index: %s-%s" % (lower_limit, upper_limit))

    train_part = df[np.logical_and(df["pmd_index"] > lower_limit, df["pmd_index"] < upper_limit)]
    test_part = df[np.logical_or(df["pmd_index"] <= lower_limit, df["pmd_index"] >= upper_limit)]

    return train_part, test_part


def create_bulk_dataset():
    args = parse_input()
    parse_format = args.parse_format

    # Read and add features
    df = pd.read_pickle(args.input_file)

    # We start with solo which are methylated in NC
    df["ccpg"] = df["sequence"].str.count("CG")
    df = df[np.logical_and(df["ccpg"] < 4, df["orig_meth"] >= 0.7)]
    df = label_bulk_based_on_meth_covariance_flat(df)

    # Split the data to train and test based on pmd
    train, test = split_bulk_dataset(df)

    # Add reverse strand
    train = double_with_reverse_strand(train)
    test = double_with_reverse_strand(test)

    nn_dataset = {"train": train, "test": test}
    files_tools.save_pickle(file_path=os.path.join(args.output_folder, args.output_name), data=nn_dataset)

    print_statistics(train, test)


if __name__ == '__main__':
    create_bulk_dataset()
