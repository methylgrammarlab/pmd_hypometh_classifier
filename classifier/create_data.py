"""
Create dataset for the nn, adding labels, splitting to test and train in a way that will consider the pmd
and add reverse strand.
This code is only for the scWGBS
"""
# TODO: can we combine with the zhou?

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import data_tools, files_tools

# scWGBS 
TRAIN_PATIENT = ["CRC01", "CRC11"]
TEST_PATIENT = ["CRC10", "CRC13"]

LABEL_PARTIAL_LOST = 0  # Group 0 - partial loss
LABEL_COMPLETELY_LOST = 1  # Group 1 - comp loss


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input_file', help='path for the input file', required=True)
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


def split_df_by_pmd(df, size_wheel=1.8):
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
    reversed_seq = data_tools.get_reverse_seq(df["sequence"])
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
    df.loc[np.logical_and(df["meth"] >= pl_min_meth, df["var"] >= pl_min_var), "label"] = LABEL_PARTIAL_LOST
    df.loc[np.logical_and(df["meth"] <= cl_max_meth, df["var"] <= cl_max_var), "label"] = \
        LABEL_COMPLETELY_LOST

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
            LABEL_COMPLETELY_LOST

        df.loc[np.logical_and(df["meth"] >= pl_min_meth, df["var"] >= pl_min_var), label] = \
            LABEL_PARTIAL_LOST

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
    return combined_df[good_index]


def create_scwgbs_nn_dataset():
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)

    # We start with solo which are methylated in NC
    df["ccpg"] = df["sequence"].str.count("CG")
    df = df[df["ccpg"] == 1]
    df = df[df["nc_avg"] >= 0.6]

    # Split the data to train and test based on pmd
    train, test = split_df_by_pmd(df)

    # Flat based on patients and add labels
    train = label_and_flat_sc_based_on_meth_var(train, TRAIN_PATIENT)
    test = label_and_flat_sc_based_on_meth_var(test, TEST_PATIENT)

    # Add reverse strand
    train = double_with_reverse_strand(train)
    test = double_with_reverse_strand(test)

    nn_dataset = {"train": train, "test": test}
    files_tools.save_pickle(file_path=os.path.join(args.output_folder, args.output_name), data=nn_dataset)

    print_statistics(train, test)


def main():
    create_scwgbs_nn_dataset()


if __name__ == '__main__':
    main()
