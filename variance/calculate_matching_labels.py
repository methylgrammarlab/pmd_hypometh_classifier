import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

OUTPUT_FILE = "classifier_data_ccpg1.pkl"
TRANSLATION_TABLE = {84: 65, 65: 84, 67: 71, 71: 67}

TRAIN_PATIENT = ["CRC11", "CRC13"]
TEST_PATIENT = ["CRC01"]

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


def add_label(df, patient):
    # Group 1 (partial loss)- variance > 0.2
    # Group 2 (completely loss)- variance < 0.1 with meth < 0.2

    label = "label%s" % patient
    df[label] = 2  # other
    df.loc[df["var%s" % patient] >= 0.2, label] = 0  # partial
    df.loc[np.logical_and(df["meth%s" % patient] <= 0.2, df["var%s" % patient] <= 0.1), label] = 1  # total

    return df


def create_three_on_three(df, p1, p2):
    label1 = "label%s" % p1
    label2 = "label%s" % p2
    df = df[[label1, label2]]

    partial_partial = np.sum(np.logical_and(df[label1] == 0, df[label2] == 0))
    total_total = np.sum(np.logical_and(df[label1] == 1, df[label2] == 1))
    other_other = np.sum(np.logical_and(df[label1] == 2, df[label2] == 2))
    partial_total = np.sum(np.logical_and(df[label1] == 0, df[label2] == 1))
    partial_other = np.sum(np.logical_and(df[label1] == 0, df[label2] == 2))
    total_partial = np.sum(np.logical_and(df[label1] == 1, df[label2] == 0))
    total_other = np.sum(np.logical_and(df[label1] == 1, df[label2] == 2))
    other_total = np.sum(np.logical_and(df[label1] == 2, df[label2] == 1))
    other_partial = np.sum(np.logical_and(df[label1] == 2, df[label2] == 0))

    fdf = pd.DataFrame(index=["crc13_partial", "crc13_total", "crc13_other"],
                       columns=["crc11_partial", "crc11_total", "crc11_other"])

    fdf.loc["crc13_partial", "crc11_partial"] = partial_partial
    fdf.loc["crc13_total", "crc11_total"] = total_total
    fdf.loc["crc13_other", "crc11_other"] = other_other
    fdf.loc["crc13_partial", "crc11_total"] = partial_total
    fdf.loc["crc13_partial", "crc11_other"] = partial_other
    fdf.loc["crc13_total", "crc11_partial"] = total_partial
    fdf.loc["crc13_total", "crc11_other"] = total_other
    fdf.loc["crc13_other", "crc11_total"] = other_total
    fdf.loc["crc13_other", "crc11_partial"] = other_partial

    with open("matching_labels.csv", "w") as output:
        output.write(fdf.to_csv().replace("\r\n", "\n"))


def main():
    args = parse_input()

    # Read and add features
    df = pd.read_pickle(args.input_file)
    df["ccpg"] = df["sequence"].str.count("CG")

    # We start with solo which are methylated in NC
    df = df[df["ccpg"] == 1]
    df = df[df["nc_avg"] > 0.5]

    df = add_label(df, "11")
    df = add_label(df, "01")

    # Filter places we don't have info
    var_columns = df[["var11", "var01"]]
    not_null = np.mean(~pd.isnull(var_columns), axis=1) == 1
    print("We have %s(%s%%) shared cpg" % (np.sum(not_null), np.sum(not_null) * 100 / df.shape[0]))

    df = df[not_null]

    create_three_on_three(df, "11", "01")


if __name__ == '__main__':
    main()
