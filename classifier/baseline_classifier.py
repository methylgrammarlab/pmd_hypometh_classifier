"""
This is a baseline classifier to try and predict the labels based on a simple flanking
This was used to get simple idea of the data compare to a NN
"""

import argparse
import os
import random
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from classifier import utils, create_data

SS_RE = re.compile("[CG]CG[CG]")
WW_RE = re.compile("[AT]CG[AT]")
SWCGWS_RE = re.compile("[CG][AT]CG[AT][CG]")
SSCGSS_RE = re.compile("[CG][CG]CG[CG][CG]")
WSCGSW_RE = re.compile("[AT][CG]CG[CG][AT]")
WWCGWW_RE = re.compile("[AT][AT]CG[AT][AT]")
CWCGWG_RE = re.compile("C[AT]CG[AT]G")
SW_RE = re.compile("[CG]CG[AT]")
WS_RE = re.compile("[AT]CG[CG]")

plt.style.use('seaborn-deep')


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', dest='data_path', type=str, help='input path')

    return parser.parse_args()


def print_results(test_name, y_test, y_pred):
    """
    Print to the screen the results of the a specific test
    :param test_name: The test name
    :param y_test: y_test labels
    :param y_pred: y_perd labels
    :return: Nothing
    """
    print("{title},{accuracy},{recall},{precision}".format(**{"title": test_name,
                                                              "accuracy": utils.accuracy(y_test, y_pred),
                                                              "recall": utils.recall_TP(y_test, y_pred),
                                                              "precision": utils.precision(y_test, y_pred)}))


def test_majority_vote(train_data, test_data):
    """
    Test the prediction using majority vote
    :param train_data: The train data, the majority vote will be calculated according to it
    :param test_data: The test data
    """
    y_train, y_test = train_data["label"], test_data["label"]

    mean = y_train.mean()
    r_mean = np.round(mean)
    y_pred = np.full_like(y_test, r_mean)
    print_results(test_name="majority vote", y_test=y_test, y_pred=y_pred)


def test_based_on_flank_num(train_data, test_data, flank_num):
    """
    Test the prediction using specific flank number and dictionary, this might take a lot of time or be
    very random if the flank number is large and we don't have many sequences matching in the train-test
    :param train_data: The train data, the majority vote will be calculated according to it
    :param test_data: The test data
    :param flank_num: Length of flank from CpG site
    """
    y_test = test_data["label"]
    seq_label = "temp%s" % flank_num
    center = int(len(train_data.iloc[0]["sequence"]) / 2)

    train_data[seq_label] = train_data["sequence"].str[center - flank_num - 1:center + 1 + flank_num]
    test_data[seq_label] = test_data["sequence"].str[center - flank_num - 1:center + 1 + flank_num]

    mean_by_seq = train_data.groupby(seq_label).mean()

    seq_dict = np.round(mean_by_seq)["label"].to_dict()
    y_pred = np.array([seq_dict.get(seq, random.randint(0, 1)) for seq in test_data[seq_label]])

    print_results(test_name="%s flank" % flank_num, y_test=y_test, y_pred=y_pred)


def test_based_on_wcgw_and_scgs(df):
    """
    Test the prediction using specific the WCGW and SCGW rules we learned, trying to see if this is enough
    or if the NN can\will find something more sophisticated
    This is based on the entire df and not only on the test, train
    :param df: the data frame with all the data
    """
    y_test = df["label"]
    x_test = list(df["sequence"])
    y_pred = []

    # Trying many regex
    for i in x_test:
        if WW_RE.search(i) is not None:
            if CWCGWG_RE.search(i) is not None:
                y_pred.append(create_data.LABEL_COMPLETELY_LOST)

            elif WWCGWW_RE.search(i) is not None:
                y_pred.append(create_data.LABEL_PARTIAL_LOST)
            else:
                y_pred.append(create_data.LABEL_COMPLETELY_LOST)

        elif SS_RE.search(i) is not None:
            if SSCGSS_RE.search(i) is not None:
                y_pred.append(create_data.LABEL_COMPLETELY_LOST)
            else:
                y_pred.append(create_data.LABEL_PARTIAL_LOST)

        elif SW_RE.search(i) is not None or WS_RE.search(i) is not None:
            y_pred.append(create_data.LABEL_PARTIAL_LOST)

        else:
            y_pred.append(create_data.LABEL_PARTIAL_LOST)

    y_pred = np.array(y_pred)
    print_results(test_name="Using WCGW\SCGW", y_test=y_test, y_pred=y_pred)


def main():
    args = format_args()
    train_data, test_data = utils.get_train_test_data(args.data_path)

    test_majority_vote(train_data, train_data)
    test_based_on_flank_num(train_data, train_data, flank_num=1)
    test_based_on_flank_num(train_data, train_data, flank_num=2)
    test_based_on_flank_num(train_data, train_data, flank_num=3)
    test_based_on_flank_num(train_data, train_data, flank_num=4)

    test_majority_vote(train_data, test_data)
    test_based_on_flank_num(train_data, test_data, flank_num=1)
    test_based_on_flank_num(train_data, test_data, flank_num=2)
    test_based_on_flank_num(train_data, test_data, flank_num=3)
    test_based_on_flank_num(train_data, test_data, flank_num=4)

    test_based_on_wcgw_and_scgs(test_data)


if __name__ == '__main__':
    main()
