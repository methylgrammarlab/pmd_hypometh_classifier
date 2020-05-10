import argparse
import os
import pickle
import random
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from classifier.utils import precision

plt.style.use('seaborn-deep')


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', dest='data_path', type=str, help='input path')

    return parser.parse_args()


def get_data(path_to_data):
    with open(path_to_data, "rb") as path_to_data_h:
        data = pickle.load(path_to_data_h)

    train_data, test_data = data["train"], data["test"]

    return train_data, test_data


def accuracy(y_true, y_pred):
    diff = y_true - y_pred
    tptn = np.sum(diff == 0)
    return tptn / y_true.shape[0]


def test_majority_vote(train_data, test_data, output):
    y_train = train_data["label"]
    y_test = test_data["label"]

    mean = y_train.mean()
    r_mean = np.round(mean)
    y_pred = np.full_like(y_test, r_mean)

    output.write("majority vote, %s, %s, %s\n"
                 % (accuracy(y_test, y_pred), recall(y_test, y_pred), precision(y_test, y_pred)))


def test_based_on_4(train_data, test_data, output):
    y_test = test_data["label"]

    mean_by_seq = train_data.groupby("seq4").mean()
    with open("mean_by_seq4_train.csv", "w") as f:
        f.write(mean_by_seq.to_csv().replace("\r\n", "\n"))

    seq_dict = np.round(mean_by_seq)["label"].to_dict()
    y_pred = [seq_dict[seq] for seq in test_data["seq4"]]
    y_pred = np.array(y_pred)

    output.write("memory based on 1 flank, %s, %s, %s\n"
                 % (accuracy(y_test, y_pred), recall(y_test, y_pred), precision(y_test, y_pred)))

    mean_by_seq = test_data.groupby("seq4").mean()
    with open("mean_by_seq4_test.csv", "w") as f:
        f.write(mean_by_seq.to_csv().replace("\r\n", "\n"))


def test_based_on_10(train_data, test_data, output):
    y_test = test_data["label"]

    mean_by_seq = train_data.groupby("seq10").mean()

    seq_dict = np.round(mean_by_seq)["label"].to_dict()
    y_pred = []
    r = 0
    for seq in test_data["seq10"]:
        if seq in seq_dict:
            v = seq_dict[seq]
        else:
            v = random.randint(0, 1)
            r += 1

        y_pred.append(v)

    y_pred = np.array(y_pred)

    output.write("memory based on 4 flank, %s, %s, %s\n"
                 % (accuracy(y_test, y_pred), recall(y_test, y_pred), precision(y_test, y_pred)))


def test_based_on_6(train_data, test_data, output):
    y_test = test_data["label"]

    train_data["seq6"] = train_data["seq10"].str[2:-2]
    test_data["seq6"] = test_data["seq10"].str[2:-2]

    mean_by_seq = train_data.groupby("seq6").mean()

    seq_dict = np.round(mean_by_seq)["label"].to_dict()
    y_pred = []
    r = 0
    for seq in test_data["seq6"]:
        if seq in seq_dict:
            v = seq_dict[seq]
        else:
            v = random.randint(0, 1)
            r += 1

        y_pred.append(v)

    y_pred = np.array(y_pred)

    output.write("memory based on 2 flank, %s, %s, %s\n"
                 % (accuracy(y_test, y_pred), recall(y_test, y_pred), precision(y_test, y_pred)))


def test_based_on_8(train_data, test_data, output):
    y_test = test_data["label"]

    train_data["seq8"] = train_data["seq10"].str[1:-1]
    test_data["seq8"] = test_data["seq10"].str[1:-1]

    mean_by_seq = train_data.groupby("seq8").mean()

    seq_dict = np.round(mean_by_seq)["label"].to_dict()
    y_pred = []
    r = 0
    for seq in test_data["seq8"]:
        if seq in seq_dict:
            v = seq_dict[seq]
        else:
            v = random.randint(0, 1)
            r += 1

        y_pred.append(v)

    y_pred = np.array(y_pred)

    output.write("memory based on 3 flank, %s, %s, %s\n"
                 % (accuracy(y_test, y_pred), recall(y_test, y_pred), precision(y_test, y_pred)))


def test_agreement_by_10(train_data):
    mean_by_seq = train_data.groupby("seq10").mean()
    values_rounded_to_2 = np.round(mean_by_seq["label"] * 100) / 100

    _ = plt.hist(values_rounded_to_2, bins='auto')
    plt.style.use('ggplot')
    plt.title("Distribution of labels on sequence of 10 ")
    plt.savefig("dist_based_on_10.png")
    plt.close()


def test_agreement_by_150(train_data):
    mean_by_seq = train_data.groupby("sequence").mean()
    values_rounded_to_2 = np.round(mean_by_seq["label"] * 100) / 100

    _ = plt.hist(values_rounded_to_2, bins='auto')
    plt.style.use('ggplot')
    plt.title("Distribution of labels on sequence of 150")
    plt.savefig("dist_based_on_150.png")
    plt.close()


def test_based_on_10_remove_uncertain(train_data, test_data):
    y_test = test_data["label"]

    mean_by_seq = train_data.groupby("seq10").mean()

    mean_av = mean_by_seq["label"]
    mean_av = mean_av[np.logical_or(mean_av < 0.4, mean_av > 0.6)]
    mean_av = np.round(mean_av)
    seq_dict = mean_av.to_dict()

    y_pred = []
    y_test_l = []
    r = 0
    for i in range(test_data.shape[0]):
        seq = test_data.iloc[i]["seq10"]
        l = test_data.iloc[i]["label"]
        if seq in seq_dict:
            v = seq_dict[seq]
        else:
            continue

        y_pred.append(v)
        y_test_l.append(l)

    y_pred = np.array(y_pred)
    y_test_l = np.array(y_test_l)

    print("Using memory based on 4 flank: Accuracy: %s. Recall: %s. Precision: %s. Random: %s(%s%%)" %
          (accuracy(y_test_l, y_pred), recall(y_test_l, y_pred), precision(y_test_l, y_pred),
           r, r / y_pred.shape[0] * 100))


def test_by_given_length(train_data, test_data, output, i, dict_4):
    y_test = test_data["label"]
    label = "seq_%s" % i

    missing_values = 0
    # majority = np.round(np.mean(train_data["label"]))

    if i == 0:
        train_data[label] = train_data["sequence"]
        test_data[label] = test_data["sequence"]

    else:
        train_data[label] = train_data["sequence"].str[i:-i]
        test_data[label] = test_data["sequence"].str[i:-i]

    mean_by_seq = train_data.groupby(label).mean()
    seq_dict = np.round(mean_by_seq)["label"].to_dict()

    y_pred = []
    y_test_l = []
    for j in range(test_data.shape[0]):
        seq = test_data.iloc[j][label]
        truth = y_test.iloc[j]

        if seq in seq_dict:
            y_pred.append(seq_dict[seq])
            y_test_l.append(truth)
        else:
            m = int(len(seq) / 2)
            min_seq = seq[m - 4:m + 4]

            y_pred.append(dict_4[min_seq])
            y_test_l.append(truth)
            missing_values += 1

    y_pred = np.array(y_pred)
    y_test_l = np.array(y_test_l)
    perc_of_data_exists = 100 - missing_values / y_test.shape[0] * 100

    output.write("%s,%s,%s\n" % (150 - 2 * i, accuracy(y_test_l, y_pred), perc_of_data_exists))


def main():
    args = format_args()
    train_data, test_data = get_data(args.data_path)

    # test_based_on_10_remove_uncertain(train_data, test_data)
    # test_agreement_by_10(train_data)
    # test_agreement_by_150(train_data)
    #
    # with open("output.csv", "w") as output:
    #     output.write("test type, accuracy, recall, percision\n")
    #     test_majority_vote(train_data, train_data, output)
    #     test_based_on_4(train_data, train_data, output)
    #     test_based_on_6(train_data, train_data, output)
    #     test_based_on_8(train_data, train_data, output)
    #     test_based_on_10(train_data, train_data, output)
    #
    #     test_majority_vote(train_data, test_data, output)
    #     test_based_on_4(train_data, test_data, output)
    #     test_based_on_6(train_data, test_data, output)
    #     test_based_on_8(train_data, test_data, output)
    #     test_based_on_10(train_data, test_data, output)

    with open("accuracy_by_length60-75.csv", "w") as output:
        output.write("length_of_seq, accuracy, percentage of seq\n")
        label = "seq_d%s" % 4
        train_data[label] = train_data["sequence"].str[71:-71]

        mean_by_seq = train_data.groupby(label).mean()
        seq_dict = np.round(mean_by_seq)["label"].to_dict()

        for i in range(60, 75):
            print(i)
            test_by_given_length(train_data, test_data, output, i, seq_dict)


if __name__ == '__main__':
    main()
