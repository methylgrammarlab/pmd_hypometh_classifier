# TODO: how to write it ? some of the code is from deepripe but some is ours
"""

"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from classifier.plotseqlogo import seqlogo_fig
from commons import files_tools

sns.set()
sns.set_style('whitegrid')


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--interpretation_file', help='path for the input file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))

    args = parser.parse_args()
    return args


def plot_one_seq(seq, output, title, yl=None):
    fig = seqlogo_fig(seq[:, :], vocab="DNA", yl=yl, figsize=(20, 4), ncol=1, plot_name=title)
    fig.savefig(output)
    plt.close()


def plot_multi_seq(sequences_dict, number_of_seq, output_folder):
    """
    Plot the multiple sequences in one figure
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param number_of_seq: number of sequences in one figure
    :param output_folder: Output folder
    """
    for k in sequences_dict:
        ex_seq = sequences_dict[k][:number_of_seq]
        fig = seqlogo_fig(np.transpose(ex_seq[:, :, :], axes=(1, 2, 0)), vocab="DNA",
                          figsize=(8, ex_seq.shape[0]), ncol=1, yl=0.1,
                          plot_name="seq for top %s of type %s" % (number_of_seq, k))

        fig.savefig(os.path.join(output_folder, "seq_for_top_%s_of_type_%s" % (number_of_seq, k)))

    plt.close()


def plot_avg_sequence(sequences_dict, output_folder):
    """
    Plot the average sequence across 30 letters and all the sequence
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param output_folder: Output folder
    """
    for k in sequences_dict:
        ex_seq = sequences_dict[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, 60:90, :], axis=0).reshape(1, 30, 4), axes=(1, 2, 0))

        name = "Completely lost" if k == "cl" else "Partially lost"
        fig = seqlogo_fig(mean_seq, vocab="DNA", figsize=(20, 4), ncol=1,
                          plot_name="Average attribution score for prediction %s" % name)

        ax = fig.axes[0]
        ax.set_title("Average sequence for prediction %s" % name, fontsize=16)

        fig.savefig(os.path.join(output_folder, "Avg_seq_for_%s30.png" % k))
        plt.close()

    for k in sequences_dict:
        ex_seq = sequences_dict[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))

        fig = seqlogo_fig(mean_seq, vocab="DNA", figsize=(20, 4), ncol=1,
                          plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output_folder, "Avg_seq_for_%s.png" % k))
        plt.close()


def plot_avg_sequence_sw(sequences_dict, output_folder):
    """
    plot the avg sequence using SW, flatten the AT to W and CG to S
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param output_folder: Output folder
    """
    for k in sequences_dict:
        ex_seq = sequences_dict[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, 60:90, :], axis=0).reshape(1, 30, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            new_seq[i][0] = mean_seq[i][0] + mean_seq[i][3]
            new_seq[i][1] = mean_seq[i][1] + mean_seq[i][2]

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(20, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output_folder, "Avg_seq_for_%s_sw30.png" % k))

    for k in sequences_dict:
        ex_seq = sequences_dict[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            new_seq[i][0] = mean_seq[i][0] + mean_seq[i][3]
            new_seq[i][1] = mean_seq[i][1] + mean_seq[i][2]

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(20, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output_folder, "Avg_seq_for_%s_sw.png" % k))

    plt.close()


def plot_avg_sequence_sw_flatten_values(sequences_dict, output_folder):
    """
    plot the avg sequence using SW, flatten the AT to W and CG to S and combining both options to get one
    number per sequence place
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param output_folder: Output folder
    """
    for k in sequences_dict:
        ex_seq = sequences_dict[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, 60:90, :], axis=0).reshape(1, 30, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            w = mean_seq[i][0] + mean_seq[i][3]
            s = mean_seq[i][1] + mean_seq[i][2]
            delta = s - w
            sw_index = 1 if delta > 0 else 0
            new_seq[i][sw_index] = abs(delta)

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(8, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output_folder, "Avg_seq_for_%s_sw30_flatten.png" % k))
        # fig.show()

    for k in sequences_dict:
        ex_seq = sequences_dict[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            w = mean_seq[i][0] + mean_seq[i][3]
            s = mean_seq[i][1] + mean_seq[i][2]
            delta = s - w
            sw_index = 1 if delta > 0 else 0
            new_seq[i][sw_index] = abs(delta)

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(20, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output_folder, "Avg_seq_for_%s_sw_flatten.png" % k))

    plt.close()


def plot_distance_weight_two_sides(sequences_dict, output_folder):
    """
    Plot the integrated gradient value of each feature based on distance from center, two ways graph(-74->74)
    We wanted to see if there are indexes and some periodicity
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param output_folder: Output folder
    """
    for k in sequences_dict:
        class_type = "Completely lost" if k == "cl" else "Partially lost"
        ex_seq = np.abs(sequences_dict[k])

        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        seq_weight = np.sum(mean_seq, axis=1)
        middle = int(seq_weight.shape[0] / 2) - 1
        seq_weight[middle] = None
        seq_weight[middle + 1] = None

        x = np.arange(-74, 1).astype(np.int)
        x = np.append(x, x[::-1] * -1)

        x_ticks = [i for i in range(-70, 80, 10)]
        plt.xticks(x_ticks)
        plt.plot(x, seq_weight, '.-')
        plt.legend()
        plt.grid(axis="y")
        plt.xlabel("Distance from CpG Site", fontsize=12)
        plt.ylabel("Attribute score", fontsize=12)
        plt.title("Attribute score base on distance from CpG site for %s" % class_type, fontsize=14)

        plt.savefig(
            os.path.join(output_folder, "distance_importance_of_flanking_letters_type_%s_two_way.png" % k))
        plt.close()


def plot_distance_weight_one_side(sequences_dict, output_folder):
    """
    Plot the integrated gradient value of each feature based on distance from center, one way graph (0->74)
    We wanted to see if there are indexes and some periodicity
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param output_folder: Output folder
   """
    for k in sequences_dict:
        class_type = "partial lost" if k == "pl" else "completely lost"
        ex_seq = np.abs(sequences_dict[k])

        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        seq_weight = np.sum(mean_seq, axis=1)
        std_seq = np.std(mean_seq, axis=1)
        middle = int(seq_weight.shape[0] / 2) - 1
        seq_to_values = np.flip(seq_weight[:middle])
        seq_from_values = seq_weight[middle + 2:]
        seq_to_std = np.flip(std_seq[:middle])
        seq_from_std = std_seq[middle + 2:]
        x = np.arange(1, seq_from_values.shape[0] + 1)

        plt.errorbar(x, seq_to_values, seq_to_std, marker='^', label="to", alpha=0.5)
        plt.errorbar(x, seq_from_values, seq_from_std, marker='^', label="from", alpha=0.5)

        plt.legend()
        x_ticks = [i for i in range(1, 5)] + [i for i in range(5, 75, 5)]
        plt.xticks(x_ticks)
        plt.xlabel("Distance from CG")
        plt.ylabel("Importance shannon values")
        plt.title("Importance of flanking letters - %s" % (class_type))

        plt.savefig(os.path.join(output_folder,
                                 "distance_importance_of_flanking_letters_type_%s_one_way.png" % k))
        plt.close()


def print_each_seq(sequences_dict, output_folder):
    """
    Plot all the sequences on after the other
    :param sequences_dict: A dictionary with pl or cl as key and the integrated values results for each
    sequence in this label
    :param output_folder: Output folder
    """
    cl_list = []
    pl_list = []

    # Remove duplicates
    seq = None

    for i in range(sequences_dict["cl"].shape[0]):
        new_seq = sequences_dict["cl"][i]
        if np.all(new_seq == seq):
            continue
        else:
            cl_list.append(new_seq)
            seq = new_seq

    seq = None
    for i in range(sequences_dict["pl"].shape[0]):
        new_seq = sequences_dict["pl"][i]
        if np.all(new_seq == seq):
            continue
        else:
            pl_list.append(new_seq)
            seq = new_seq

    for i in range(1000):
        plot_one_seq(seq=cl_list[i], output=os.path.join(output_folder, "cl_seq_%s.png" % i),
                     title="CL seq num %s" % i, yl=0.1)

    for i in range(1000):
        plot_one_seq(seq=pl_list[i], output=os.path.join(output_folder, "pl_seq_%s.png" % i),
                     title="PL seq num %s" % i, yl=0.1)


def main():
    args = parse_input()

    ex_seq_d = files_tools.load_pickle(args.interpretation_file)
    new_d = {"cl": ex_seq_d["cl"], "pl": ex_seq_d["pl"]}

    plot_distance_weight_one_side(new_d, args.output_folder)
    plot_distance_weight_two_sides(new_d, args.output_folder)
    plot_multi_seq(new_d, 1000, args.output_folder)
    plot_avg_sequence(new_d, args.output_folder)
    plot_avg_sequence_sw(new_d, args.output_folder)
    plot_avg_sequence_sw_flatten_values(new_d, args.output_folder)
    print_each_seq(new_d, args.output_folder)


if __name__ == '__main__':
    main()
