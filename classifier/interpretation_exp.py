import argparse
import os
import pickle
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from classifier.plotseqlogo import seqlogo_fig
from classifier.utils import vecs2motif
from commons.data_tools import counter_to_csv

sns.set()
sns.set_style('whitegrid')


NCOL = 2


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--interpretation_file', help='path for the input file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))

    args = parser.parse_args()
    return args


def plot_seq_alone(seq, id, output, k, yl=None):
    fig = seqlogo_fig(seq[:, :], vocab="DNA", yl=yl,
                      figsize=(20, 4), ncol=1,
                      plot_name="Seq for top %s of type %s" % (id, k))

    fig.savefig(os.path.join(output, "Seq_for_top_%s_of_type_%s" % (id, k)))
    plt.close()


def plot_sequences(ex_seq_d, number_of_seq, output):
    for k in ex_seq_d:
        ex_seq = ex_seq_d[k][:number_of_seq]
        fig = seqlogo_fig(np.transpose(ex_seq[:, :, :], axes=(1, 2, 0)), vocab="DNA",
                          figsize=(8, ex_seq.shape[0]), ncol=1, yl=0.1,
                          plot_name="Seq for top %s of type %s" % (number_of_seq, k))

        fig.savefig(os.path.join(output, "Seq_for_top_%s_of_type_%s" % (number_of_seq, k)))
        # fig.show()
    plt.close()


def plot_avg_sequence(ex_seq_d, output):
    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, 60:90, :], axis=0).reshape(1, 30, 4), axes=(1, 2, 0))

        name = "Completely lost" if k == "cl" else "Partially lost"
        fig = seqlogo_fig(mean_seq, vocab="DNA", figsize=(20, 4), ncol=1,
                          plot_name="Average attribution score for prediction %s" % name)

        ax = fig.axes[0]
        # x = np.arange(-1 * int(mean_seq.shape[0] / 2), 1).astype(np.int)
        # x = np.append(x,x[::-1] * -1)
        # ax.set_xticks(x)
        # ax.set_xlabel("Sequence")
        # ax.set_ylabel("Attribute score (shannon values)")
        ax.set_title("Average sequence for prediction %s" % name, fontsize=16)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s30.png" % k))
        plt.close()

    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))

        fig = seqlogo_fig(mean_seq, vocab="DNA", figsize=(20, 4), ncol=1,
                          plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s.png" % k))
        plt.close()


def plot_avg_sequence_sw(ex_seq_d, output):
    """
    plot the avg sequence using SW, flatten the AT to W and CG to S
    """
    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, 60:90, :], axis=0).reshape(1, 30, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            new_seq[i][0] = mean_seq[i][0] + mean_seq[i][3]
            new_seq[i][1] = mean_seq[i][1] + mean_seq[i][2]

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(20, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s_sw30.png" % k))
        # fig.show()

    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            new_seq[i][0] = mean_seq[i][0] + mean_seq[i][3]
            new_seq[i][1] = mean_seq[i][1] + mean_seq[i][2]

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(20, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s_sw.png" % k))
        # fig.show()
    plt.close()


def plot_avg_sequence_sw_flatten_values(ex_seq_d, output):
    """
    plot the avg sequence using SW, flatten the AT to W and CG to S
    """
    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, 60:90, :], axis=0).reshape(1, 30, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            w = mean_seq[i][0] + mean_seq[i][3]
            s = mean_seq[i][1] + mean_seq[i][2]
            delta = s - w
            sw_index = 1 if delta > 0 else 0
            new_seq[i][sw_index] = abs(delta)

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(8, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s_sw30_flatten.png" % k))
        # fig.show()

    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        new_seq = np.zeros_like(mean_seq)
        for i in range(mean_seq.shape[0]):
            w = mean_seq[i][0] + mean_seq[i][3]
            s = mean_seq[i][1] + mean_seq[i][2]
            delta = s - w
            sw_index = 1 if delta > 0 else 0
            new_seq[i][sw_index] = abs(delta)

        fig = seqlogo_fig(new_seq, vocab="DNAWS", figsize=(20, 4), ncol=1, plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s_sw_flatten.png" % k))
        # fig.show()
    plt.close()


def hist_3d(ex_seq_d, number_of_seq, output):
    for k in ex_seq_d:
        if 'failed' in k:
            continue
        bins = 10
        ex_seq = ex_seq_d[k][:number_of_seq]
        tr = np.transpose(ex_seq[:, 73:77, :], axes=(1, 2, 0))  # change
        for loc_ind, loc in enumerate(tr):
            loc_df = pd.DataFrame(loc, index=["A", "C", "G", "T"]).T
            # loc_df.plot.hist(alpha=0.9, histtype='step')  # plots just the outline
            # loc_df.plot.hist(alpha=0.3, histtype='stepfilled')  # plots just the inside
            # loc_df.plot.hist(alpha=0.3, histtype='stepfilled', ec="k")  # plots with black outline
            plt.hist([loc_df.loc[:, 'A'], loc_df.loc[:, 'C'], loc_df.loc[:, 'G'], loc_df.loc[:, 'T']],
                     label=["A", "C", "G", "T"], alpha=0.5, align="mid")
            plt.legend(loc='upper left')
            plt.yscale('log')
            ind = str(int(loc_ind - ((len(tr) / 2) - 1))) if loc_ind <= ((len(tr) / 2) - 1) else '+' + str(
                int(loc_ind - len(tr) / 2))
            plt.title("Seq for location %s of type %s with %s seqs" % (ind, k, number_of_seq))
            path = os.path.join('exp_output', k, "Seq_for_location_%s_of_type_%s_%s_seqs" % (ind, k, number_of_seq))
            plt.savefig(path)
            plt.close()
        pass


def plot_distance_weight_combine(ex_seq_d, output):
    for k in ex_seq_d:
        class_type = "Completely lost" if k == "cl" else "Partially lost"
        ex_seq = np.abs(ex_seq_d[k])

        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        seq_weight = np.sum(mean_seq, axis=1)
        std_seq = np.std(mean_seq, axis=1)
        middle = int(seq_weight.shape[0] / 2) - 1
        seq_weight[middle] = None
        seq_weight[middle + 1] = None

        # seq_to_values = np.flip(seq_weight[:middle])
        # seq_from_values = seq_weight[middle + 2:]
        # seq_to_std = np.flip(std_seq[:middle])
        # seq_from_std = std_seq[middle + 2:]
        x = np.arange(-74, 1).astype(np.int)
        x = np.append(x, x[::-1] * -1)

        # plt.errorbar(x, seq_to_values, seq_to_std, marker='^', label="to", alpha=0.5)
        # plt.errorbar(x, seq_from_values, seq_from_std, marker='^', label="from", alpha=0.5)
        # # plt.errorbar(x, np.mean([seq_to_values, seq_from_values],axis=0),
        # #              np.mean([seq_to_std, seq_from_std], axis=0),  marker='^', label="combine", alpha=0.3)
        x_ticks = [i for i in range(-70, 80, 10)]
        plt.xticks(x_ticks)
        plt.plot(x, seq_weight, '.-')
        plt.legend()
        plt.grid(axis="y")
        plt.xlabel("Distance from CpG Site", fontsize=12)
        plt.ylabel("Attribute score", fontsize=12)
        plt.title("Attribute score base on distance from CpG site for %s" % class_type, fontsize=14)

        plt.savefig(os.path.join(output, "pres_distance_importance_of_flanking_letters_type_%s.png" % k))
        plt.close()


def plot_distance_weight(ex_seq_d, output):
    """
    trying to plot the importance of each flanking vlue
    """
    for k in ex_seq_d:
        class_type = "partial lost" if k == "pl" else "completely lost"
        ex_seq = np.abs(ex_seq_d[k])

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
        # plt.errorbar(x, np.mean([seq_to_values, seq_from_values],axis=0),
        #              np.mean([seq_to_std, seq_from_std], axis=0),  marker='^', label="combine", alpha=0.3)
        plt.legend()
        x_ticks = [i for i in range(1, 5)] + [i for i in range(5, 75, 5)]
        plt.xticks(x_ticks)
        plt.xlabel("Distance from CG")
        plt.ylabel("Importance shannon values")
        plt.title("Importance of flanking letters - %s" % (class_type))

        plt.savefig(os.path.join(output, "distance_importance_of_flanking_letters_type_%s.png" % k))
        plt.close()


def print_each_seq(ex_seq_d, output_folder):
    cl_list = []
    pl_list = []

    # Remove duplicates
    seq = None

    for i in range(ex_seq_d["cl"].shape[0]):
        new_seq = ex_seq_d["cl"][i]
        if np.all(new_seq == seq):
            continue
        else:
            cl_list.append(new_seq)
            seq = new_seq

        #
    seq = None
    for i in range(ex_seq_d["pl"].shape[0]):
        new_seq = ex_seq_d["pl"][i]
        if np.all(new_seq == seq):
            continue
        else:
            pl_list.append(new_seq)
            seq = new_seq

    for i in range(1000):
        plot_seq_alone(cl_list[i], i, os.path.join(output_folder, "cl"), "cl", yl=0.1)

    for i in range(1000):
        plot_seq_alone(pl_list[i], i, os.path.join(output_folder, "pl"), "pl", yl=0.1)
    plt.close()


def find_motifs(new_d, output_folder):
    for k in new_d:
        seqs = vecs2motif(new_d[k])
        seqs_splitted = []
        for seq in seqs:
            seq_s = seq.split("_")
            seq_l = [i for i in seq_s if i != "" and len(i) > 1]
            seqs_splitted.append(list(set(seq_l)))

        motif_counter = Counter()
        for seq in seqs_splitted:
            motif_counter.update(seq)

        print(len(seqs_splitted))
        counter_to_csv(motif_counter, os.path.join(output_folder, k + ".csv"))


def plot_distance_weight_combine_single_letter(ex_seq_d, output):
    np.nanmean(np.where(b >= 0, b, np.nan), axis=1)
    for k in ex_seq_d:
        class_type = "Completely lost" if k == "cl" else "Partially lost"
        ex_seq = np.abs(ex_seq_d[k])

        mean_seq = np.transpose(np.mean(ex_seq[:, :, :], axis=0).reshape(1, 150, 4), axes=(1, 2, 0))
        seq_weight = np.sum(mean_seq, axis=1)
        std_seq = np.std(mean_seq, axis=1)
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
        plt.xlabel("Distance from CpG Site")
        plt.ylabel("Attribute score")
        plt.title("Attribute score base on distance from CpG site for %s" % class_type)

        plt.savefig(os.path.join(output, "pres_distance_importance_of_flanking_letters_type_%s.png" % k))
        plt.close()

def main():
    args = parse_input()

    with open(args.interpretation_file, "rb") as f:
        ex_seq_d = pickle.load(f)

    new_d = {"cl": ex_seq_d["cl"], "pl": ex_seq_d["pl"]}

    # hist_3d(ex_seq_d, 1000, args.output_folder)
    # plot_distance_weight(new_d, args.output_folder)
    # plot_distance_weight_combine(new_d, args.output_folder)
    # plot_distance_weight_combine_single_letter(new_d, args.output_folder)
    # plot_sequences(new_d, 1000, args.output_folder)
    plot_avg_sequence(new_d, args.output_folder)
    # plot_avg_sequence_sw(new_d, args.output_folder)
    # plot_avg_sequence_sw_flatten_values(new_d, args.output_folder)
    # print_each_seq(new_d, args.output_folder)
    # find_motifs(new_d, args.output_folder)


if __name__ == '__main__':
    main()
