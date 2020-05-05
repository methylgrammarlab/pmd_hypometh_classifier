import argparse
import os
import pickle
import sys

import numpy as np

from classifier.plotseqlogo import seqlogo_fig

NCOL = 2


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--interpretation_file', help='path for the input file', required=True)

    args = parser.parse_args()
    return args


def plot_sequences(ex_seq_d, number_of_seq, output):
    for k in ex_seq_d:
        ex_seq = ex_seq_d[k][:number_of_seq]
        fig = seqlogo_fig(np.transpose(ex_seq[:, 60:90, :4], axes=(1, 2, 0)), vocab="DNA",
                          figsize=(8, ex_seq.shape[0] / NCOL + 2), ncol=NCOL,
                          plot_name="Seq for top %s of type %s" % (number_of_seq, k))

        fig.savefig(os.path.join(output, "Seq_for_top_%s_of_type_%s" % (number_of_seq, k)))


def plot_avg_sequence(ex_seq_d, output):
    for k in ex_seq_d:
        ex_seq = ex_seq_d[k]
        fig = seqlogo_fig(np.transpose(np.mean(ex_seq[:, 60:90, :4], axis=0).reshape(1, 30, 4),
                                       axes=(1, 2, 0)), vocab="DNA", figsize=(8, 4), ncol=1,
                          plot_name="Avg seq for %s" % k)

        fig.savefig(os.path.join(output, "Avg_seq_for_%s.png" % k))


def main():
    args = parse_input()

    with open(args.interpretation_file, "rb") as f:
        ex_seq_d = pickle.load(f)

    plot_sequences(ex_seq_d, 10, args.output_folder)
    plot_avg_sequence(ex_seq_d, args.output_folder)


if __name__ == '__main__':
    main()
