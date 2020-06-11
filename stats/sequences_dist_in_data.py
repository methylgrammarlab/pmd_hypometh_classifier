import argparse
import os
import sys

import matplotlib as mpl

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TRANSLATION_TABLE = {84: 65, 65: 84, 67: 71, 71: 67}
LETTER_DICT = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4}
# CMAP = mpl.colors.ListedColormap(['blue', 'orange', 'red', 'green'])
CMAP = mpl.colors.ListedColormap(['blue', 'orange', 'orange', 'blue'])


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input_file', help='path for the input file', required=True)

    args = parser.parse_args()
    return args


def get_reverse(sequence):
    flipped = sequence.translate(TRANSLATION_TABLE)[::-1]
    return flipped


def get_info(df, sequence):
    seq_type = sequence.replace("A", "W")
    seq_type = seq_type.replace("T", "W")
    seq_type = seq_type.replace("C", "S")
    seq_type = seq_type.replace("G", "S")

    reversed = get_reverse(sequence)
    seq_df = df[df["sequence"].str.count(sequence) == 1]
    total = df.shape[0]
    num_of_sequences = seq_df.shape[0]
    p_seq = num_of_sequences / total
    num_of_pl_out_of_seq = np.sum(seq_df["label"] == 0)
    p_pl_out_of_seq = num_of_pl_out_of_seq / num_of_sequences
    num_of_cl_out_of_seq = np.sum(seq_df["label"] == 1)
    p_cl_out_of_seq = num_of_cl_out_of_seq / num_of_sequences

    p_pl_out_of_pl = np.sum(seq_df["label"] == 0) / np.sum(df["label"] == 0)
    p_cl_out_of_cl = np.sum(seq_df["label"] == 1) / np.sum(df["label"] == 1)

    #  format ,#sequnces, %sequences, #pl from format seq, %pl from format seq, #cl from format seq,
    # %cl from format seq, %format pl out of pl,%format cl out of cl
    return "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (seq_type, sequence, reversed, num_of_sequences, p_seq,
                                                   num_of_pl_out_of_seq,
                                                   num_of_cl_out_of_seq, p_pl_out_of_seq, p_cl_out_of_seq,
                                                   p_pl_out_of_pl, p_cl_out_of_cl)


def create_info():
    args = parse_input()
    letters = ["A", "C", "G", "T"]

    with open(args.input_file, "rb") as input_file:
        data = pickle.load(input_file)

    df = pd.concat([data["train"], data["test"]])
    df = df[["sequence", "label"]]

    with open(r"sequence_dist\2.csv", "w") as dist2:
        dist2_l = []
        sequences = set()

        for i in letters:
            for j in letters:
                seq = i + "CG" + j
                reversed = get_reverse(seq)
                if reversed in sequences:
                    continue

                sequences.add(reversed)
                sequences.add(seq)
                dist2_l.append(get_info(df, seq))

        dist2.write("seq_type, format ,reversed ,#sequnces, %sequences, #pl from format seq, "
                    "#cl from format seq, "
                    "%pl from format seq, "
                    "%cl from format seq, %format pl out of pl,%format cl out of cl\n")
        dist2.writelines(dist2_l)

    print("start3")
    with open(r"sequence_dist\3.csv", "w") as dist2:
        dist2_l = []
        sequences = set()
        for i in letters:
            for j in letters:
                for i1 in letters:
                    for j1 in letters:
                        seq = i + i1 + "CG" + j + j1

                        if seq.count("CG") > 1:
                            continue

                        reversed = get_reverse(seq)
                        if reversed in sequences:
                            continue

                        sequences.add(reversed)
                        sequences.add(seq)

                        dist2_l.append(get_info(df, seq))

        dist2.write(
            "seq_type, format ,reversed ,#sequnces, %sequences, #pl from format seq, #cl from format seq, "
            "%pl from format seq, "
            "%cl from format seq, %format pl out of pl,%format cl out of cl\n")
        dist2.writelines(dist2_l)

    print("start4")
    with open(r"sequence_dist\4.csv", "w") as dist2:
        dist2_l = []
        sequences = set()
        for i in letters:
            for j in letters:
                for i1 in letters:
                    for j1 in letters:
                        for i2 in letters:
                            for j2 in letters:
                                seq = i + i1 + i2 + "CG" + j + j1 + j2

                                if seq.count("CG") > 1:
                                    continue

                                reversed = get_reverse(seq)
                                if reversed in sequences:
                                    continue

                                sequences.add(reversed)
                                sequences.add(seq)

                                dist2_l.append(get_info(df, seq))

        dist2.write(
            "seq_type, format ,reversed  ,#sequnces, %sequences, #pl from format seq, #cl from format seq, "
            "%pl from format seq, "
            "%cl from format seq, %format pl out of pl,%format cl out of cl\n")
        dist2.writelines(dist2_l)


def split_word(words):
    output = []
    for word in words:
        i = []
        try:
            for s in word:
                i.append(LETTER_DICT[s])
        except Exception:
            pass
        output.append(i)
    return output


def get_ticks(word):
    length = len(word)
    start = -1 * length / 2 + 1
    a = np.arange(start, 1).astype(np.int)
    return np.append(a, a[::-1] * -1)


def plot_info():
    args = parse_input()
    data = pd.read_csv(r"H:\Study\university\Computational-Biology\Year "
                       r"3\Projects\proj_scwgbs\stats\sequence_dist\flank3.csv")
    data.loc[data["ratio"] == "#NUM!", "ratio"] = 1000
    data["ratio"] = data["ratio"].astype(np.float)
    data = data.sort_values('ratio')
    words = split_word(data[' format '])
    v = data["ratio"] * 100
    v = v.astype(np.int) / 100
    v_labels = []
    found = False
    for i in v:
        if i > 0 and not found:
            v_labels.append(0)
            found = True
        else:
            v_labels.append("")

    fig, ax = plt.subplots(figsize=(6, 10))
    c = ax.pcolor(words, cmap=CMAP)
    ax.set_title("Ordered by log(%CLlabel/%PLlabel) S/W")

    plt.xticks(np.arange(len(words[0])), get_ticks(words[0]), **{"ha": "center", "va": "center"})
    # plt.xlabel("A:b, C:o, G:r, T:g")
    plt.xlabel("W:b, S:o")

    plt.yticks(np.arange(len(words)), v_labels, **{"va": "bottom"})

    fig.tight_layout()
    plt.savefig(r"sequence_dist\flank3.png")


if __name__ == '__main__':
    create_info()
    # plot_info()
