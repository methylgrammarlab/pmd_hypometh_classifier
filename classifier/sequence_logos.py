from weblogo import *
import seqlogo
import logomaker as lm
import argparse
import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools

PL = 0
CL = 1


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequence_file', help='path for the sequence prediction file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))

    args = parser.parse_args()
    return args


def get_data(score, sequence):
    ind = sequence.find('CG')
    pl_cl = PL if float(score) < 0.5 else CL
    minus_2 = sequence[ind - 1]
    minus_1 = sequence[ind - 1]
    plus_1 = sequence[ind + 2]
    plus_2 = sequence[ind + 3]
    return pl_cl, minus_2, minus_1, plus_1, plus_2, sequence


def get_instances(sequence_df, pl_cl, minus_2_op=('A', 'C', 'G', 'T'), minus_1_op=('A', 'C', 'G', 'T'),
                  plus_1_op=('A', 'C', 'G', 'T'), plus_2_op=('A', 'C', 'G', 'T')):
    instances = sequence_df.loc[
        (sequence_df['PL_CL'] == pl_cl) &
        (sequence_df['minus_2'].isin(minus_2_op)) &
        (sequence_df['minus_1'].isin(minus_1_op)) &
        (sequence_df['plus_1'].isin(plus_1_op)) &
        (sequence_df['plus_2'].isin(plus_2_op)),
        'sequence']
    return instances


def create_seq_data(seq_list):
    seqs = SeqList(seq_list)
    seqs.alphabet = Alphabet('ACGT')
    return seqs


def create_ppm(seq_df):
    ppm = seq_df.apply(list).apply(pd.Series).apply(pd.Series.value_counts).T
    return ppm.fillna(0) / seq_df.shape[0]


def plot_logomaker(seq_list):
    counts_mat = lm.alignment_to_matrix(list(seq_list)).iloc[50:100, :]
    lm.Logo(counts_mat)
    plt.show()
    info_mat = lm.transform_matrix(counts_mat,
                                   from_type='counts',
                                   to_type='information')
    lm.Logo(info_mat)
    plt.show()
    prob_mat = lm.transform_matrix(counts_mat,
                                   from_type='counts',
                                   to_type='probability')
    lm.Logo(prob_mat)
    plt.show()


def plot_seqlogo(seq_list):
    # This is probability actualy :(
    summary_ppm = create_ppm(seq_list)
    ppm = seqlogo.Ppm(summary_ppm)
    seqlogo.seqlogo(ppm, ic_scale=False, format='eps', size='xlarge', filename='seqlogo_logo.eps')


def plot_weblogo(seq_list):
    seqs = create_seq_data(seq_list)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "A Logo Title"
    logoformat = LogoFormat(logodata, logooptions)
    eps = eps_formatter(logodata, logoformat)
    with open('weblogo_logo.eps', "wb") as f:
        f.write(eps)


def main():
    args = parse_input()
    sequences = files_tools.load_compressed_pickle_not_zlib(args.sequence_file)

    sequence_df = pd.DataFrame(columns=['PL_CL', 'minus_2', 'minus_1', 'plus_1', 'plus_2', 'sequence'],
                               index=range(sequences.shape[1]))
    for i in tqdm(range(sequences.shape[1])):
        sequence_df.iloc[i, :] = get_data(sequences[0, i], sequences[1, i])

    instances = get_instances(sequence_df, PL, minus_1_op=('C', 'G'), plus_1_op=('C', 'G'))

    seq_list = list(instances)
    plot_logomaker(seq_list)
    plot_seqlogo(instances)
    plot_weblogo(seq_list)


if __name__ == '__main__':
    main()

# instances = df['binding'] #just input the list of DNA sequences
# m = motifs.create(instances)
# m.weblogo('logo.png')
