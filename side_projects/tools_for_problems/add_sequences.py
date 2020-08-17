"""
Allow to append sequence to a dataframe if the sequence failed
"""
import argparse
import os
import sys

import pandas as pd
import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import sequence_tools


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', help='Path for input df', required=True)
    args = parser.parse_args()
    return args


def get_seq_info(df):
    sequences = []

    for i in tqdm.trange(df.shape[0]):
        row = df.iloc[i]
        seq = sequence_tools.get_cpg_sequence(chr_num=row["chromosome"][3:], cpg_index=row["location"])
        sequences.append(seq)

    return sequences


def main():
    args = parse_input()
    input_file = args.input_file

    df = pd.read_pickle(input_file)

    output_folder = os.path.dirname(input_file)
    output_file = os.path.basename(input_file) + "with_seq.pkl"

    df["sequence"] = get_seq_info(df)
    try:
        df.to_pickle(os.path.join(output_folder, output_file))
    except Exception as ex:
        import pdb;
        pdb.set_trace()


if __name__ == '__main__':
    main()
