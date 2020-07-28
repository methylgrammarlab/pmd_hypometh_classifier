# This was added only because we don't have pyfaidx on the cluster
import os
import sys

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from zhou_dataset.get_valid_cpg import *
import tqdm


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', help='Path foir input df', required=True)
    args = parser.parse_args()
    return args


def get_seq_info(df, ):
    sequences = []

    for i in tqdm.trange(df.shape[0]):
        row = df.iloc[i]
        seq = get_seq_for_cpg(chr_num=row["chromosome"], i=row["location"], seq_size=SEQ_SIZE)
        sequences.append(seq)

    return sequences


def add_chr_number(df):
    df["chr_num"] = 0
    for chromosome in pd.unique(df["chromosome"]):
        df.loc[df["chromosome"] == chromosome, "chr_num"] = INT_RE.findall(chromosome)[0]

    return df


def main():
    args = parse_input()
    input_file = args.input_file

    df = pd.read_pickle(input_file)
    dff = add_chr_number(df)

    output_folder = os.path.dirname(input_file)
    output_file = os.path.basename(input_file) + "with_seq.pkl"

    df["sequence"] = get_seq_info(df)
    df.to_pickle(os.path.join(output_folder, output_file))


if __name__ == '__main__':
    main()
