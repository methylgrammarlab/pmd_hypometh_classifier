# This was added only because we don't have pyfaidx on the cluster
import os
import sys

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from zhou_dataset.get_valid_cpg import *


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', help='Path for input df', required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_input()
    input_file = args.input_file
    df = pd.read_pickle(input_file)

    output_folder = os.path.dirname(input_file)
    output_file = os.path.join(output_folder, "context_dict.pkl")

    chr_dict = {}
    for chromosome in pd.unique(df["chromosome"]):
        temp_df = df[df["chromosome"] == chromosome]
        chromosome_dict = temp_df[["location", "sequence"]].set_index("location").to_dict()
        chr_dict[chromosome] = chromosome_dict["sequence"]

    with open(output_file, "wb") as f:
        pickle.dump(chr_dict, f)


if __name__ == '__main__':
    main()
