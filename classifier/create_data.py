import argparse
import glob
import os
import re
import sys

import pandas as pd
import pyfaidx

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import consts

COVARIANCE_INFO_FILE_FORMAT = "covariance_info_chr_*.pkl"
COVARIANCE_CHR_FORMAT_RE = re.compile("covariance_info_chr_(\d{1,2}).pkl")
OUTPUT_FILE = "classifier_data_chr%s.pkl"

genome = pyfaidx.Fasta(consts.GENOME_FILE_LOCAL_DROR, sequence_always_upper=True, as_raw=True)


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--covariance_info_folder', help='path for the covariance info folder', required=True)
    parser.add_argument('--seq_size', help="the length of the sequence", required=False, default=150,
                        type=int)

    args = parser.parse_args()
    return args


def update_file(cov_file, output_folder, seq_size):
    chr_num = COVARIANCE_CHR_FORMAT_RE.findall(cov_file)[0]
    output_path = os.path.join(output_folder, OUTPUT_FILE % chr_num)
    cov_data = pd.read_pickle(cov_file)
    seq = []

    for i in cov_data["cpg_index"]:
        seq = get_seq_for_cpg(chr_num, genome, i, seq_size)
        seq.append(seq)

    cov_data["seq"] = seq
    cov_data.to_pickle(output_path)


def get_seq_for_cpg(chr_num, i, seq_size):
    chr_info = genome[chr_num]
    # assert chr_info[i] == "G"
    seq_size = int((seq_size - 2) / 2)
    return chr_info[i - seq_size - 1:i + seq_size + 1]


def main():
    args = parse_input()

    covariance_files = glob.glob(os.path.join(args.covariance_info_folder, COVARIANCE_INFO_FILE_FORMAT))

    for cov_file in covariance_files:
        update_file(cov_file, args.output_folder, args.seq_size)


if __name__ == '__main__':
    main()
