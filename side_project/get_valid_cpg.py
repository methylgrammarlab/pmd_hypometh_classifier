# extract all chr16 cpg - without CpGI

CHR_PATH = r"H:\Study\university\Computational-Biology\Year 3\Projects\proj_scwgbs\resource\genomic_data\filtered_by_bl_and_cpgi\full_cpg_seq_chr16.hg19.pickle.zlib"

import os
import re
import sys

import pandas as pd
import pyfaidx

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import files_tools, consts
from variance import get_valid_cpg

MEAN = 1
VAR = 2
METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

SEQ_SIZE = 150

genome = pyfaidx.Fasta(consts.GENOME_FILE_LOCAL_DROR, sequence_always_upper=True, as_raw=True)


def main():
    chr_data = files_tools.load_compressed_pickle(CHR_PATH)
    all_cpg = chr_data[:, 0]
    sequences = get_valid_cpg.get_seq_info(all_cpg, '16')

    sum_df = pd.DataFrame(columns=["location"])
    sum_df["location"] = all_cpg
    sum_df["sequence"] = sequences
    sum_df.to_pickle("side_project_chr16_cpg.pkl")


if __name__ == '__main__':
    main()
