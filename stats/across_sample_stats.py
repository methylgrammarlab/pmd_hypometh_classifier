import argparse
import collections
import glob
import os
import re
import sys

from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

import numpy as np
import pandas as pd

FOLDER_SUFFIX = os.path.join("*", "all_cpg_ratios_CRC*_chr*.dummy.pkl.zip")
CHR_RE = re.compile(".*chr(\d+).*")
COLUMNS = ["read_number", "counter"]


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()

    return args.files_folder, args.output_folder


def main():
    input_folder, output_folder = format_args()
    all_files = glob.glob(os.path.join(input_folder, FOLDER_SUFFIX))
    chr_cpg_dict = {}
    not_nans_dict = {}

    covered_counter = collections.Counter()

    for file_path in tqdm(all_files):
        chr = CHR_RE.findall(file_path)[0]
        data = pd.read_pickle(file_path)
        not_nans = data.count()

        if chr not in chr_cpg_dict:
            chr_cpg_dict[chr] = not_nans.index

        if chr not in not_nans_dict:
            not_nans_dict[chr] = []

        not_nans_dict[chr].append(not_nans.values)
        covered_counter.update(not_nans.values / data.shape[0])

    # Combine all the chr to one
    for chr in chr_cpg_dict:
        table = np.array(not_nans_dict[chr])
        df = pd.DataFrame(data=table.astype(np.int), columns=chr_cpg_dict[chr].astype(np.int))
        df.to_csv(os.path.join(output_folder, "chr%s_full_mapping.csv" % chr))

    counter_df = pd.DataFrame.from_dict(covered_counter, orient='index').reset_index()
    counter_df.columns = COLUMNS
    counter_df.to_csv(os.path.join(output_folder, "covered_samples_counter.csv"), "w")


if __name__ == '__main__':
    main()
