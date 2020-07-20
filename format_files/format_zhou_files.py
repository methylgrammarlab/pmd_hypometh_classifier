import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
import commons.files_tools as tools

OUTPUT_FILE_FORMAT = "all_cpg_ratios_%s_hg19.dummy.pkl.zip"
FILE_SUFFIX = "*.bedgraph.gz"
FILE_RE = re.compile("(.*)\.bedgraph\.gz")


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)

    args = parser.parse_args()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    files_path = os.path.join(args.raw, FILE_SUFFIX)

    return output, files_path


def get_data_from_file(bedgraph_path):
    data = pd.read_csv(bedgraph_path, sep="\t", names=["chr", "start", "end", "ratio", "reads"],
                       usecols=["chr", "start", "ratio", "reads"])
    data["start"] = data["start"] + 1
    return data


def main():
    output, files_path = format_args()
    files_to_process = glob.glob(files_path)
    all_cpg_locations = tools.get_cpg_context_map(only_locations=True)

    chr_data = {}
    for chr in all_cpg_locations:
        chr_data[chr] = []

    for f in tqdm.tqdm(files_to_process, desc="Reading files"):
        data = get_data_from_file(f)
        sample_name = FILE_RE.findall(os.path.basename(f))[0]

        for chr in all_cpg_locations:
            chr_data[chr].append((sample_name, data[data["chr"] == chr]))

    for chr in tqdm.tqdm(all_cpg_locations, desc="Creating df"):
        chr_full_cpg = all_cpg_locations[chr]

        chr_all_cells = np.full((len(chr_data[chr]) + 1, chr_full_cpg.size), None, dtype=np.float)
        chr_all_cells[0, :] = chr_full_cpg
        cell_names = []

        i = 1
        for sample in chr_data[chr]:
            cell_names.append(sample[0])
            cell = sample[1]
            match = np.isin(chr_full_cpg, cell["start"])
            chr_all_cells[i, match] = cell["ratio"]

            i += 1

        df = pd.DataFrame(data=chr_all_cells[1:, :], columns=chr_all_cells[0, :].astype(np.int),
                          index=cell_names,
                          dtype=np.float16)
        output_path = os.path.join(output, OUTPUT_FILE_FORMAT % chr)
        df.to_pickle(output_path, compression='zip')


if __name__ == '__main__':
    main()
