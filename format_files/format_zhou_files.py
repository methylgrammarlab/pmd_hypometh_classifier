import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
import commons.files_tools as tools
import commons.consts as consts

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
    return data[data["reads"] >= 5]


def main():
    output, files_path = format_args()
    files_to_process = glob.glob(files_path)
    all_cpg_locations = tools.get_cpg_context_map(only_locations=True,
                                                  load_with_path=consts.CONTEXT_MAP_FULL)

    chr_data_dict = {}
    for chr in all_cpg_locations:
        chr_full_cpg = all_cpg_locations[chr]
        chr_all_cells = np.full((len(files_to_process) + 1, all_cpg_locations[chr].size), None,
                                dtype=np.float)
        chr_all_cells[0, :] = chr_full_cpg
        chr_data_dict[chr] = chr_all_cells

    i = 1
    cell_names = []
    for f in files_to_process:
        sample_name = FILE_RE.findall(os.path.basename(f))[0]
        try:
            data = get_data_from_file(f)
        except:
            print(sample_name)
        cell_names.append(sample_name)

        for chr in chr_data_dict:
            chr_full_cpg = all_cpg_locations[chr]
            chr_data = data[data["chr"] == chr]
            match = np.isin(chr_full_cpg, chr_data["start"])
            chr_data_dict[chr][i, match] = chr_data["ratio"]

        i += 1

    for chr in chr_data_dict:
        chr_all_cells = chr_data_dict[chr]
        df = pd.DataFrame(data=chr_all_cells[1:, :], columns=chr_all_cells[0, :].astype(np.int),
                          index=cell_names, dtype=np.float16)
        output_path = os.path.join(output, OUTPUT_FILE_FORMAT % chr)
        df.to_pickle(output_path, compression='zip')


if __name__ == '__main__':
    main()
