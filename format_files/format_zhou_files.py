"""
Format the zhou samples for the pipelines for files:
- One file for each chromosome, each row is a sample and each column is a CpG
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
import commons.files_tools as tools
import commons.consts as consts

INPUT_FILE_GLOB_FORMAT = "*.bedgraph.gz"
INPUT_FILE_INFO_RE = re.compile("(.*)\.bedgraph\.gz")

MIN_READS = 5

def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', help='Path to raw files of scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--min_reads', help='min read for CpG', required=False, default=MIN_READS)

    args = parser.parse_args()

    files_path = os.path.join(args.raw, INPUT_FILE_GLOB_FORMAT)

    return args.output_folder, files_path, args.min_reads


def get_data_from_samples(bedgraph_path, min_reads=MIN_READS):
    """
    Extract information from a file, only take data from places we had more then min_reads
    :param bedgraph_path: The path for the file
    :param min_reads: the minimum number of reads for a cpg to include it
    :return: The bedgraph data - chromosome, cpg index and ratio of methylated
    """
    data = pd.read_csv(bedgraph_path, sep="\t", names=["chr", "location", "end", "ratio", "reads"],
                       usecols=["chr", "location", "ratio", "reads"])
    data["location"] = data["location"] + 1  # to get the CpG index
    return data[data["reads"] >= min_reads]


def main():
    output, files_path, min_reads = format_args()
    files_to_process = glob.glob(files_path)
    all_cpg_locations = tools.get_cpg_context_map(only_locations=True,
                                                  load_with_path=consts.CONTEXT_MAP_FULL)

    # Create a dictionary for chromosomse, each value is a np array with the size of the chromosome CpG x
    # num of samples
    chr_data_dict = {}
    for chromosome in all_cpg_locations:
        chr_full_cpg = all_cpg_locations[chromosome]
        chr_all_cells = np.full((len(files_to_process) + 1, all_cpg_locations[chromosome].size), None,
                                dtype=np.float)
        chr_all_cells[0, :] = chr_full_cpg
        chr_data_dict[chromosome] = chr_all_cells

    i = 1
    cell_names = []
    for f in files_to_process:
        sample_name = INPUT_FILE_INFO_RE.findall(os.path.basename(f))[0]
        cell_names.append(sample_name)

        data = get_data_from_samples(f)

        for chromosome in chr_data_dict:
            chr_full_cpg = all_cpg_locations[chromosome]
            chr_data = data[data["chromosome"] == chromosome]
            match = np.isin(chr_full_cpg, chr_data["location"])
            chr_data_dict[chromosome][i, match] = chr_data["ratio"]

        i += 1

    # Convert each np to pandas data frame
    for chromosome in chr_data_dict:
        chr_all_cells = chr_data_dict[chromosome]
        df = pd.DataFrame(data=chr_all_cells[1:, :], columns=chr_all_cells[0, :].astype(np.int),
                          index=cell_names, dtype=np.float16)
        output_path = os.path.join(output, consts.BULK_FILE_FORMAT % chromosome)
        df.to_pickle(output_path, compression='zip')


if __name__ == '__main__':
    main()
