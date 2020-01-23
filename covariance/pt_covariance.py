import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import glob
import os
import re
import sys

# sys.path.append(os.path.dirname(os.getcwd()))

BEDGRAPH_LINE_FORMAT = "{chr_name}\t{start}\t{end}\t{number}\n"
BEDGRPH_OUTPUT_FILE = "covariance_between_CpG_%s.bedgraph"
BEDGRAPH_OUTPUT_FILE_FORMAT = "average_covariance_between_cpg_%s_chr_%s_group_%s.bedgraph"
WINDOWS_SIZE = 5000

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_chr(\d+).dummy.pkl.zip")

GROUPS = ['PT']


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_group_bedgraph(file, patient, chromosome, group, output):
    to_read = file
    df = pd.read_pickle(to_read)
    pt_index = [pt for pt in df.index if pt.startswith(group)]
    pt_df = df.loc[pt_index, :]
    num_of_cpg = 10000
    # num_of_cpg = pt_df.shape[1]

    with open(os.path.join(output, BEDGRAPH_OUTPUT_FILE_FORMAT % (patient, chromosome, group)), "w") as output_file:
        for i in tqdm(range(0, num_of_cpg, WINDOWS_SIZE)):
            pmd_columns = pt_df.columns[i:min(i + WINDOWS_SIZE, num_of_cpg)]
            covariance = df.loc[:, pmd_columns].cov()
            average_covariance = covariance.mean()
            for cpg in pmd_columns:
                start = cpg
                end = cpg + 1
                number = average_covariance[cpg]
                if not np.isnan(number):
                    line = BEDGRAPH_LINE_FORMAT.format(chr_name=chromosome, start=start, end=end,
                                                       number=number)
                    output_file.write(line)


def create_bedgraphs(file, patient, chromosome, output):
    for group in GROUPS:
        create_group_bedgraph(file, patient, chromosome, group, output)


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    chr = args.chr

    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        if chr:
            cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % chr)

        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)
    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    for file in all_cpg_format_file_paths:
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file)[0]
        create_bedgraphs(file, patient, chromosome, output)


if __name__ == '__main__':
    main()
