import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.getcwd())
from format_files import format_sublineage_info
from commons import consts

BEDGRAPH_LINE_FORMAT = "{chr_name}\t{start}\t{end}\t{number}\n"
BEDGRPH_OUTPUT_FILE = "covariance_between_CpG_%s.bedgraph"
BEDGRAPH_OUTPUT_FILE_FORMAT = "average_covariance_between_cpg_%s_chr_%s_region_%s.bedgraph"
WINDOWS_SIZE = 500

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_chr(\d+).dummy.pkl.zip")

ALL = 'ALL'
REGIONS = [ALL, 'PT']


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_region_bedgraph(file, patient, chromosome, region, region_name, output):
    """
    Creates a bedgraph of the mean covariance of cPg's, using data only from one region at a time, adding the NC. The
    covariance is calculated in windows.
    :param file: The filename of the file with the parsed scWGBS data
    :param patient: Patient ID of the current data
    :param chromosome: The current chromosome
    :param region: The region ID
    :param output: The output folder
    """
    df = pd.read_pickle(file)
    if region is not ALL:
        region_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]

        for sample in region:
            region_cell_ids.extend(cell_id for cell_id in df.index if cell_id.startswith(sample))

        region_df = df.loc[region_cell_ids, :]
    else:
        region_df = df
    # num_of_cpg = 100
    num_of_cpg = region_df.shape[1]
    output_filename = os.path.join(output, BEDGRAPH_OUTPUT_FILE_FORMAT % (patient, chromosome, 'NCand%s' %
                                                                          region_name))

    # TODO covariance with min_periods=10
    with open(output_filename, "w") as output_file:
        for i in tqdm(range(0, num_of_cpg, WINDOWS_SIZE)):
            cur_columns_inds = region_df.columns[i:min(i + WINDOWS_SIZE, num_of_cpg)]
            covariance = region_df.loc[:, cur_columns_inds].cov()
            average_covariance = covariance.mean()
            for cpg in cur_columns_inds:
                start = cpg
                end = cpg + 1
                number = average_covariance[cpg]
                if not np.isnan(number):
                    line = BEDGRAPH_LINE_FORMAT.format(chr_name=chromosome, start=start, end=end,
                                                       number=number)
                    output_file.write(line)


def create_bedgraphs(file, patient, chromosome, output):
    """
    Goes over all the regions and calls the function to create a bedgraph for each one.
    :param file: The filename of the file with the parsed scWGBS data
    :param patient: Patient ID of the current data
    :param chromosome: The current chromosome
    :param output: The output folder
    """
    sublineage_info = format_sublineage_info.get_sublineage_info(consts.SUBLINEAGE_FILE_LOCAL_DROR)
    patient_info = sublineage_info[patient]
    for sublineage in patient_info:
        create_region_bedgraph(file, patient, chromosome, patient_info[sublineage], sublineage, output)


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
