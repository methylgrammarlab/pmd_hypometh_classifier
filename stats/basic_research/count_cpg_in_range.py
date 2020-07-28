"""
Count the amount of CpG in a specific range for a specific chromosome
This was used to get a grip on teh amount of CpG we will have in each window and thorughout the work
"""


import argparse
import os
import sys

import numpy as np

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
import commons.files_tools as tools


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chr', help='name of chr - just number', required=True)
    parser.add_argument('--start', help='start position', required=True, type=int)
    parser.add_argument('--end', help='end position', required=True, type=int)
    args = parser.parse_args()

    return args.chr, args.start, args.end


def count_cpg_in_range(chromosome_num, start_index, end_index):
    """
    Count the amount of CpG in a specifici range
    :param chromosome_num: The chromosome number
    :param start_index: The start index to check
    :param end_index: The end index to check
    :return:
    """
    chr_locations_dict = tools.get_cpg_context_map(only_locations=True)
    chr_locations = chr_locations_dict["chr%s" % chromosome_num]
    print(np.where(np.logical_and(chr_locations > start_index, chr_locations < end_index))[0].size)


if __name__ == '__main__':
    chr_name, start, end = format_args()
    count_cpg_in_range(chr_name, start, end)
