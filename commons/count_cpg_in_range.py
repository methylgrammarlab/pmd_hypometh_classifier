import argparse
import os
import sys

import numpy as np

import commons.files_tools as tools

sys.path.append(os.path.dirname(os.getcwd()))


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chr', help='name of chr - just number', required=True)
    parser.add_argument('--start', help='start position', required=True, type=int)
    parser.add_argument('--end', help='end position', required=True, type=int)
    args = parser.parse_args()

    return args.chr, args.start, args.end


def count_cpg_in_range(chr_name, start, end):
    chr_locations_dict = tools.get_cpg_context_map()
    chr_locations = chr_locations_dict[chr_name]
    print(np.where(np.logical_and(chr_locations > start, chr_locations < end))[0].size)


if __name__ == '__main__':
    chr_name, start, end = format_args()
    count_cpg_in_range(chr_name, start, end)
