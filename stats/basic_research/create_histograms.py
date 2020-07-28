# TODO: lior do you remember why ?


import argparse
import glob
import os
import re
import sys

from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools, data_tools

DICT_FILE_FORMAT = "dict_*_%s.pickle"
DICT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).pickle")
TITLE_FORMAT = "histogram %s %s"
HISTOGRAM_FORMAT = "histogram_%s_%s"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dict_files', help='Path to folder or file of dictionary', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_histogram(file, output):
    patient, chromosome = DICT_FILE_RE.findall(file)[0]
    chromosome_dict = files_tools.load_pickle(file)
    data_tools.dict_to_histogram(chromosome_dict, xlabel='number of cells covering both location pairs',
                                 title=TITLE_FORMAT % (patient, chromosome),
                                 save_path=os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome)))


def main():
    args = parse_input()

    output_folder = args.output_folder
    chromosome = args.chr

    if os.path.isdir(args.dict_files):
        cpg_format_file_path = os.path.join(args.dict_files, DICT_FILE_FORMAT % '*')
        if chromosome:
            cpg_format_file_path = os.path.join(args.dict_files, DICT_FILE_FORMAT % chromosome)

        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)
    else:
        all_cpg_format_file_paths = [args.dict_files]

    for file in tqdm(all_cpg_format_file_paths, desc='files'):
        create_histogram(file, output_folder)


if __name__ == '__main__':
    main()
