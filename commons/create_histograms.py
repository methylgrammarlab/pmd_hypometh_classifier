import argparse
import glob
import matplotlib.pyplot as plt
import re
import sys
import os

sys.path.append(os.path.dirname(os.getcwd()))
from commons import files_tools

DICT_FILE_FORMAT = "dict_*_%s.pickle"  # TODO remove the mini
DICT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).pickle")
TITLE_FORMAT = "histogram %s %s"
HISTOGRAM_FORMAT = "histogram_%s_%s"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dict_files', help='Path to folder or file of dictionary', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


# def create_histogram(series, patient, chromosome, num_of_bins, output):
def create_histogram(file, output):
    patient, chromosome = DICT_FILE_RE.findall(file)[0]
    dict = files_tools.load_compressed_pickle_not_zlib(file)
    plt.bar(list(dict.keys()), dict.values(), color='#0E74E3')
    plt.xlabel('number of cells covering both location pairs')
    plt.title(TITLE_FORMAT % (patient, chromosome))
    plt.savefig(os.path.join(output, HISTOGRAM_FORMAT % (patient, chromosome)))
    plt.close()


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    chr = args.chr

    if os.path.isdir(args.dict_files):
        cpg_format_file_path = os.path.join(args.dict_files, DICT_FILE_FORMAT % '*')
        if chr:
            cpg_format_file_path = os.path.join(args.dict_files, DICT_FILE_FORMAT % chr)

        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)
    else:
        all_cpg_format_file_paths = [args.dict_files]

    # for file in tqdm(all_cpg_format_file_paths, desc='files'):
    for file in all_cpg_format_file_paths:
        create_histogram(file, output)


if __name__ == '__main__':
    main()
