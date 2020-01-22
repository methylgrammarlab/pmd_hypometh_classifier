import argparse
import glob
import os
import re
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
import commons.files_tools as tools

CPGI_CHR_INDEX = 0
CPGI_START_INDEX = 1
CPGI_END_INDEX = 2
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', help='Path to files', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=True)
    parser.add_argument('--blacklist', help='Path to the blacklist file', required=False)
    parser.add_argument('--cpgi', help='Path to cpg island file', required=False)
    args = parser.parse_args()

    input_files = args.files
    if not os.path.exists(input_files):
        print("Couldn't find the input path")
        sys.exit(-1)

    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        print("Couldn't find the output path")
        sys.exit(-1)

    boundaries_paths = []
    blacklist_path = args.blacklist
    cpgi_path = args.cpgi

    if cpgi_path:
        if not os.path.exists(cpgi_path):
            print("Coludn't find the cpgi path")
            sys.exit(-1)
        boundaries_paths.append(cpgi_path)

    if blacklist_path:
        if not os.path.exists(blacklist_path):
            print("Coludn't find the blacklist path")
            sys.exit(-1)
        boundaries_paths.append(blacklist_path)

    return input_files, output_folder, boundaries_paths


def skim_cpg(files, boundaries_paths, output_folder):
    # Get a mapping of chr and sites to remove
    cpg_locations_mask = {}
    cpgi_boundaries_dict = {}
    for boundaries_path in boundaries_paths:
        with open(boundaries_path, "r") as boundaries_file:
            for line in tqdm(boundaries_file.readlines(), desc="CpGI file: %s" % os.path.basename(boundaries_path)):
                s_line = line.split()
                chr = s_line[CPGI_CHR_INDEX]
                start = int(s_line[CPGI_START_INDEX])
                end = int(s_line[CPGI_END_INDEX])

                if chr not in cpgi_boundaries_dict:
                    cpgi_boundaries_dict[chr] = []

                cpgi_boundaries_dict[chr].append((start, end))

    all_cpg_locations = tools.get_cpg_context_map(only_locations=True)
    for chromosome in tqdm(all_cpg_locations, desc="CpGI create chr mask"):
        chr_location = all_cpg_locations[chromosome]
        for boundary in cpgi_boundaries_dict[chromosome]:
            chr_location = chr_location[(chr_location < boundary[0]) | (chr_location > boundary[1])]

        cpg_locations_mask[chromosome] = chr_location

    for cpg_format_file in tqdm(files, desc="CpGI convert"):
        df = pd.read_pickle(cpg_format_file)
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(cpg_format_file)[0]
        mask = cpg_locations_mask[chromosome]

        output_path = get_output_path_for_crc(cpg_format_file, output_folder)
        pd.to_pickle(df[mask], output_path)


def get_output_path_for_crc(cpg_format_file, output_folder):
    output_path = os.path.join(output_folder, os.path.basename(os.path.dirname(cpg_format_file)),
                               os.path.basename(cpg_format_file))
    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    return output_path


def main():
    input_files, output_folder, boundaries_paths = format_args()

    all_cpg_format_file_paths = glob.glob(os.path.join(input_files, "CRC*", "*.dummy.pkl.zip"))
    skim_cpg(all_cpg_format_file_paths, boundaries_paths, output_folder)


if __name__ == '__main__':
    main()
