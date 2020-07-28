"""
Remove CpG from a dataframe which match some kind of blacklist of CpG island file
We used blacklist file: hg19-blacklist.v2.bed
We used CpG island file: Irizarry2009-model-based-cpg-islands-hg19.bed
"""

import argparse
import glob
import os
import sys

import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
import commons.files_tools as tools
from commons import consts

CPGI_CHR_INDEX = 0
CPGI_START_INDEX = 1
CPGI_END_INDEX = 2


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', help='Path to files', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=True)
    parser.add_argument('--blacklist', help='Path to the blacklist file', required=False)
    parser.add_argument('--cpgi', help='Path to cpg island file', required=False)
    args = parser.parse_args()

    input_files = args.files
    output_folder = args.output_folder

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


def skim_cpg(cpg_files_to_filter, boundaries_paths, output_folder, data_file_re):
    """
    Skim the cpg in all the files based on the boundaries paths
    :param cpg_files_to_filter: The files to work on
    :param boundaries_paths: list of boundaries files 
    :param output_folder: The output folder 
    :param data_file_re: The regex for the cpg files 
    """
    all_cpg_locations = tools.get_cpg_context_map(only_locations=True)

    # Get a mapping of chr and sites to remove
    cpg_boundaries_dict = get_boundaries(boundaries_paths)
    cpg_locations_mask = get_valid_location_per_chr(all_cpg_locations, cpg_boundaries_dict)

    for cpg_format_file in tqdm(cpg_files_to_filter):
        df = pd.read_pickle(cpg_format_file)

        if data_file_re == consts.DATA_FILE_SCWGBS_RE:
            patient, chromosome = data_file_re.findall(cpg_format_file)[0]
            output_path = get_output_path_for_crc(cpg_format_file, output_folder)

        else:
            chromosome = data_file_re.findall(cpg_format_file)[0]
            output_path = os.path.join(output_folder, os.path.basename(cpg_format_file))

        mask = cpg_locations_mask[chromosome]
        pd.to_pickle(df[mask], output_path)


def get_valid_location_per_chr(all_cpg_locations, cpg_boundaries_dict):
    """
    Get valid location per chromosome
    :param all_cpg_locations: All the cpg we want to look at
    :param cpg_boundaries_dict: The CpG boundaries dict
    :return: A dict of all valid locations after removing the boundaries
    """
    cpg_locations_mask = {}
    # Remove locations based on the boundaries
    for chromosome in tqdm(all_cpg_locations):
        chr_location = all_cpg_locations[chromosome]
        for boundary in cpg_boundaries_dict[chromosome]:
            chr_location = chr_location[(chr_location < boundary[0]) | (chr_location > boundary[1])]

        cpg_locations_mask[chromosome] = chr_location
    return cpg_locations_mask


def get_boundaries(boundaries_file_paths):
    """
    Get the boundaries indexes from teh files
    :param boundaries_file_paths: List with paths of the boundaries files
    :return: The dict with information per chromsome about the boundaries
    """
    cpgi_boundaries_dict = {}
    for boundaries_path in boundaries_file_paths:
        with open(boundaries_path, "r") as boundaries_file:
            for line in tqdm(boundaries_file.readlines()):
                s_line = line.split()
                chr = s_line[CPGI_CHR_INDEX]
                start = int(s_line[CPGI_START_INDEX])
                end = int(s_line[CPGI_END_INDEX])

                if chr not in cpgi_boundaries_dict:
                    cpgi_boundaries_dict[chr] = []

                cpgi_boundaries_dict[chr].append((start, end))
    return cpgi_boundaries_dict


def get_output_path_for_crc(cpg_format_file, output_folder):
    output_path = os.path.join(output_folder, os.path.basename(os.path.dirname(cpg_format_file)),
                               os.path.basename(cpg_format_file))
    if not os.path.exists(os.path.dirname(output_path)):
        os.mkdir(os.path.dirname(output_path))
    return output_path


def convert_sc_files():
    input_files, output_folder, boundaries_paths = format_args()

    all_cpg_format_file_paths = glob.glob(os.path.join(input_files, "CRC*", "*.dummy.pkl.zip"))
    skim_cpg(all_cpg_format_file_paths, boundaries_paths, output_folder,
             data_file_re=consts.DATA_FILE_SCWGBS_RE)


def convert_bulk_files():
    input_files, output_folder, boundaries_paths = format_args()

    all_cpg_format_file_paths = glob.glob(os.path.join(input_files, "*.dummy.pkl.zip"))
    skim_cpg(all_cpg_format_file_paths, boundaries_paths, output_folder,
             data_file_re=consts.DATA_FILE_BULK_RE)


if __name__ == '__main__':
    # convert_sc_files()
    convert_bulk_files()
