"""
Get information about the converage across the genome to csv and bedgraph
"""

import argparse
import collections
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from tqdm import tqdm

matplotlib.use("Agg")

# Needed for imports
sys.path.append(os.path.dirname(os.getcwd()))
from commons import files_tools, consts, data_tools

BEDGRAPH_LINE_FORMAT = "{chr_name}\t{start}\t{end}\t{number}\n"
BEDGRPH_OUTPUT_FILE = "cover_across_samples_%s.bedgraph"


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))

    return parser


def get_file_info(file_path):
    """
    Read the file, extract the information about the patient and chr by the name
    :param file_path: The path of the file
    :return: The relevant data, the pateint name and the chromosome name
    """
    patient, chromosome = consts.PATIENT_CHR_NAME_RE(file_path)[0]

    data = pd.read_pickle(file_path)
    return data, patient, chromosome


def extract_coverage_global(input_folder, output_folder):
    """
    Coverage of location across chromosome for all patient and all chr
    :param input_folder: The path of all the files
    :param output_folder:  The output path
    """
    global_counter = collections.Counter()

    all_files = files_tools.get_files_to_work(input_folder, consts.SCWGBS_FILE_FORMAT % ("*", "*"))

    for file_path in tqdm(all_files):
        data, patient, chr_name = get_file_info(file_path)
        coverage_across_samples = data.count()
        global_counter.update(coverage_across_samples.values / data.shape[0])

    commons.data_tools.counter_to_csv(global_counter,
                                      os.path.join(output_folder, "covered_samples_counter.csv"))


def extract_coverage_per_chr(input_folder, output_folder):
    """
    Coverage of location across chromosome for all patient: same as before but trying to figure out if
    there is something different between chromosome, we can also export this information to image
    :param input_folder: The path of all the files
    :param output_folder:  The output path
    """
    chr_cpg_dict = files_tools.get_cpg_context_map(only_locations=True)
    coverage_across_chr_dict = {}

    all_files = files_tools.get_files_to_work(input_folder, consts.SCWGBS_FILE_FORMAT % ("*", "*"))

    for file_path in tqdm(all_files):
        data, patient, chr_name = get_file_info(file_path)
        coverage_across_samples = data.count()

        # Init the coverage_across_chr_dict dict
        if chr_name not in chr_cpg_dict:
            coverage_across_chr_dict[chr_name] = []

        coverage_across_chr_dict[chr_name].append(coverage_across_samples.values / data.shape[0])

    # File per chr - all patients
    for chr_name in chr_cpg_dict:
        table = np.array(coverage_across_chr_dict[chr_name])

        df = pd.DataFrame(data=table.astype(np.int), columns=chr_cpg_dict[chr_name].astype(np.int))
        df.to_csv(os.path.join(output_folder, "%s_full_mapping.csv" % chr_name))


def extract_bedgraph_information(input_folder, output_folder):
    """
    Create a bedgraph file for each patient with the coverage across different samples
    :param input_folder: The folder with all the input files
    :param output_folder: The ouput folder for the bedgraphs
    """
    chr_cpg_dict = files_tools.get_cpg_context_map(only_locations=True)
    all_files = files_tools.get_files_to_work(input_folder, consts.SCWGBS_FILE_FORMAT % ("*", "*"))

    coverage_per_patient_dict = {}

    for file_path in tqdm(all_files):
        data, patient, chr_name = get_file_info(file_path)
        not_nans = data.count()

        if patient not in coverage_per_patient_dict:
            coverage_per_patient_dict[patient] = {}  # Chr dict

        coverage_per_patient_dict[patient][chr_name] = not_nans.values / data.shape[0]

    create_bedgraph_file(chr_cpg_dict, coverage_per_patient_dict, output_folder)


def create_bedgraph_file(chr_cpg_dict, coverage_per_patient_dict, output_folder):
    """
    Create a bedgraph files based on a specific patient
    :param chr_cpg_dict: The position of the CpG
    :param coverage_per_patient_dict: Information about the coverage for each patient
    :param output_folder: The output folder
    """
    for patient in coverage_per_patient_dict:
        output_path = os.path.join(output_folder, BEDGRPH_OUTPUT_FILE % patient)
        with open(output_path, "w") as output_file:
            for chromosome in coverage_per_patient_dict[patient]:
                for i in range(len(coverage_per_patient_dict[patient][chromosome])):
                    start = chr_cpg_dict[chromosome][i]
                    end = start + 1
                    number = coverage_per_patient_dict[patient][chromosome][i]
                    line = BEDGRAPH_LINE_FORMAT.format(chr_name=chromosome, start=start, end=end,
                                                       number=number)
                    output_file.write(line)


def extract_coverage_per_patient(input_folder, output_folder):
    """
    Coverage of locations per patient across all chr: for each location count how many samples contains
    something from this locations. Add this information to counter.
    :param input_folder: The path of all the files
    :param output_folder:  The output path
    """
    coverage_across_patient_dict = {}

    all_files = files_tools.get_files_to_work(input_folder, consts.SCWGBS_FILE_FORMAT % ("*", "*"))

    for file_path in tqdm(all_files):
        data, patient, chr_name = get_file_info(file_path)
        coverage_across_samples = data.count()

        # Init the coverage_across_patient_dict dict
        if patient not in coverage_across_patient_dict:
            coverage_across_patient_dict[patient] = collections.Counter()

        coverage_across_patient_dict[patient].update(coverage_across_samples.values / data.shape[0])

    # Covered of each patient, global on chr
    for patient in coverage_across_patient_dict:
        data_tools.counter_to_csv(coverage_across_patient_dict[patient],
                                  os.path.join(output_folder, "%s_covered_samples_counter.csv" % patient))


def main():
    parser = format_args()
    args = parser.parse_args()
    input_folder, output_folder = args.files_folder, args.output_folder

    extract_bedgraph_information(input_folder, output_folder)
    extract_coverage_per_patient(input_folder, output_folder)
    extract_coverage_per_chr(input_folder, output_folder)
    extract_coverage_global(input_folder, output_folder)


if __name__ == '__main__':
    main()
