import argparse
import collections
import glob
import os
import re
import sys

import matplotlib
import numpy as np
import pandas as pd
from tqdm import tqdm

# Needed to run from shall
import commons.data_tools
import commons.slurm_tools

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Needed for imports
sys.path.append(os.path.dirname(os.getcwd()))
from commons.consts import CENETROMERE_DICT
from commons import files_tools

FOLDER_SUFFIX = os.path.join("*", "all_cpg_ratios_CRC*_chr*.dummy.pkl.zip")

CHR_RE = re.compile(".*chr(\d+).*")
PATIENT_RE = re.compile("CRC\d\d")

BEDGRAPH_LINE_FORMAT = "chr{chr_name}\t{start}\t{end}\t{number}\n"
BEDGRPH_OUTPUT_FILE = "cover_across_samples_%s.bedgraph"


def fix_counter(file_name):
    """
    This should actually be used, only used once to convert the fractions to integars after rounding to
    show them better on excel
    :param file_name:
    :return:
    """
    new_file = file_name.replace(".csv", "_new.csv")
    data = open(file_name, "r").readlines()

    c = collections.Counter()
    for line in data[1:]:
        _, number, counter = line.strip().split("w")
        c.update({round(float(number) * 100): int(counter)})

    commons.data_tools.counter_to_csv(c, new_file)


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--create_chromosome_img', help='Should create chromosome img', required=False,
                        default=False, action='store_true')

    return parser


def get_file_info(file_path):
    """
    Read the file, extract the information about the patient and chr by the name
    :param file_path:
    :return:
    """
    patient = PATIENT_RE.findall(os.path.basename(file_path))[0]
    chr_name = CHR_RE.findall(file_path)[0]

    data = pd.read_pickle(file_path)
    return data, patient, chr_name


def extract_coverage_global(input_folder, output_folder):
    """
    Coverage of location across chromosome for all patient and all chr
    """
    chr_cpg_dict = files_tools.get_all_cpg_locations_across_chr()
    global_counter = collections.Counter()

    all_files = glob.glob(os.path.join(input_folder, FOLDER_SUFFIX))

    for file_path in tqdm(all_files):
        data, patient, chr_name = get_file_info(file_path)
        coverage_across_samples = data.count()
        global_counter.update(coverage_across_samples.values / data.shape[0])

    commons.data_tools.counter_to_csv(global_counter, os.path.join(output_folder, "covered_samples_counter.csv"))


def extract_coverage_per_chr(input_folder, output_folder):
    """
    Coverage of location across chromosome for all patient: same as before but trying to figure out if
    there is something different between chromosome, we can also export this information to image
    """
    parser = format_args()
    args = parser.parse_args()
    create_chromosome_img = args.create_chromosome_img

    chr_cpg_dict = files_tools.get_all_cpg_locations_across_chr()
    coverage_across_chr_dict = {}

    all_files = glob.glob(os.path.join(input_folder, FOLDER_SUFFIX))

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

        if create_chromosome_img:
            create_chromosome_img(positions=chr_cpg_dict[chr_name], reads=table, chr_number=chr_name,
                                  output_folder=output_folder)

        df = pd.DataFrame(data=table.astype(np.int), columns=chr_cpg_dict[chr_name].astype(np.int))
        df.to_csv(os.path.join(output_folder, "chr%s_full_mapping.csv" % chr_name))


def create_chromosome_img(positions, reads, chr_number, output_folder):
    """
    Create an image representing the read coverage across the chromosome, each point will show 50 CpG,
    will also write on the image the locations of the Centromere and several others locations across the
    genome
    """
    # Get average across samples and split to bins + take median for every bean for so every point will be
    # representing 50 CpG
    reads_avg = np.average(reads, 0) * 100
    number_of_points_per_bin = np.round(len(reads_avg) / 50)
    splitted_reads = np.array_split(reads_avg, number_of_points_per_bin)
    splitted_location = np.array_split(positions, number_of_points_per_bin)
    med_reads = [np.median(i) for i in splitted_reads]
    med_pos = [np.median(i) for i in splitted_location]

    # Get the location value for each 1000 values (look good at some graphs) and add the centromem
    # this is used as labels
    all_pos = med_pos[::1000]
    all_pos += CENETROMERE_DICT[chr_number]

    all_labels = [str(int(i)) for i in all_pos]
    all_labels[-2] = "Centromere"
    all_labels[-1] = "Centromere"

    # Plot the figure
    plt.figure(figsize=(50, 8))
    ax_scatter = plt.axes()
    ax_scatter.scatter(med_pos, med_reads, c=med_reads)
    ax_scatter.set_ylim(top=max(med_reads))
    ax_scatter.set_xticks(all_pos)
    ax_scatter.set_xticklabels(all_labels)

    plt.savefig(os.path.join(output_folder, "coverage_along_chr_%s.png" % chr_number))
    plt.close()


def extract_bedgraph_information(input_folder, output_folder):
    """
    Create a bedgraph file for each patient with the coverage across different samples
    """
    chr_cpg_dict = files_tools.get_all_cpg_locations_across_chr()
    all_files = glob.glob(os.path.join(input_folder, FOLDER_SUFFIX))
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
    Run on all the patient and create bedgraph file
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
    """
    chr_cpg_dict = files_tools.get_all_cpg_locations_across_chr()
    coverage_across_patient_dict = {}

    all_files = glob.glob(os.path.join(input_folder, FOLDER_SUFFIX))

    for file_path in tqdm(all_files):
        data, patient, chr_name = get_file_info(file_path)
        coverage_across_samples = data.count()

        # Init the coverage_across_patient_dict dict
        if patient not in coverage_across_patient_dict:
            coverage_across_patient_dict[patient] = collections.Counter()

        coverage_across_patient_dict[patient].update(coverage_across_samples.values / data.shape[0])

    # Covered of each patient, global on chr
    for patient in coverage_across_patient_dict:
        commons.data_tools.counter_to_csv(coverage_across_patient_dict[patient],
                                          os.path.join(output_folder, "%s_covered_samples_counter.csv" % patient))


def main():
    parser = format_args()
    args = parser.parse_args()
    input_folder, output_folder = args.files_folder, args.output_folder

    extract_bedgraph_information(input_folder, output_folder)
    #extract_coverage_per_patient(input_folder, output_folder)
    #extract_coverage_per_chr(input_folder, output_folder)
    #extract_coverage_global(input_folder, output_folder)


if __name__ == '__main__':
    commons.slurm_tools.init_slurm(main)
