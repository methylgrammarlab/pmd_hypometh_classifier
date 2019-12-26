import argparse
import collections
import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

import numpy as np
import pandas as pd

FOLDER_SUFFIX = os.path.join("*", "all_cpg_ratios_CRC*_chr*.dummy.pkl.zip")

CHR_RE = re.compile(".*chr(\d+).*")
PATIENT_RE = re.compile("CRC\d\d")
COLUMNS = ["read_number", "counter"]

CENETROMERE_DICT = {
    "1": [121311574, 129523877],
    "2": [91278726, 97911435],
    "3": [88724336, 93610603],
    "4": [48396285, 52374061],
    "5": [45111338, 50515299],
    "6": [57556887, 63334797],
    "7": [57041911, 63035444],
    "8": [43719124, 48091035],
    "9": [47315670, 51717126],
    "10": [38196156, 42596634],
    "11": [52073942, 56281937],
    "12": [33028390, 38938733],
    "13": [16004126, 20341692],
    "14": [16172139, 19796928],
    "15": [15579446, 20905751],
    "16": [34499088, 39310184],
    "17": [22776839, 25729391],
    "18": [15615450, 19164415],
    "19": [24342712, 28719791],
    "20": [25701316, 29630179],
    "21": [10688588, 14563981],
    "22": [12326422, 17790024],
}


def format_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_folder', help='Path of the files folder', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()

    return args.files_folder, args.output_folder


def create_stats():
    input_folder, output_folder = format_args()
    all_files = glob.glob(os.path.join(input_folder, FOLDER_SUFFIX))
    chr_cpg_dict = {}
    not_nans_dict = {}
    covered_dict_per_patient = {}

    covered_counter = collections.Counter()

    for file_path in tqdm(all_files):
        patient = PATIENT_RE.findall(os.path.basename(file_path))[0]
        chr_name = CHR_RE.findall(file_path)[0]
        data = pd.read_pickle(file_path)
        not_nans = data.count()

        if chr_name not in chr_cpg_dict:
            not_nans_dict[chr_name] = []
            chr_cpg_dict[chr_name] = not_nans.index

        if patient not in covered_dict_per_patient:
            covered_dict_per_patient[patient] = collections.Counter()

        not_nans_dict[chr_name].append(not_nans.values / data.shape[0])
        covered_counter.update(not_nans.values / data.shape[0])
        covered_dict_per_patient[patient].update(not_nans.values / data.shape[0])

    # Combine all the chr to one
    for chr_name in chr_cpg_dict:
        table = np.array(not_nans_dict[chr_name])
        create_chromosome_img(positions=chr_cpg_dict[chr_name], reads=table, chr_number=chr_name,
                              output_folder=output_folder)

        df = pd.DataFrame(data=table.astype(np.int), columns=chr_cpg_dict[chr_name].astype(np.int))
        df.to_csv(os.path.join(output_folder, "chr%s_full_mapping.csv" % chr_name))

    for patient in covered_dict_per_patient:
        counter_df = pd.DataFrame.from_dict(covered_dict_per_patient[patient], orient='index').reset_index()
        counter_df.columns = COLUMNS
        counter_df.to_csv(os.path.join(output_folder, "%s_covered_samples_counter.csv" % patient), "w")

    counter_df = pd.DataFrame.from_dict(covered_counter, orient='index').reset_index()
    counter_df.columns = COLUMNS
    counter_df.to_csv(os.path.join(output_folder, "covered_samples_counter.csv"), "w")


def create_chromosome_img(positions, reads, chr_number, output_folder):
    reads_avg = np.average(reads, 0) * 100

    number_of_points_per_bin = np.round(len(reads_avg) / 50)
    splitted_reads = np.array_split(reads_avg, number_of_points_per_bin)
    splitted_location = np.array_split(positions, number_of_points_per_bin)

    med_reads = [np.median(i) for i in splitted_reads]
    med_pos = [np.median(i) for i in splitted_location]

    plt.figure(figsize=(50, 8))
    ax_scatter = plt.axes()

    # the scatter plot:
    ax_scatter.scatter(med_pos, med_reads, c=med_reads)
    ax_scatter.set_ylim(top=max(med_reads))

    all_pos = med_pos[::1000]
    all_pos += CENETROMERE_DICT[chr_number]

    ax_scatter.set_xticks(all_pos)
    all_labels = [str(int(i)) for i in all_pos]
    all_labels[-2] = "Centromere"
    all_labels[-1] = "Centromere"
    ax_scatter.set_xticklabels(all_labels)

    plt.savefig(os.path.join(output_folder, "coverage_along_chr_%s.png" % chr_number))


def fix_counter(file_name):
    new_file = file_name.replace(".csv", "_new.csv")
    data = open(file_name, "r").readlines()

    c = collections.Counter()
    for line in data[1:]:
        _, number, counter = line.strip().split("w")
        c.update({round(float(number) * 100): int(counter)})

    counter_df = pd.DataFrame.from_dict(c, orient='index').reset_index()
    counter_df.columns = COLUMNS
    counter_df.to_csv(new_file, "w")


def main():
    create_stats()


if __name__ == '__main__':
    main()
