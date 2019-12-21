import argparse
import collections
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))

import numpy as np
import pandas as pd

FOLDER_SUFFIX = os.path.join("CRC09", "all_cpg_ratios_CRC*_chr*.dummy.pkl.zip")
CHR_RE = re.compile(".*chr(\d+).*")
PATIENT_RE = re.compile("CRC\d\d")
COLUMNS = ["read_number", "counter"]


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
        create_chromosome_img(positions=chr_cpg_dict[chr_name], reads=table, chr_number=chr_name)

        df = pd.DataFrame(data=table.astype(np.int), columns=chr_cpg_dict[chr_name].astype(np.int))
        df.to_csv(os.path.join(output_folder, "chr%s_full_mapping.csv" % chr_name))

    for patient in covered_dict_per_patient:
        counter_df = pd.DataFrame.from_dict(covered_dict_per_patient[patient], orient='index').reset_index()
        counter_df.columns = COLUMNS
        counter_df.to_csv(os.path.join(output_folder, "%s_covered_samples_counter.csv" % patient), "w")

    counter_df = pd.DataFrame.from_dict(covered_counter, orient='index').reset_index()
    counter_df.columns = COLUMNS
    counter_df.to_csv(os.path.join(output_folder, "covered_samples_counter.csv"), "w")


def create_chromosome_img(positions, reads, chr_number):
    reads_avg = np.average(reads, 0) * 100

    splitted_reads = np.array_split(reads_avg, 50000)
    splitted_location = np.array_split(positions, 50000)

    med_reads = [np.median(i) for i in splitted_reads]
    med_pos = [np.median(i) for i in splitted_location]

    plt.figure(figsize=(50, 8))
    ax_scatter = plt.axes()

    # the scatter plot:
    ax_scatter.scatter(med_pos, med_reads, c=med_reads)
    ax_scatter.set_ylim(top=max(med_reads))

    all_pos = med_pos[::2000]
    all_pos.append(10766630)
    all_pos.append(14633713)

    ax_scatter.set_xticks(all_pos)
    all_labels = [str(i) for i in all_pos]
    all_labels[-2] = "Sgap"
    all_labels[-1] = "Egap"
    ax_scatter.set_xticklabels(all_labels)

    plt.savefig("coverage_along_chr_%s.png" % chr_number)
    plt.show()


def main():
    create_stats()


if __name__ == '__main__':
    main()
