import glob
import os
import pickle
import re

import numpy as np
import pandas as pd

CHR_NAME = re.compile("(chr\d+)")
PMD_LABEL = "commonPMD"


def create_pmd_dict(bed_file, output_path):
    all_data = pd.read_csv(bed_file, sep="\t")._values
    chr_dict = {}

    chromose = np.unique(all_data[:, 0])
    for chro in chromose:
        chr_lines = all_data[all_data[:, 0] == chro]
        pmd_lines = chr_lines[chr_lines[:, -1] == PMD_LABEL]
        chr_dict[chro] = []

        global_start = 100000
        global_end = 199999
        start_new = False
        for line in pmd_lines[1:]:
            start = line[1]
            end = line[2]
            if start_new:
                start_new = False
                global_start = start
                global_end = end
                continue

            elif global_end + 1 == start:
                global_end = end
                continue
            else:
                chr_dict[chro].append((global_start, global_end))
                start_new = True

    with open(os.path.join(output_path, "pmd_dict.pickle"), "wb") as pmd_dict_file:
        pickle.dump(chr_dict, pmd_dict_file)


def split_files(input_path, output_path, dict_path, patients):
    with open(dict_path, "rb") as pmd_dict_f:
        pmd_dict = pickle.load(pmd_dict_f)

    for patient in patients:
        files = glob.glob(os.path.join(input_path, patients, "*_all_data.csv.gzip"))
        for f in files:
            chr_name = CHR_NAME.findall(f)[0]
            csv_data = pd.read_csv(f, compression="gzip")
            pmd_list = pmd_dict[chr_name]
            for pmd_tuple in pmd_list:
                start, end = pmd_tuple
                pmd_mask = (csv_data['locations'] > start) & (csv_data['locations'] < end)
                pmd = csv_data.loc[pmd_mask, :]
                output_data = os.path.join(output_path, patient, "%s_pmd_%s_%s.csv.gzip" % (chr_name,
                                                                                            start, end))
                if not os.path.exists(os.path.dirname(output_data)):
                    os.mkdir(os.path.dirname(output_data))

                pmd.to_csv(output_data, compression="gzip")
