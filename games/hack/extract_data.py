import glob
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd

CHR_NAME = re.compile("(chr\d+)")
DATA = "commonPMD"
BED = r"H:\Study\university\Computational-Biology\Year 3\Semester A\CBIO\CBioHackathon\PMD_coordinates_hg19.bedgraph"
FILES_PATH = r"H:\Study\university\Computational-Biology\Year 3\Semester A\CBIO\CBioHackathon\data"
OUTPUT_PATH = r"H:\Study\university\Computational-Biology\Year 3\Semester A\CBIO\CBioHackathon\data_pmds"


def create_pmd_dict():
    all_data = pd.read_csv(BED, sep="\t")._values
    chr_dict = {}

    chromose = np.unique(all_data[:, 0])
    for chro in chromose:
        chr_lines = all_data[all_data[:, 0] == chro]
        pmd_lines = chr_lines[chr_lines[:, -1] == DATA]
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

    with open("pmd_dict.pickle", "wb") as pmd_dict_file:
        pickle.dump(chr_dict, pmd_dict_file)


def split_files():
    with open("pmd_dict.pickle", "rb") as pmd_dict_f:
        pmd_dict = pickle.load(pmd_dict_f)

    files = glob.glob(os.path.join(FILES_PATH, sys.argv[1], "*_all_data.csv.gzip"))
    for f in files:
        i = 0
        patient = os.path.basename(os.path.dirname(f))
        chr_name = CHR_NAME.findall(f)[0]
        if chr_name == "chr1":
            pass
        csv_data = pd.read_csv(f, compression="gzip")
        pmd_list = pmd_dict[chr_name]
        for pmd_tuple in pmd_list:
            start, end = pmd_tuple
            pmd_mask = (csv_data['locations'] > start) & (csv_data['locations'] < end)
            pmd = csv_data.loc[pmd_mask, :]
            output_data = os.path.join(OUTPUT_PATH, patient, "pmd_%s_%s.csv.gzip" % (chr_name, i))
            if not os.path.exists(os.path.dirname(output_data)):
                os.mkdir(os.path.dirname(output_data))

            if pmd.size == 0:
                continue

            pmd.to_csv(output_data, compression="gzip")
            i += 1


if __name__ == '__main__':
    split_files()
    create_pmd_dict()
