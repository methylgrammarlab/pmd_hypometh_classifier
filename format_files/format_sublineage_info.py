"""
Convert a sublineage files to a dictionary we know how to work with
The input file format should be:

Patient_ID	Sample	Sublineage
CRC01	CRC01_LN1_172	A0
"""

import os
import sys

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import consts
from commons import files_tools

OUTPUT_FILE = "patient_sublineage_dict.pickle"


def get_sublineage_info(file_path=consts.SUBLINEAGE_FILE):
    return files_tools.load_pickle(file_path)


def format_sublineage_file(sublineage_file):
    sublineage_dict = {}
    with open(sublineage_file, "rb") as sublineage:
        _ = sublineage.readline()
        for line in sublineage:
            patient, sample, lineage = line.split()
            patient = patient.decode("ascii")
            sample = sample.decode("ascii")
            lineage = lineage.decode("ascii")

            sample = "_".join(sample.split("_")[1:])

            if patient not in sublineage_dict:
                sublineage_dict[patient] = {}

            if lineage not in sublineage_dict[patient]:
                sublineage_dict[patient][lineage] = []

            sublineage_dict[patient][lineage].append(sample)

    files_tools.save_pickle(OUTPUT_FILE, sublineage_dict)


def main():
    format_sublineage_file(sys.argv[1])


if __name__ == '__main__':
    main()
