import os
import pickle
import sys

sys.path.append(os.getcwd())
from commons import consts


def get_sublineage_info(file_path=consts.SUBLINEAGE_FILE):
    with open(file_path, "rb") as f:
        return pickle.load(f)


def format_file(sublineage_file):
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

    with open("patient_sublineage_dict.pickle", "wb") as f:
        pickle.dump(sublineage_dict, f)


def main():
    format_file(sys.argv[1])


if __name__ == '__main__':
    main()
