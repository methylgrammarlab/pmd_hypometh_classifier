import os
import pickle
import sys

sys.path.append(os.getcwd())
from commons import consts


def get_sublinage_info(file_path=consts.SUBLINAGE_FILE):
    with open(file_path, "rb") as f:
        return pickle.load(f)


def format_file(sublinage_file):
    sublinage_dict = {}
    with open(sublinage_file, "rb") as sublinage:
        _ = sublinage.readline()
        for line in sublinage:
            patient, sample, linage = line.split()

            if patient not in sublinage_dict:
                sublinage_dict[patient] = {}

            if linage not in sublinage_dict[patient]:
                sublinage_dict[patient][linage] = []

            sublinage_dict[patient][linage].append(sample)

    with open("patient_sublinage_dict.pickle", "wb") as f:
        pickle.dump(sublinage_dict, f)


def main():
    format_file(sys.argv[1])


if __name__ == '__main__':
    main()
