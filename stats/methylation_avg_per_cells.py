import argparse
import glob
import os
import re
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from covariance import covariance_to_bedgraph

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
SEQ_SIZE = 150
TOP_LOW_PERCENTAGE_TO_REMOVE = 5


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_bedgraph_files(files):
    if os.path.isdir(files):
        file_path = os.path.join(files, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [files]

    return all_file_paths


def get_patient_dict(all_file_paths):
    d = {}
    for file_path in all_file_paths:
        try:
            patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
        except:
            continue

        if patient not in d:
            d[patient] = []

        d[patient].append((chromosome, file_path))

    return d


def get_cancer_methylation_of_patient(patient, chromosome, indexes, methylation_folder):
    methylation_file_path = os.path.join(methylation_folder, patient,
                                         METHYLATION_FILE_FORMAT % (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    _, df = covariance_to_bedgraph.get_region_df(df, sublineage_cells=[],
                                                 sublineage_name=covariance_to_bedgraph.ONLY_PT)
    mdf = df.mean()
    return mdf.loc[indexes].values


def main():
    args = parse_input()

    methylation_folder = args.methylation_folder
    all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                global_windows_data)

    pa_cells = set([])
    for patient in patients_dict:
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            pa_cells |= df.index

        pa_df_meth = pd.DataFrame(index=list(pa_cells))
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            pa_df_meth.loc[df.index, chromosome] = np.mean(df, axis=1)

        values = pa_df_meth.median(axis=1).sort_values()
        x = [i for i in range(values.shape[0])]
        plt.bar(x, values)
        plt.xlabel("tumor cells (self sorted)")
        plt.ylabel("avg meth")
        plt.title("cell methylation means for %s" % patient)
        plt.savefig("cell_methylation_means_for_%s.png" % patient)
        plt.close()


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
