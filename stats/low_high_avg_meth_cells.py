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
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    methylation_folder = args.methylation_folder
    patients = ["CRC01", "CRC04", "CRC10", "CRC13", "CRC14"]
    path = os.path.join(methylation_folder, "%s", "all*.pkl.zip")
    all_files_dict = get_patient_dict(
        [path for sublist in [glob.glob(path % patient) for patient in patients] for path in sublist])
    # all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                global_windows_data)
    bins_patients_dict = {}
    for patient in patients_dict:
        cell_mean = pd.concat(patients_dict[patient], axis=1).mean(axis=1).sort_values()
        bottom = min(len(cell_mean) - 1, int((len(cell_mean) / 100) * 15))
        top = int((len(cell_mean) / 100) * 15)
        low_15 = list(cell_mean.iloc[:bottom].index)
        high_15 = list(cell_mean.iloc[-top - 1:-1].index)
        bins_patients_dict[patient] = {"low": low_15, "high": high_15}
    output_path = os.path.join(output, "15perc_low_high_avg_methylation_cells.pickle.zlib")
    files_tools.save_as_compressed_pickle(output_path, bins_patients_dict)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
    # for patient in ["CRC01", "CRC04", "CRC10", "CRC13", "CRC14"]:
    #     dict = files_tools.load_compressed_pickle(
    #         R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\stats\top_bottom\%s_15perc_low_high_avg_methylation_cells.pickle.zlib" % patient)
    #     high_list = []
    #     low_list = []
    #     for chr in dict:
    #         high_list += dict[chr]["high"]
    #         low_list += dict[chr]["low"]
    #     pd.Series(high_list).value_counts().plot('bar')
    #     plt.title("%s high" % patient)
    #     plt.savefig("top_bottom\%s_high" % patient)
    #     plt.close()
    #     pd.Series(low_list).value_counts().plot('bar')
    #     plt.title("%s low" % patient)
    #     plt.savefig("top_bottom\%s_low" % patient)
    #     plt.close()
