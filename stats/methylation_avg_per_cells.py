import argparse
import glob
import os
import re
import sys
import warnings

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from tqdm import tqdm

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

# colors where A, B, C are each in similar shades
CONVERT_COLOURS = {'A0': '#223b7c', 'A1': '#26428b', 'A2': '#2a499a', 'A3': '#2e51aa', 'A4': '#3358b9', 'A5': '#375fc8',
                   'A6': '#476ccd', 'A7': '#5678d1', 'A8': '#6584d5', 'A9': '#7591d9', 'B': '#42d142', 'B0': '#52d552',
                   'B1': '#61d961', 'B2': '#71dc71', 'B3': '#81e081', 'C0': '#d251cd', 'C1': '#d660d2', 'C2': '#da70d6',
                   'C3': '#de80da', 'C4': '#e28fdf', 'C5': '#e69fe3'}
NAN_COLOUR = '#9a9a9a'

# colors used in paper
CONVERT_COLOURS = {'A0': '#F8B4C0', 'A1': '#E0E346', 'A2': '#6F51A1', 'A3': '#F89C31', 'A4': '#EF292A',
                   'A5': '#A45AA4',
                   'A6': '#993232', 'A7': '#2256A6', 'A8': '#BC84BA', 'A9': '#ED3095', 'B': '#3B86C6',
                   'B0': '#3B86C6',
                   'B1': '#66CAD4', 'B2': '#6D8841', 'B3': '#28898A', 'C0': '#E7CFA0', 'C1': '#DDBD7E',
                   'C2': '#D1A34A',
                   'C3': '#B89979', 'C4': '#AA845D', 'C5': '#8C6E4A', 'undefined': '#BBBABC'}


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


def get_sublineage(sublineage_info, patient, cell):
    try:
        return sublineage_info[patient + '_' + cell]
    except:
        return "undefined"


def get_colours(sublineage_info, patient, cell):
    try:
        return CONVERT_COLOURS[sublineage_info[patient + '_' + cell]]
    except:
        return CONVERT_COLOURS["undefined"]


def main():
    args = parse_input()

    methylation_folder = args.methylation_folder
    all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                global_windows_data)
    sublineage = files_tools.load_compressed_pickle(
        R"/cs/usr/liorf/PycharmProjects/proj_scwgbs/stats/convert_sublineage.pickle.zlib")
    # sublineage = files_tools.load_compressed_pickle(
    #     R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\stats\top_bottom\convert_sublineage.pickle.zlib")

    pa_cells = set([])
    for patient in tqdm(patients_dict):
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            pa_cells |= set(df.index)

        pa_df_meth = pd.DataFrame(index=list(pa_cells))
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            pa_df_meth.loc[df.index, chromosome] = np.mean(df, axis=1)

        values = pa_df_meth.median(axis=1).sort_values()
        colors = [get_colours(sublineage, patient, cell) for cell in values.index]
        x = [i for i in range(values.shape[0])]
        legend_names = sorted({get_sublineage(sublineage, patient, cell) for cell in values.index})
        custom_lines = [Line2D([0], [0], color=CONVERT_COLOURS[i], ls=' ', marker='s', ms=12) for i in
                        legend_names]
        plt.bar(x, values, color=colors)
        plt.xlabel("tumor cells (self sorted)")
        plt.ylabel("avg meth")
        lgd = plt.legend(custom_lines, legend_names, loc='center left', bbox_to_anchor=(1.0, 0.5),
                         title='Genetic\nsub-lineages:')
        plt.title("cell methylation means for %s" % patient)
        plt.savefig("cell_methylation_means_for_%s_paperc.png" % patient, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
