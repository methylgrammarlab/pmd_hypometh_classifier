# todo: lior
import argparse
import glob
import os
import re
import sys
import warnings

import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from tqdm import tqdm

import commons.data_tools
import commons.utils
from commons.utils import get_nc_avg

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from covariance import covariance_to_bedgraph
from format_files import format_cpg_context_map
import commons.consts as consts

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
PMD_DICT_FILE_FORMAT = "%s_pmd_dicts.pkl.zip"
SEQ_SIZE = 150
TOP_LOW_PERCENTAGE_TO_REMOVE = 5
SOLO = 0

# colors where A, B, C are each in similar shades
CONVERT_COLOURS = {'A0': '#223b7c', 'A1': '#26428b', 'A2': '#2a499a', 'A3': '#2e51aa', 'A4': '#3358b9', 'A5': '#375fc8',
                   'A6': '#476ccd', 'A7': '#5678d1', 'A8': '#6584d5', 'A9': '#7591d9', 'B': '#42d142', 'B0': '#52d552',
                   'B1': '#61d961', 'B2': '#71dc71', 'B3': '#81e081', 'C0': '#d251cd', 'C1': '#d660d2', 'C2': '#da70d6',
                   'C3': '#de80da', 'C4': '#e28fdf', 'C5': '#e69fe3'}
NAN_COLOUR = '#9a9a9a'

# colors used in paper
SUBLINEAGE_COLOURS = {'NC': '#6FBE44', 'A0': '#F8B4C0', 'A1': '#E0E346', 'A2': '#6F51A1', 'A3': '#F89C31',
                      'A4': '#EF292A', 'A5': '#A45AA4',
                      'A6': '#993232', 'A7': '#2256A6', 'A8': '#BC84BA', 'A9': '#ED3095', 'B': '#3B86C6',
                      'B0': '#3B86C6',
                      'B1': '#66CAD4', 'B2': '#6D8841', 'B3': '#28898A', 'C0': '#E7CFA0', 'C1': '#DDBD7E',
                      'C2': '#D1A34A',
                      'C3': '#B89979', 'C4': '#AA845D', 'C5': '#8C6E4A', 'undefined': '#BBBABC'}
REGION_COLOURS = {'NC': '#6FBE44', 'PT1': '#F57E32', 'PT2': '#3C8A45', 'PT3': '#e4373b', 'PT4': '#53ABDA',
                  'PT5': '#EBE94F', 'PT6': '#EFB0BC', 'LN1': '#B6B6BA', 'LN2': '#B28A8E', 'LN3': '#866c6f',
                  'LN4': '#67532B', 'LN5': '#514321',
                  'ML1': '#94D6E4', 'ML2': '#4872B7', 'ML3': '#1C49A0', 'ML4': '#333463', 'ML5': '#464B7D',
                  'ML6': '#3C3E69', 'MP1': '#CFB0D3',
                  'MP2': '#BC85BB', 'MP3': '#8254A2', 'MP4': '#842F8D', 'MP5': '#632E65', 'MO1': '#E7CF9F',
                  'MO2': '#E2C68E', 'MO3': '#E2C68E', 'MO4': '#D9B36B', 'MO5': '#D5AB5B', 'MO6': '#D1A349'}
LESION = {"NC": "#6FBE44", "PT": "#F380A3", "LN": "#834D40", "MO": "#D1A349", 'ML': '#1D489E', "MP": "#815BA6"}
STRONG_WEAK = {"strong": "SCGS", "weak": "WCGW"}
WACGAW = {"AACGAA": '$\\bf{A}$ACGA$\\bf{A}$', "AACGAT": '$\\bf{A}$ACGA$\\bf{T}$',
          "TACGAA": '$\\bf{T}$ACGA$\\bf{A}$', "TACGAT": '$\\bf{T}$ACGA$\\bf{T}$'}

# PATTERNS_TO_CALC = ["ACGAA", "ACGAT", "ACGTA", "ACGTT", "AACGA", "ATCGA", "TACGA", "TTCGA"]
PATTERNS_TO_CALC = ['AACGAA', 'AACGAT', 'AACGTA', 'AACGTT', 'ATCGAA', 'ATCGAT', 'ATCGTA', 'ATCGTT', 'TACGAA', 'TACGAT',
                    'TACGTA', 'TACGTT', 'TTCGAA', 'TTCGAT', 'TTCGTA', 'TTCGTT']  # WWCGWW
PATTERNS_TO_CALC = ['ACGA', 'ACGC', 'ACGG', 'ACGT',
                    'CCGA', 'CCGC', 'CCGG', 'CCGT',
                    'GCGA', 'GCGC', 'GCGG', 'GCGT',
                    'TCGA', 'TCGC', 'TCGG', 'TCGT']  # [SW]CG[SW]
PATTERNS_TO_CALC = ['ANCGNA', 'ANCGNC', 'ANCGNG', 'ANCGNT',
                    'CNCGNA', 'CNCGNC', 'CNCGNG', 'CNCGNT',
                    'GNCGNA', 'GNCGNC', 'GNCGNG', 'GNCGNT',
                    'TNCGNA', 'TNCGNC', 'TNCGNG', 'TNCGNT']

CONVERT_PATTERN = {
    'ANCGNA': ['AACGAA', 'AACGCA', 'AACGGA', 'AACGTA', 'ACCGAA', 'ACCGCA', 'ACCGGA', 'ACCGTA', 'AGCGAA', 'AGCGCA',
               'AGCGGA', 'AGCGTA', 'ATCGAA', 'ATCGCA', 'ATCGGA', 'ATCGTA'],
    'ANCGNC': ['AACGAC', 'AACGCC', 'AACGGC', 'AACGTC', 'ACCGAC', 'ACCGCC', 'ACCGGC', 'ACCGTC', 'AGCGAC', 'AGCGCC',
               'AGCGGC', 'AGCGTC', 'ATCGAC', 'ATCGCC', 'ATCGGC', 'ATCGTC'],
    'ANCGNG': ['AACGAG', 'AACGCG', 'AACGGG', 'AACGTG', 'ACCGAG', 'ACCGCG', 'ACCGGG', 'ACCGTG', 'AGCGAG', 'AGCGCG',
               'AGCGGG', 'AGCGTG', 'ATCGAG', 'ATCGCG', 'ATCGGG', 'ATCGTG'],
    'ANCGNT': ['AACGAT', 'AACGCT', 'AACGGT', 'AACGTT', 'ACCGAT', 'ACCGCT', 'ACCGGT', 'ACCGTT', 'AGCGAT', 'AGCGCT',
               'AGCGGT', 'AGCGTT', 'ATCGAT', 'ATCGCT', 'ATCGGT', 'ATCGTT'],
    'CNCGNA': ['CACGAA', 'CACGCA', 'CACGGA', 'CACGTA', 'CCCGAA', 'CCCGCA', 'CCCGGA', 'CCCGTA', 'CGCGAA', 'CGCGCA',
               'CGCGGA', 'CGCGTA', 'CTCGAA', 'CTCGCA', 'CTCGGA', 'CTCGTA'],
    'CNCGNC': ['CACGAC', 'CACGCC', 'CACGGC', 'CACGTC', 'CCCGAC', 'CCCGCC', 'CCCGGC', 'CCCGTC', 'CGCGAC', 'CGCGCC',
               'CGCGGC', 'CGCGTC', 'CTCGAC', 'CTCGCC', 'CTCGGC', 'CTCGTC'],
    'CNCGNG': ['CACGAG', 'CACGCG', 'CACGGG', 'CACGTG', 'CCCGAG', 'CCCGCG', 'CCCGGG', 'CCCGTG', 'CGCGAG', 'CGCGCG',
               'CGCGGG', 'CGCGTG', 'CTCGAG', 'CTCGCG', 'CTCGGG', 'CTCGTG'],
    'CNCGNT': ['CACGAT', 'CACGCT', 'CACGGT', 'CACGTT', 'CCCGAT', 'CCCGCT', 'CCCGGT', 'CCCGTT', 'CGCGAT', 'CGCGCT',
               'CGCGGT', 'CGCGTT', 'CTCGAT', 'CTCGCT', 'CTCGGT', 'CTCGTT'],
    'GNCGNA': ['GACGAA', 'GACGCA', 'GACGGA', 'GACGTA', 'GCCGAA', 'GCCGCA', 'GCCGGA', 'GCCGTA', 'GGCGAA', 'GGCGCA',
               'GGCGGA', 'GGCGTA', 'GTCGAA', 'GTCGCA', 'GTCGGA', 'GTCGTA'],
    'GNCGNC': ['GACGAC', 'GACGCC', 'GACGGC', 'GACGTC', 'GCCGAC', 'GCCGCC', 'GCCGGC', 'GCCGTC', 'GGCGAC', 'GGCGCC',
               'GGCGGC', 'GGCGTC', 'GTCGAC', 'GTCGCC', 'GTCGGC', 'GTCGTC'],
    'GNCGNG': ['GACGAG', 'GACGCG', 'GACGGG', 'GACGTG', 'GCCGAG', 'GCCGCG', 'GCCGGG', 'GCCGTG', 'GGCGAG', 'GGCGCG',
               'GGCGGG', 'GGCGTG', 'GTCGAG', 'GTCGCG', 'GTCGGG', 'GTCGTG'],
    'GNCGNT': ['GACGAT', 'GACGCT', 'GACGGT', 'GACGTT', 'GCCGAT', 'GCCGCT', 'GCCGGT', 'GCCGTT', 'GGCGAT', 'GGCGCT',
               'GGCGGT', 'GGCGTT', 'GTCGAT', 'GTCGCT', 'GTCGGT', 'GTCGTT'],
    'TNCGNA': ['TACGAA', 'TACGCA', 'TACGGA', 'TACGTA', 'TCCGAA', 'TCCGCA', 'TCCGGA', 'TCCGTA', 'TGCGAA', 'TGCGCA',
               'TGCGGA', 'TGCGTA', 'TTCGAA', 'TTCGCA', 'TTCGGA', 'TTCGTA'],
    'TNCGNC': ['TACGAC', 'TACGCC', 'TACGGC', 'TACGTC', 'TCCGAC', 'TCCGCC', 'TCCGGC', 'TCCGTC', 'TGCGAC', 'TGCGCC',
               'TGCGGC', 'TGCGTC', 'TTCGAC', 'TTCGCC', 'TTCGGC', 'TTCGTC'],
    'TNCGNG': ['TACGAG', 'TACGCG', 'TACGGG', 'TACGTG', 'TCCGAG', 'TCCGCG', 'TCCGGG', 'TCCGTG', 'TGCGAG', 'TGCGCG',
               'TGCGGG', 'TGCGTG', 'TTCGAG', 'TTCGCG', 'TTCGGG', 'TTCGTG'],
    'TNCGNT': ['TACGAT', 'TACGCT', 'TACGGT', 'TACGTT', 'TCCGAT', 'TCCGCT', 'TCCGGT', 'TCCGTT', 'TGCGAT', 'TGCGCT',
               'TGCGGT', 'TGCGTT', 'TTCGAT', 'TTCGCT', 'TTCGGT', 'TTCGTT']}

FIRST_SUPERIMPOSED_BLUE = '#4872B7'
ALPHA_1_BLUE = '#ECF0F7'
ALPHA_3_BLUE = '#C7D4E9'
SECOND_SUPERIMPOSED_PURPLE = '#8254A2'
ALPHA_1_PURPLE = '#F2EDF5'
ALPHA_3_PURPLE = '#D9CBE2'
THIRD_SUPERIMPOSED_GREEN = '#3C8A45'
ALPHA_1_GREEN = '#EBF3EC'
ALPHA_3_GREEN = '#C4DBC6'
FOURTH_SUPERIMPOSED_PINK = '#f376b9'
ALPHA_1_PINK = '#FDF1F7'
ALPHA_3_PINK = '#FBD5E9'

LESION_ABBREVIATION = {'ML': 'Liver Metastasis', "PT": "Primary Tumor", "LN": "Lymph Node\nMetastasis",
                       "MP": "Post-treatment\nLiver Metastasis", "NC": "Normal Cell", "MO": "Omental Metastasis"}


# PATTERNS_TO_CALC = ['CACGAC', 'CACGAG', 'CACGTC', 'CACGTG', 'CTCGAC', 'CTCGAG', 'CTCGTC', 'CTCGTG', 'GACGAC', 'GACGAG',
#                     'GACGTC', 'GACGTG', 'GTCGAC', 'GTCGAG', 'GTCGTC', 'GTCGTG']  # SWCGWS
# PATTERNS_TO_CALC = ['CACGAC', 'CACGAG', 'CACGAA', 'CACGAT', 'GACGAC', 'GACGAG', 'GACGAA', 'GACGAT', 'AACGAC', 'AACGAG',
#                     'AACGAA', 'AACGAT', 'TACGAC', 'TACGAG', 'TACGAA', 'TACGAT']  # [SW]ACGA[SW]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--patients_data', help='Path to files with the calculated data', required=False)
    parser.add_argument('--pmd_dicts', help='Path to location with the pre-calculated PMD dicts', required=False)
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


# def get_patient_dict(all_file_paths):
#     d = {}
#     for file_path in all_file_paths:
#         try:
#             patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
#         except:
#             continue
#
#         if patient not in d:
#             d[patient] = []
#
#         d[patient].append((chromosome, file_path))
#
#     return d


def get_cancer_methylation_of_patient(patient, chromosome, indexes, methylation_folder):
    methylation_file_path = os.path.join(methylation_folder, patient,
                                         METHYLATION_FILE_FORMAT % (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    df = commons.utils.filter_df_based_on_region_name(df, region_name=commons.utils.ONLY_PT)
    mdf = df.mean()
    return mdf.loc[indexes].values


def get_sublineage(cell, sublineage_info, patient):
    try:
        return sublineage_info[patient + '_' + cell]
    except:
        return "undefined"


def get_region(cell):
    return cell.split('_')[0]


def get_colours(key, col_name):
    if col_name == 'sublineage':
        return SUBLINEAGE_COLOURS[key]
    elif col_name == 'region':
        return REGION_COLOURS[key]
    elif col_name == 'lesion':
        return LESION[key]


def get_legend_title(col_name):
    if col_name == 'sublineage':
        return 'Genetic\nsub-lineages:'
    elif col_name == 'region':
        return 'Sampling regions:'


def get_label_name(label):
    try:
        return STRONG_WEAK[label]
    except:
        try:
            return WACGAW[label]
        except:
            return label


def create_multiple_plots(df, color_by_superimpose_first, to_plot_superimpose_second, superimpose=False, third=None,
                          fourth=None, plot_all=None, ylen=14, boxplot=False):
    fig, axs = plt.subplots(3, 4, figsize=(22, ylen))
    layout_dict = {'CRC01': (0, 0), 'CRC02': (0, 1), 'CRC04': (0, 2), 'CRC09': (0, 3), 'CRC10': (1, 0), 'CRC11': (1, 1),
                   'CRC12': (1, 2), 'CRC13': (1, 3), 'CRC14': (2, 1), 'CRC15': (2, 2)}
    for patient in tqdm(df.patient.unique()):
        if superimpose:
            create_superimposed_graphs(df, patient, axs[layout_dict[patient]], color_by_superimpose_first,
                                       to_plot_superimpose_second, third, fourth)
        elif plot_all:
            plot_16(df, patient, axs[layout_dict[patient]])
        elif boxplot:
            plot_boxplot(df, patient, axs[layout_dict[patient]], color_by_superimpose_first)
        else:
            create_colored_graphs(df, color_by_superimpose_first, patient, axs[layout_dict[patient]],
                                  to_plot_superimpose_second)

    for ax in axs.flat:
        ax.set(xlabel="tumor cells (self sorted)", ylabel="median of average methylation")
        if boxplot:
            ax.set(xlabel="Lesion", ylabel="median of average methylation")
    for ax in axs.flat:
        ax.label_outer()
    axs[2, 0].set_visible(False)
    axs[2, 3].set_visible(False)

    if superimpose:
        fig_title = 'cell methylation %s vs %s' % (
            get_label_name(color_by_superimpose_first), get_label_name(to_plot_superimpose_second))
        path = "solo_nc/cell_methylation_means_superimposed_%s_%s_solo_nc" % (
            color_by_superimpose_first, to_plot_superimpose_second)
        labels = [get_label_name(color_by_superimpose_first), get_label_name(color_by_superimpose_first) + "-NC",
                  get_label_name(to_plot_superimpose_second), get_label_name(to_plot_superimpose_second) + "-NC"]
        custom_lines = [
            Line2D([0], [0], color=ALPHA_1_BLUE, ls=' ', marker='s', ms=14, mec=FIRST_SUPERIMPOSED_BLUE, mew=1.0),
            Line2D([0], [0], color=ALPHA_3_BLUE, ls=' ', marker='s', ms=14, mec=FIRST_SUPERIMPOSED_BLUE, mew=1.0),
            Line2D([0], [0], color=ALPHA_1_PURPLE, ls=' ', marker='s', ms=14, mec=SECOND_SUPERIMPOSED_PURPLE, mew=1.0),
            Line2D([0], [0], color=ALPHA_3_PURPLE, ls=' ', marker='s', ms=14, mec=SECOND_SUPERIMPOSED_PURPLE, mew=1.0)]
        if third:
            fig_title += " vs %s" % third
            path += "_%s" % third
            labels.append(get_label_name(third))
            labels.append(get_label_name(third) + '-NC')
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_1_GREEN, ls=' ', marker='s', ms=14, mec=THIRD_SUPERIMPOSED_GREEN, mew=1.0))
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_3_GREEN, ls=' ', marker='s', ms=14, mec=THIRD_SUPERIMPOSED_GREEN, mew=1.0))
        if fourth:
            fig_title += " vs %s" % fourth
            path += "_%s" % fourth
            labels.append(get_label_name(fourth))
            labels.append(get_label_name(fourth) + '-NC')
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_1_PINK, ls=' ', marker='s', ms=14, mec=FOURTH_SUPERIMPOSED_PINK, mew=1.0))
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_3_PINK, ls=' ', marker='s', ms=14, mec=FOURTH_SUPERIMPOSED_PINK, mew=1.0))
        fig.suptitle(fig_title, fontsize=24)
        lgd = fig.legend(custom_lines, labels, loc='center right', fontsize=14)
        plt.savefig(path)
    elif plot_all:
        fig.suptitle('Cell methylation wwcgww smooth', fontsize=24)
        handles, labels = axs[0, 0].get_legend_handles_labels()
        order = list(df.loc[:, PATTERNS_TO_CALC].mean(axis=0).sort_values(ascending=False).index)
        lgd = fig.legend([handles[labels.index(lab)] for lab in order], order, loc='center right',
                         fontsize=14)
        for line in lgd.get_lines():
            line.set_linewidth(2)
        plt.savefig("solo_nc/wwcgww_solo_nc.png")
    elif boxplot:
        fig.suptitle('Cell methylation', fontsize=24)
        legend_names = sorted(LESION.keys())
        custom_lines = [Line2D([0], [0], color=get_colours(key, 'lesion'), ls=' ', marker='s', ms=14)
                        for key in
                        legend_names]
        # legend_names = [LESION_ABBREVIATION[key] for key in legend_names]
        title = "Lesion:"
        lgd = fig.legend(custom_lines, legend_names, loc='center right', fontsize=14)
        lgd.set_title(title, prop={'size': 14})
        plt.savefig("solo_nc/boxplot_solo_nc.png")
    else:
        fig.suptitle('cell methylation for %s colored by %s' % (to_plot_superimpose_second, color_by_superimpose_first),
                     fontsize=24)
        legend_names = sorted(df.loc[:, color_by_superimpose_first].unique())
        custom_lines = [Line2D([0], [0], color=get_colours(key, color_by_superimpose_first), ls=' ', marker='s', ms=14)
                        for key in
                        legend_names]
        title = get_legend_title(color_by_superimpose_first)
        lgd = fig.legend(custom_lines, legend_names, loc='center right', fontsize=14)
        lgd.set_title(title, prop={'size': 14})
        plt.savefig("solo_nc/cell_methylation_means_colored_by_%s_for_%s_solo_nc.png" % (
            color_by_superimpose_first, to_plot_superimpose_second))
    plt.show()


def plot_16(df, patient, axs):
    df = df.sort_values(by=['median'])
    unique_colours = matplotlib.colors.ListedColormap(["#000000", "#FF0000", "#A0A0A0", "#FF7F7F",
                                                       "#2A4BD7", "#1D6914", "#814A19", "#8126C0",
                                                       "#9DAFFF", "#81C57A", "#E9DEBB", "#AD2323",
                                                       "#29D0D0", "#FFEE33", "#FF9233", "#FF30D5"])

    # df2 = df.loc[df.patient == patient, PATTERNS_TO_CALC].reset_index().interpolate(method='cubic')
    # df2.plot(ax=axs, legend=False, linewidth=1, xticks=[], colormap=unique_colours)

    df.loc[df.patient == patient, PATTERNS_TO_CALC].plot(ax=axs, legend=False, linewidth=1, xticks=[],
                                                         colormap=unique_colours)
    axs.set_title("%s" % patient)
    axs.set_ylim([0, 1])


def create_superimposed_graphs(df, patient, axs, first, second, third=None, fourth=None):
    single = not axs
    if single:
        figs, axs = plt.subplots()
    df = df.sort_values(by=['median'])
    add_superimposition_to_graph(axs, df, first, patient, FIRST_SUPERIMPOSED_BLUE)
    add_superimposition_to_graph(axs, df, second, patient, SECOND_SUPERIMPOSED_PURPLE)
    if third:
        add_superimposition_to_graph(axs, df, third, patient, THIRD_SUPERIMPOSED_GREEN)
    if fourth:
        add_superimposition_to_graph(axs, df, fourth, patient, FOURTH_SUPERIMPOSED_PINK)

    axs.set_title("%s" % patient)
    axs.set_ylim([0, 1])
    if single:
        # axs.set_title("Cell Methylation %s vs %s, %s" % (get_label_name(first), get_label_name(second), patient))
        axs.set_title("Cell Methylation $\\bf{W}$ACGA$\\bf{W}$, %s" % patient)
        labels = [get_label_name(first), get_label_name(first) + "-NC", get_label_name(second),
                  get_label_name(second) + "-NC"]
        custom_lines = [
            Line2D([0], [0], color=ALPHA_1_BLUE, ls=' ', marker='s', ms=8, mec=FIRST_SUPERIMPOSED_BLUE, mew=1.0),
            Line2D([0], [0], color=ALPHA_3_BLUE, ls=' ', marker='s', ms=8, mec=FIRST_SUPERIMPOSED_BLUE, mew=1.0),
            Line2D([0], [0], color=ALPHA_1_PURPLE, ls=' ', marker='s', ms=8, mec=SECOND_SUPERIMPOSED_PURPLE, mew=1.0),
            Line2D([0], [0], color=ALPHA_3_PURPLE, ls=' ', marker='s', ms=8, mec=SECOND_SUPERIMPOSED_PURPLE, mew=1.0)]
        if third:
            # fig_title += " vs %s" % third
            # path += "_%s" % third
            labels.append(get_label_name(third))
            labels.append(get_label_name(third) + '-NC')
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_1_GREEN, ls=' ', marker='s', ms=8, mec=THIRD_SUPERIMPOSED_GREEN, mew=1.0))
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_3_GREEN, ls=' ', marker='s', ms=8, mec=THIRD_SUPERIMPOSED_GREEN, mew=1.0))
        if fourth:
            # fig_title += " vs %s" % fourth
            # path += "_%s" % fourth
            labels.append(get_label_name(fourth))
            labels.append(get_label_name(fourth) + '-NC')
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_1_PINK, ls=' ', marker='s', ms=8, mec=FOURTH_SUPERIMPOSED_PINK, mew=1.0))
            custom_lines.append(
                Line2D([0], [0], color=ALPHA_3_PINK, ls=' ', marker='s', ms=8, mec=FOURTH_SUPERIMPOSED_PINK, mew=1.0))
        lgd = figs.legend(custom_lines, labels, loc='center right', bbox_to_anchor=(1.15, 0.5))
        axs.set(xlabel="tumor cells (self sorted)", ylabel="median of average methylation")
        # path = "solo_nc/cell_methylation_means_superimposed_%s_%s_solo_nc_patient_%s" % (first, second, patient)
        path = "solo_nc/cell_methylation_means_superimposed_WACGAW_solo_nc_patient_%s" % patient
        plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight', pad_inches=0.2)
        plt.show()


def add_superimposition_to_graph(axs, df, col_name, patient, color):
    values = df.loc[df.patient == patient, col_name]
    x = [i for i in range(values.shape[0])]
    colors = [color if key == "NC" else "#ffffff" for key in df.loc[df.patient == patient, "region"]]
    axs.plot(x, values, color=color, linewidth=1, drawstyle='steps-mid')
    axs.bar(x, values, width=1, alpha=0.3, color=colors)
    axs.bar(x, values, width=1, alpha=0.1, color=color)


def create_colored_graphs(df, color_by, patient, axs, to_plot):
    df = df.sort_values(by=['median'])
    values = df.loc[df.patient == patient, to_plot]
    colors = [get_colours(key, color_by) for key in df.loc[df.patient == patient, color_by]]
    x = [i for i in range(values.shape[0])]
    axs.bar(x, values, color=colors, width=1)
    axs.set_title("%s" % patient)
    axs.set_ylim([0, 1])

    ### For single graphs ###

    # plt.bar(x, values, color=colors, width=1)
    # plt.xlabel("tumor cells (self sorted)")
    # plt.ylabel("avg meth")
    # title = get_legend_title(col_name)
    # lgd = plt.legend(custom_lines, legend_names, loc='center left', bbox_to_anchor=(1.0, 0.5),
    #                  title=title)
    # plt.title("cell methylation means for %s" % patient)
    # plt.savefig("cell_methylation_means_for_%s_colored_%s.png" % (patient, col_name), bbox_extra_artists=(lgd,),
    #             bbox_inches='tight')
    # plt.close()


def plot_boxplot(df, patient, axs, split_by, type='box'):
    single = not axs
    if single:
        figs, axs = plt.subplots()
    data = {}
    for i in df.loc[df.patient == patient, split_by].apply(lambda x: x[:2]).unique():
        data[i] = df.loc[(df.patient == patient) & (df.loc[:, split_by].apply(lambda x: x[:2]) == i), "median"].dropna()
    df.loc[:, "lesion"] = df.loc[:, split_by].apply(lambda x: x[:2])

    data = {key: data[key] for key in LESION.keys() if key in data and data[key].size > 0}

    if type == 'box':
        bplot = axs.boxplot(data.values(), patch_artist=True)
        #
        colors = [get_colours(key, 'lesion') for key in data.keys()]
        for median, color in zip(bplot['medians'], colors):
            # median.set(color=color, linewidth=2)
            median.set(color='black', linewidth=2)

        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
    import seaborn as sns

    order = [key for key in LESION.keys() if key in data]
    if type == 'swarm':
        ax = sns.swarmplot(x="lesion", y="median", data=df.loc[df.patient == patient, :].dropna(), palette=LESION,
                           alpha=0.7, order=order)
    if type == 'strip':
        ax = sns.stripplot(x="lesion", y="median", data=df.loc[df.patient == patient, :].dropna(), palette=LESION,
                           alpha=0.7, order=order)
    if type == 'violin':
        ax = sns.violinplot(x="lesion", y="median", data=df.loc[df.patient == patient, :].dropna(), palette=LESION,
                            order=order)
    if type == 'boxpoint':
        ax = sns.swarmplot(x="lesion", y="median", data=df.loc[df.patient == patient, :].dropna(), palette=LESION,
                           alpha=0.7, order=order)
        ax = sns.boxplot(x="lesion", y="median", data=df.loc[df.patient == patient, :].dropna(), palette=LESION,
                         order=order,
                         showcaps=False, boxprops={'facecolor': 'None'},
                         showfliers=False, whiskerprops={'linewidth': 0})

    if type == 'swarm' or type == 'strip':
        # distance across the "X" or "Y" stipplot column to span, in this case 40%
        median_width = 0.4

        for tick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
            sample_name = text.get_text()  # lesion name
            median_val = df.loc[(df.patient == patient) & (df['lesion'] == sample_name), "median"].median()
            # plot horizontal lines across the column, centered on the tick
            ax.plot([tick - median_width / 2, tick + median_width / 2], [median_val, median_val],
                    lw=4, color='k')

    plt.title("%s" % patient)
    if single:
        axs.set_title("%s Lesions" % patient)
        if type == 'box':
            axs.set_xticklabels(data.keys())
        # legend_names = sorted(data.keys())
        legend_names = data.keys()
        custom_lines = [Line2D([0], [0], color=get_colours(key, 'lesion'), ls=' ', marker='s', ms=8)
                        for key in
                        legend_names]
        legend_names = [LESION_ABBREVIATION[key] for key in legend_names]
        title = "Lesion:"
        lgd = figs.legend(custom_lines, legend_names, loc='center right', bbox_to_anchor=(1.18, 0.5))
        lgd.set_title(title)
        axs.set(xlabel="Lesion", ylabel="median of average methylation")
        plt.savefig("solo_nc/%s_%splot_regions.png" % (patient, type), bbox_extra_artists=(lgd,),
                    bbox_inches='tight', pad_inches=0.2)
        plt.show()


def plot_histogram(df, patient, axs, split_by):
    single = not axs
    if single:
        figs, axs = plt.subplots()
    data = {}
    for i in df.loc[df.patient == patient, split_by].apply(lambda x: x[:2]).unique():
        data[i] = df.loc[(df.patient == patient) & (df.loc[:, split_by].apply(lambda x: x[:2]) == i), "median"].dropna()

    # for i in df.loc[df.patient == patient, split_by].unique():
    #     data[i] = df.loc[(df.patient == patient) & (df.loc[:, split_by] == i), "median"].fillna(0)
    import scipy.stats as st
    colors = [get_colours(key, 'lesion') for key in data.keys()]
    plt.hist(data.values(), label=data.keys(), alpha=0.7, density=True, color=colors)
    mn, mx = plt.xlim()
    plt.show()
    for (i, color) in zip(data, colors):
        plt.xlim(mn, mx)
        kde_xs = np.linspace(mn, mx, 301)
        kde = st.gaussian_kde(data[i])
        plt.plot(kde_xs, kde.pdf(kde_xs), label=i, color=color)
    plt.legend(loc='upper left')
    plt.title("CRC01 Regions")
    plt.savefig("solo_nc/CRC01_hist_regions.png")
    plt.show()

    import seaborn as sns
    for (i, color) in zip(data, colors):
        sns.distplot(data[i], hist=False, kde=True,
                     kde_kws={'shade': False, 'linewidth': 1},
                     label=i, color=color)
    plt.title("CRC01 seaborn")
    plt.show()

    if single:
        plt.show()


def create_graphs(path):
    df = pd.read_pickle(path)
    df.loc[:, "lesion"] = df.loc[:, 'region'].apply(lambda x: x[:2])
    cols = ['patient', 'region', 'lesion', 'sublineage', 'median', 'strong', 'weak', 'AACGAA', 'AACGAT', 'AACGTA',
            'AACGTT', 'ATCGAA', 'ATCGAT', 'ATCGTA', 'ATCGTT', 'TACGAA', 'TACGAT', 'TACGTA', 'TACGTT', 'TTCGAA',
            'TTCGAT', 'TTCGTA', 'TTCGTT']
    df = df[cols]

    # for patient in ["CRC01", "CRC02", "CRC04", "CRC09", "CRC10", "CRC11", "CRC12", "CRC13", "CRC14", "CRC14"]:
    #     weak_d = np.max(df.loc[df.patient == patient, "weak"]) - np.min(df.loc[df.patient == patient, "weak"])
    #     strong_d = np.max(df.loc[df.patient == patient, "strong"]) - np.min(df.loc[df.patient == patient, "strong"])
    #     print("delta for", patient, "weak:", weak_d, "strong:", strong_d)
    # print("weak_max", np.max(df.loc[df.patient == patient, "weak"]), 'strong_max', np.max(df.loc[df.patient == patient, "strong"]))
    df.to_csv("solo_nc/all_wwcgww_solo_nc.csv", index_label='cell_name')
    # create_multiple_plots(df, 'sublineage', 'median')
    # create_multiple_plots(df, 'region', 'median')
    # create_multiple_plots(df, 'sublineage', 'weak')
    # create_multiple_plots(df, 'region', 'weak')
    # create_multiple_plots(df, 'sublineage', 'strong')
    # create_multiple_plots(df, 'region', 'strong')
    # create_multiple_plots(df, 'ACGAA', 'ACGAT', superimpose=True, third='ACGTA', forth='ACGTT')
    # create_multiple_plots(df, 'AACGA', 'ATCGA', superimpose=True, third='TACGA', forth='TTCGA')
    # create_multiple_plots(df, 'sublineage', 'acgtt')
    # create_multiple_plots(df, 'region', 'acgtt')
    # create_multiple_plots(df, 'sublineage', 'acgtv')
    # create_multiple_plots(df, 'region', 'acgtv')
    # create_multiple_plots(df, 'acgtt', 'acgtv', superimpose=True)
    # create_multiple_plots(df, 'acgat', 'acgav', superimpose=True)
    # create_multiple_pdf.loc[(df.patient == patient) & (df.loc[:, split_by] == i), "median"]lots(df, 'acgaa', 'acgab', superimpose=True)
    # create_multiple_plots(df, 'strong', 'weak', superimpose=True)
    # create_multiple_plots(df, "", "", plot_all=True, ylen=30)
    # create_multiple_plots(df, "region", "", boxplot=True)
    # plot_boxplot(df, "CRC13", None, "region", 'box')
    # plot_boxplot(df, "CRC11", None, "region", 'swarm')
    # plot_boxplot(df, "CRC13", None, "region", 'strip')
    # plot_boxplot(df, "CRC13", None, "region", 'violin')
    # plot_boxplot(df, "CRC13", None, "region", 'boxpoint')
    # plot_boxplot(df, "CRC01", None, "region")
    # plot_histogram(df, "CRC01", None, "region")
    # create_multiple_plots(df, 'AACGAA', 'AACGAT', superimpose=True, third='TACGAA', forth='TACGAT')
    # create_superimposed_graphs(df, "CRC01", None, 'strong', 'weak')
    # create_superimposed_graphs(df, "CRC11", None, 'TACGAT', 'AACGAA', 'TACGAA', 'AACGAT')

    # for patient in ['CRC01', 'CRC10', 'CRC11', 'CRC13']:
    # for patient in df.patient.unique():
    #     create_colored_graphs(df, 'sublineage', patient)
    #     create_colored_graphs(df, 'region', patient)


def translate_pattern(pattern):
    if 'N' in pattern:
        return CONVERT_PATTERN[pattern]
    else:
        return [pattern]


def get_all_pattern_cells_from_df(df, context_info, pattern):
    pattern_cells = set([])
    all_patterns = translate_pattern(pattern)
    for p in all_patterns:
        pattern_cells |= set(context_info[context_info.str.contains(p)].index)
    return np.mean(df.loc[:, df.columns & pattern_cells], axis=1)


def get_patient_df_dict(all_file_paths):
    patient_df_dict = {}
    for file_path in tqdm(all_file_paths):
        patient, chromosome = consts.DATA_FILE_SCWGBS_RE.findall(file_path)[0]
        patient_df_dict[patient] = {}

        df = pd.read_pickle(file_path)  # this gets really slow, todo dror - ideas?

        patient_df_dict[patient][chromosome] = df

    return patient_df_dict


def filter_patient_df_dict(patient_df_dict, boundaries_data, perc_of_cpg_to_remove_based_on_coverage=5):
    """
    For each of the 22 chromosomes for every patient:
    - filter out non PMD CpGs
    - filter based on boundaries
    - filter out low and high coverage CpGs
    :param patient_df_dict: Dictionary that for each patient holds a didtionary of the DFs of all it's chromosomes
    :param boundaries_data: The boundaries data
    :param perc_of_cpg_to_remove_based_on_coverage: Percentage of CpGs to remove from top and low coverage
    :return: A new dictionary of the same format, with the filtered DFs
    """
    filtered_dict = {patient: {} for patient in patient_df_dict}
    for patient in patient_df_dict:
        for chromosome in patient_df_dict[patient]:
            df = patient_df_dict[patient][chromosome]
            # filter oyt non-PMDs
            filtered_df = handle_pmds.filtered_out_non_pmd(df, chromosome, add_pmd_index=False, pmd_file=consts.PMD_FILE_LOCAL_LIOR)
            # filter based on given boundaries
            filtered_df = commons.utils.filter_df_based_on_tuple_list(filtered_df, boundaries_data[chromosome[3:]])
            # filter out low and high coverage cpg
            filtered_df = commons.utils.remove_extreme_cpgs_by_coverage(
                filtered_df, top_low_level_to_remove=perc_of_cpg_to_remove_based_on_coverage)
            filtered_dict[patient][chromosome] = filtered_df
    return filtered_dict


def main():
    args = parse_input()

    # if args.patients_data:
    #     create_graphs(args.patients_data)
    #     return

    nc_files = consts.NC_LIOR

    # methylation_folder = args.methylation_folder
    # all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    # global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    # patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_coverage_filter(all_files_dict,
    #                                                                                      global_windows_data)

    all_files = files_tools.get_files_to_work(args.methylation_folder, pattern=os.path.join("CRC09", "*.pkl.zip"))
    boundaries_data = files_tools.load_compressed_pickle(args.windows_file)

    patients_dict = get_patient_df_dict(all_files)
    patients_dict = filter_patient_df_dict(patients_dict, boundaries_data)

    sublineage = files_tools.load_compressed_pickle(consts.CONVERT_SUBLINEAGE_LIOR)

    # # Collect global sequence information from the cpg dict
    cpg_dict = files_tools.get_cpg_context_map(load_with_path=consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI_LIOR)
    strong = {}
    weak = {}
    context_info = {}
    cpg75flank = {}  # todo: dror
    nc_meth = {}
    for chr in range(1, 23):
        chr_info = cpg_dict['chr%d' % chr]
        weak['%d' % chr] = chr_info[chr_info[:, -2] == 1][:, 0]
        strong['%d' % chr] = chr_info[chr_info[:, -1] == 1][:, 0]
        context_info['%d' % chr] = pd.Series(chr_info[:, -3], index=chr_info[:, 0]).apply(
            format_cpg_context_map.convert_context_int_to_str)
        cpg75flank["%d" % chr] = chr_info[chr_info[:, -7] == SOLO][:, 0]
        nc_meth_avg = get_nc_avg(chr, nc_files)
        nc_meth["%d" % chr] = nc_meth_avg[nc_meth_avg > 0.5].index
        pass

    pa_cells = set([])
    all_patients_mean = []
    all_patients_median = []
    pmd_filter = 0
    solo_filter = 0
    solo_nc_filter = 0
    for patient in tqdm(patients_dict, desc="get data"):
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            pmd_filter += np.sum(df.mean() != 0)
            solo_filter += np.sum(df.loc[:, df.columns & cpg75flank[chromosome]].mean() != 0)
            df = df.loc[:, df.columns & cpg75flank[chromosome] & nc_meth[chromosome]]
            solo_nc_filter += np.sum(df.mean() != 0)

            pa_cells |= set(df.index)
        print("patient:", patient, "pmd:", pmd_filter, "solo:", solo_filter, "nc_solo:", solo_nc_filter)

        pa_df_meth = pd.DataFrame(index=list(pa_cells))
        pa_df_meth_strong = pd.DataFrame(index=list(pa_cells))
        pa_df_meth_weak = pd.DataFrame(index=list(pa_cells))
        pattern_cells = {}
        for pattern in PATTERNS_TO_CALC:
            pattern_cells[pattern] = pd.DataFrame(index=list(pa_cells))
        for chromosome in tqdm(patients_dict[patient], desc="chromosome data calc %s" % patient):
            df = patients_dict[patient][chromosome]
            pa_df_meth.loc[df.index, chromosome] = np.mean(df, axis=1)
            pa_df_meth_strong.loc[df.index, chromosome] = np.mean(df.loc[:, df.columns & strong[chromosome]], axis=1)
            pa_df_meth_weak.loc[df.index, chromosome] = np.mean(df.loc[:, df.columns & weak[chromosome]], axis=1)
            for pattern in PATTERNS_TO_CALC:
                pattern_cells[pattern].loc[df.index, chromosome] = get_all_pattern_cells_from_df(df,
                                                                                                 context_info[
                                                                                                     chromosome],
                                                                                                 pattern)

        pa_df_meth.reset_index().melt(id_vars='index', var_name='chr', value_name='tot_avg')

        ### save median data ###
        patient_df_median = pd.DataFrame(index=pa_df_meth.index)
        patient_df_median.loc[:, 'median'] = pa_df_meth.median(axis=1)
        patient_df_median.loc[:, 'strong'] = pa_df_meth_strong.median(axis=1)
        patient_df_median.loc[:, 'weak'] = pa_df_meth_weak.median(axis=1)
        patient_df_median.loc[:, 'patient'] = patient
        patient_df_median.loc[:, 'region'] = pa_df_meth.index
        patient_df_median.loc[:, 'region'] = patient_df_median.loc[:, 'region'].apply(get_region)
        patient_df_median.loc[:, 'sublineage'] = pa_df_meth.index
        patient_df_median.loc[:, 'sublineage'] = patient_df_median.loc[:, 'sublineage'].apply(get_sublineage,
                                                                                              args=(
                                                                                                  sublineage, patient))

        ### save mean data ###
        patient_df_mean = pd.DataFrame(index=pa_df_meth.index)
        patient_df_mean.loc[:, 'mean'] = pa_df_meth.mean(axis=1)
        patient_df_mean.loc[:, 'strong'] = pa_df_meth_strong.mean(axis=1)
        patient_df_mean.loc[:, 'weak'] = pa_df_meth_weak.mean(axis=1)
        patient_df_mean.loc[:, 'patient'] = patient
        patient_df_mean.loc[:, 'region'] = pa_df_meth.index
        patient_df_mean.loc[:, 'region'] = patient_df_mean.loc[:, 'region'].apply(get_region)
        patient_df_mean.loc[:, 'sublineage'] = pa_df_meth.index
        patient_df_mean.loc[:, 'sublineage'] = patient_df_mean.loc[:, 'sublineage'].apply(get_sublineage,
                                                                                          args=(sublineage, patient))

        ### add the data for the pattern ###
        for pattern in PATTERNS_TO_CALC:
            patient_df_median.loc[:, pattern] = pattern_cells[pattern].median(axis=1)
            patient_df_mean.loc[:, pattern] = pattern_cells[pattern].mean(axis=1)

        all_patients_median.append(patient_df_median)
        all_patients_mean.append(patient_df_mean)

    ### save to file ###
    all_patients_median_df = pd.concat(all_patients_median)
    all_patients_mean_df = pd.concat(all_patients_mean)
    median_path = "avg_data_all_NNCGNN_solo_nc_median.pickle.zlib"
    mean_path = "avg_data_all_NNCGNN_solo_nc_mean.pickle.zlib"
    # median_path = "all_patients_avg_data_all_wwcgww_solo_nc_median.pickle.zlib"
    # mean_path = "all_patients_avg_data_all_wwcgww_solo_nc_mean.pickle.zlib"
    all_patients_median_df.to_pickle(median_path)
    all_patients_mean_df.to_pickle(mean_path)
    all_patients_median_df.to_csv("avg_data_all_NNCGNN_solo_nc_median.csv")
    # create_graphs(median_path)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
