import argparse
import glob
import os
import re
import sys
import warnings

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from tqdm import tqdm

from variance.get_valid_cpg import get_nc_avg

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

# colors where A, B, C are each in similar shades
CONVERT_COLOURS = {'A0': '#223b7c', 'A1': '#26428b', 'A2': '#2a499a', 'A3': '#2e51aa', 'A4': '#3358b9',
                   'A5': '#375fc8',
                   'A6': '#476ccd', 'A7': '#5678d1', 'A8': '#6584d5', 'A9': '#7591d9', 'B': '#42d142',
                   'B0': '#52d552',
                   'B1': '#61d961', 'B2': '#71dc71', 'B3': '#81e081', 'C0': '#d251cd', 'C1': '#d660d2',
                   'C2': '#da70d6',
                   'C3': '#de80da', 'C4': '#e28fdf', 'C5': '#e69fe3'}
NAN_COLOUR = '#9a9a9a'

# colors used in paper
SUBLINEAGE_COLOURS = {'A0': '#F8B4C0', 'A1': '#E0E346', 'A2': '#6F51A1', 'A3': '#F89C31', 'A4': '#EF292A',
                      'A5': '#A45AA4',
                      'A6': '#993232', 'A7': '#2256A6', 'A8': '#BC84BA', 'A9': '#ED3095', 'B': '#3B86C6',
                      'B0': '#3B86C6',
                      'B1': '#66CAD4', 'B2': '#6D8841', 'B3': '#28898A', 'C0': '#E7CFA0', 'C1': '#DDBD7E',
                      'C2': '#D1A34A',
                      'C3': '#B89979', 'C4': '#AA845D', 'C5': '#8C6E4A', 'undefined': '#BBBABC'}
REGION_COLOURS = {'NC': '#6FBE44', 'PT1': '#F57E32', 'PT2': '#3C8A45', 'PT3': '#3C8A45', 'PT4': '#53ABDA',
                  'PT5': '#EBE94F', 'PT6': '#EFB0BC', 'LN1': '#B6B6BA', 'LN2': '#B28A8E', 'LN3': '#B28A8E',
                  'LN4': '#67532B', 'LN5': '#514321',
                  'ML1': '#94D6E4', 'ML2': '#4872B7', 'ML3': '#1C49A0', 'ML4': '#333463', 'ML5': '#464B7D',
                  'ML6': '#3C3E69', 'MP1': '#CFB0D3',
                  'MP2': '#BC85BB', 'MP3': '#8254A2', 'MP4': '#842F8D', 'MP5': '#632E65', 'MO1': '#E7CF9F',
                  'MO2': '#E2C68E', 'MO3': '#E2C68E', 'MO4': '#D9B36B', 'MO5': '#D5AB5B', 'MO6': '#D1A349'}

# PATTERNS_TO_CALC = ["ACGAA", "ACGAT", "ACGTA", "ACGTT", "AACGA", "ATCGA", "TACGA", "TTCGA"]
PATTERNS_TO_CALC = ['AACGAA', 'AACGAT', 'AACGTA', 'AACGTT', 'ATCGAA', 'ATCGAT', 'ATCGTA', 'ATCGTT', 'TACGAA',
                    'TACGAT',
                    'TACGTA', 'TACGTT', 'TTCGAA', 'TTCGAT', 'TTCGTA', 'TTCGTT']  # WWCGWW
PATTERNS_TO_CALC = ['CACGAC', 'CACGAG', 'CACGTC', 'CACGTG', 'CTCGAC', 'CTCGAG', 'CTCGTC', 'CTCGTG', 'GACGAC',
                    'GACGAG',
                    'GACGTC', 'GACGTG', 'GTCGAC', 'GTCGAG', 'GTCGTC', 'GTCGTG']  # SWCGWS
PATTERNS_TO_CALC = ['CACGAC', 'CACGAG', 'CACGAA', 'CACGAT', 'GACGAC', 'GACGAG', 'GACGAA', 'GACGAT', 'AACGAC',
                    'AACGAG',
                    'AACGAA', 'AACGAT', 'TACGAC', 'TACGAG', 'TACGAA', 'TACGAT']  # [SW]ACGA[SW]
PATTERNS_TO_CALC = ['CACGAC', 'CACGAG', 'CACGAA', 'CACGAT', 'GACGAC', 'GACGAG', 'GACGAA', 'GACGAT', 'AACGAC',
                    'AACGAG',
                    'AACGAA', 'AACGAT', 'TACGAC', 'TACGAG', 'TACGAA', 'TACGAT']  # [SW]ACGA[SW]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--patients_data', help='Path to files with the calculated data', required=False)
    parser.add_argument('--pmd_dicts', help='Path to location with the pre-calculated PMD dicts',
                        required=False)
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


def get_legend_title(col_name):
    if col_name == 'sublineage':
        return 'Genetic\nsub-lineages:'
    elif col_name == 'region':
        return 'Sampling regions:'


def create_multiple_plots(df, color_by_superimpose_first, to_plot_superimpose_second, superimpose=False,
                          third=None,
                          forth=None, plot_all=None, ylen=14):
    fig, axs = plt.subplots(3, 4, figsize=(22, ylen))
    layout_dict = {'CRC01': (0, 0), 'CRC02': (0, 1), 'CRC04': (0, 2), 'CRC09': (0, 3), 'CRC10': (1, 0),
                   'CRC11': (1, 1),
                   'CRC12': (1, 2), 'CRC13': (1, 3), 'CRC14': (2, 1), 'CRC15': (2, 2)}
    for patient in tqdm(df.patient.unique()):
        if superimpose:
            create_superimposed_graphs(df, patient, axs[layout_dict[patient]], color_by_superimpose_first,
                                       to_plot_superimpose_second, third, forth)
        elif plot_all:
            plot_16(df, patient, axs[layout_dict[patient]])
        else:
            create_colored_graphs(df, color_by_superimpose_first, patient, axs[layout_dict[patient]],
                                  to_plot_superimpose_second)

    for ax in axs.flat:
        ax.set(xlabel="tumor cells (self sorted)", ylabel="avg meth")
    for ax in axs.flat:
        ax.label_outer()
    axs[2, 0].set_visible(False)
    axs[2, 3].set_visible(False)

    if superimpose:
        fig_title = 'cell methylation %s vs %s' % (color_by_superimpose_first, to_plot_superimpose_second)
        path = "cell_methylation_means_superimposed_%s_%s" % (
            color_by_superimpose_first, to_plot_superimpose_second)
        labels = [color_by_superimpose_first, to_plot_superimpose_second]
        custom_lines = [Line2D([0], [0], color='#4872B7'), Line2D([0], [0], color='#8254A2')]
        if third:
            fig_title += " vs %s" % third
            path += "_%s" % third
            labels.append(third)
            custom_lines.append(Line2D([0], [0], color='#3C8A45'))
        if forth:
            fig_title += " vs %s" % forth
            path += "_%s" % forth
            labels.append(forth)
            custom_lines.append(Line2D([0], [0], color='#f376b9'))
        fig.suptitle(fig_title, fontsize=24)
        lgd = fig.legend(custom_lines, labels, loc='center right', fontsize=14)
        plt.savefig(path)
    elif plot_all:
        fig.suptitle('Cell methylation WACGAW', fontsize=24)
        handles, labels = axs[0, 0].get_legend_handles_labels()
        order = list(df.loc[:, PATTERNS_TO_CALC].mean(axis=0).sort_values(ascending=False).index)
        lgd = fig.legend([handles[labels.index(lab)] for lab in order], order, loc='center right',
                         fontsize=14)
        for line in lgd.get_lines():
            line.set_linewidth(2)
        plt.savefig("all_WACGAW.png")
    else:
        fig.suptitle('cell methylation for %s colored by %s' % (
            to_plot_superimpose_second, color_by_superimpose_first),
                     fontsize=24)
        legend_names = sorted(df.loc[:, color_by_superimpose_first].unique())
        custom_lines = [
            Line2D([0], [0], color=get_colours(key, color_by_superimpose_first), ls=' ', marker='s', ms=14)
            for key in
            legend_names]
        title = get_legend_title(color_by_superimpose_first)
        lgd = fig.legend(custom_lines, legend_names, loc='center right', fontsize=14)
        lgd.set_title(title, prop={'size': 14})
        plt.savefig("cell_methylation_means_colored_by_%s_for_%s.png" % (
            color_by_superimpose_first, to_plot_superimpose_second))
    plt.show()


def plot_16(df, patient, axs):
    df = df.sort_values(by=['median'])
    unique_colours = matplotlib.colors.ListedColormap(["#000000", "#FF0000", "#A0A0A0", "#FF7F7F",
                                                       "#2A4BD7", "#1D6914", "#814A19", "#8126C0",
                                                       "#9DAFFF", "#81C57A", "#E9DEBB", "#AD2323",
                                                       "#29D0D0", "#FFEE33", "#FF9233", "#FF30D5"])
    df.loc[df.patient == patient, PATTERNS_TO_CALC].plot(ax=axs, legend=False, linewidth=1, xticks=[],
                                                         colormap=unique_colours)
    axs.set_title("%s" % patient)
    axs.set_ylim([0, 1])


def create_superimposed_graphs(df, patient, axs, first, second, third=None, fourth=None):
    df = df.sort_values(by=['median'])
    add_superimposition_to_graph(axs, df, first, patient, '#4872B7')
    add_superimposition_to_graph(axs, df, second, patient, '#8254A2')
    if third:
        add_superimposition_to_graph(axs, df, third, patient, '#3C8A45')
    if fourth:
        add_superimposition_to_graph(axs, df, fourth, patient, '#f376b9')

    axs.set_title("%s" % patient)
    axs.set_ylim([0, 1])


def add_superimposition_to_graph(axs, df, col_name, patient, color):
    values = df.loc[df.patient == patient, col_name]
    x = [i for i in range(values.shape[0])]
    axs.plot(x, values, color=color, linewidth=1, drawstyle='steps-mid')
    axs.bar(x, values, width=1, alpha=0.1, color='#4872B7')


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


def create_graphs(path):
    df = pd.read_pickle(path)
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
    # create_multiple_plots(df, 'acgaa', 'acgab', superimpose=True)
    # create_multiple_plots(df, 'strong', 'weak', superimpose=True)
    # create_multiple_plots(df, "", "", plot_all=True, ylen=30)
    create_multiple_plots(df, 'AACGAA', 'AACGAT', superimpose=True, third='TACGAA', forth='TACGAT')

    # for patient in ['CRC01', 'CRC10', 'CRC11', 'CRC13']:
    # for patient in df.patient.unique():
    #     create_colored_graphs(df, 'sublineage', patient)
    #     create_colored_graphs(df, 'region', patient)


def get_all_pattern_cells_from_df(df, context_info, pattern):
    pattern_cells = context_info[context_info.str.contains(pattern)]
    return np.mean(df.loc[:, df.columns & pattern_cells.index], axis=1)

def main():
    args = parse_input()

    if args.patients_data:
        create_graphs(args.patients_data)
        return

    nc_files = consts.NC_DROR  # todo: dror

    patients_dict = {}
    if args.pmd_dicts:  # files where to heavy to save, possibly will save in smaller bunches and then use this
        for patient in ["CRC01", "CRC02", "CRC04", "CRC09", "CRC10", "CRC11", "CRC12", "CRC13", "CRC14",
                        "CRC14"]:
            patients_dict[patient] = files_tools.load_compressed_pickle(
                os.path.join(args.pmd_dicts, PMD_DICT_FILE_FORMAT % patient))
    else:
        methylation_folder = args.methylation_folder
        # all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
        all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "CRC01", "*.pkl.zip")))
        global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
        patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                    global_windows_data)

    sublineage = files_tools.load_compressed_pickle(consts.CONVERT_SUBLINEAGE_DROR)


    # # Collect global sequence information from the cpg dict
    cpg_dict = files_tools.get_cpg_context_map(load_with_path=consts.CONTEXT_MAP_FILTERED_LOCAL_DROR)
    strong = {}
    weak = {}
    context_info = {}
    cpg75flank = {}  # todo: dror
    for chr in range(1, 23):
        chr_info = cpg_dict['chr%d' % chr]
        weak['%d' % chr] = chr_info[chr_info[:, -2] == 1][:, 0]
        strong['%d' % chr] = chr_info[chr_info[:, -1] == 1][:, 0]
        context_info['%d' % chr] = pd.Series(chr_info[:, -3], index=chr_info[:, 0]).apply(
            format_cpg_context_map.convert_context_int_to_str)
        cpg75flank["%d" % chr] = chr_info[:, -7]  # todo: dror

    pa_cells = set([])
    all_patients_mean = []
    all_patients_median = []
    for patient in tqdm(patients_dict, desc="get data"):
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            pa_cells |= set(df.index)

        pa_df_meth = pd.DataFrame(index=list(pa_cells))
        pa_df_meth_strong = pd.DataFrame(index=list(pa_cells))
        pa_df_meth_weak = pd.DataFrame(index=list(pa_cells))
        pattern_cells = {}
        for pattern in PATTERNS_TO_CALC:
            pattern_cells[pattern] = pd.DataFrame(index=list(pa_cells))
        for chromosome in tqdm(patients_dict[patient], desc="chromosome data calc %s" % patient):
            df = patients_dict[patient][chromosome]
            pa_df_meth.loc[df.index, chromosome] = np.mean(df, axis=1)
            pa_df_meth_strong.loc[df.index, chromosome] = np.mean(df.loc[:, df.columns & strong[chromosome]],
                                                                  axis=1)
            pa_df_meth_weak.loc[df.index, chromosome] = np.mean(df.loc[:, df.columns & weak[chromosome]],
                                                                axis=1)
            for pattern in PATTERNS_TO_CALC:
                pattern_cells[pattern].loc[df.index, chromosome] = get_all_pattern_cells_from_df(df,
                                                                                                 context_info[
                                                                                                     chromosome],
                                                                                                 pattern)

            nc_meth_avg = get_nc_avg(chromosome, df.index, nc_files)  # todo: dror

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
                                                                                                  sublineage,
                                                                                                  patient))

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
                                                                                          args=(sublineage,
                                                                                                patient))

        ### add the data for the pattern ###
        for pattern in PATTERNS_TO_CALC:
            patient_df_median.loc[:, pattern] = pattern_cells[pattern].median(axis=1)
            patient_df_mean.loc[:, pattern] = pattern_cells[pattern].mean(axis=1)

        all_patients_median.append(patient_df_median)
        all_patients_mean.append(patient_df_mean)

    ### save to file ###
    all_patients_median_df = pd.concat(all_patients_median)
    all_patients_mean_df = pd.concat(all_patients_mean)
    median_path = "all_patients_avg_data_all_dror_median.pickle.zlib"
    mean_path = "all_patients_avg_data_all_dror_mean.pickle.zlib"
    all_patients_median_df.to_pickle(median_path)
    all_patients_mean_df.to_pickle(mean_path)
    create_graphs(median_path)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
