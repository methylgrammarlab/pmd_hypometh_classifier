import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm

import commons.data_tools
import commons.utils
from commons.utils import get_nc_avg

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from format_files import format_cpg_context_map
import commons.consts as consts

SOLO = 0

# PATTERNS_TO_CALC = ["ACGAA", "ACGAT", "ACGTA", "ACGTT", "AACGA", "ATCGA", "TACGA", "TTCGA"]
# PATTERNS_TO_CALC = ['AACGAA', 'AACGAT', 'AACGTA', 'AACGTT', 'ATCGAA', 'ATCGAT', 'ATCGTA', 'ATCGTT', 'TACGAA', 'TACGAT',
#                     'TACGTA', 'TACGTT', 'TTCGAA', 'TTCGAT', 'TTCGTA', 'TTCGTT']  # WWCGWW
# PATTERNS_TO_CALC = ['ACGA', 'ACGC', 'ACGG', 'ACGT',
#                     'CCGA', 'CCGC', 'CCGG', 'CCGT',
#                     'GCGA', 'GCGC', 'GCGG', 'GCGT',
#                     'TCGA', 'TCGC', 'TCGG', 'TCGT']  # [SW]CG[SW]
PATTERNS_TO_CALC = []


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    # parser.add_argument('--patients_data', help='Path to files with the calculated data', required=False)
    parser.add_argument('--pmd_dicts', help='Path to location with the pre-calculated PMD dicts', required=False)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def get_sublineage(cell, sublineage_info, patient):
    """
    :param cell: The cell name
    :param sublineage_info: Dictionary of cell names -> sublineage
    :param patient: Patient name
    :return: The cells sublineage if there is, undefined otherwise
    """
    try:
        return sublineage_info[patient + '_' + cell]
    except:
        return "undefined"


def get_region(cell):
    """
    :param cell: The cell name
    :return: The cells region
    """
    return cell.split('_')[0]


def get_all_pattern_cells_from_df(patient_chromosome_df, context_info, pattern):
    """
    For every cell in the DF, calculates the mean of the cels that match the given pattern
    :param patient_chromosome_df: Patients methylation df (for chromosome)
    :param context_info: df with context information
    :param pattern: The pattern to look for
    :return: df with the mean for every cell of the CpGs that fit the pattern
    """
    pattern_cells = set([])
    pattern_cells |= set(context_info[context_info.str.contains(pattern)].index)
    return np.mean(patient_chromosome_df.loc[:, patient_chromosome_df.columns & pattern_cells], axis=1)


def get_patient_df_dict(all_file_paths):
    """
    Loads all the methylation DF and returns them as a dictionary
    :param all_file_paths: List of all the methylation files
    :return: A dictionary were the keys are the patient and the values a dictionary of chromosome and the methylation df
    """
    patient_df_dict = {}
    for file_path in tqdm(all_file_paths):
        patient, chromosome = consts.DATA_FILE_SCWGBS_RE.findall(file_path)[0]
        if patient not in patient_df_dict:
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
    :param patient_df_dict: Dictionary that for each patient holds a dictionary of the DFs of all it's chromosomes
    :param boundaries_data: The boundaries data
    :param perc_of_cpg_to_remove_based_on_coverage: Percentage of CpGs to remove from top and low coverage
    :return: A new dictionary of the same format, with the filtered DFs
    """
    filtered_dict = {patient: {} for patient in patient_df_dict}
    for patient in patient_df_dict:
        for chromosome in patient_df_dict[patient]:
            df = patient_df_dict[patient][chromosome]
            # filter oyt non-PMDs
            filtered_df = handle_pmds.filtered_out_non_pmd(df, chromosome, add_pmd_index=False,
                                                           pmd_file=consts.PMD_FILE_LOCAL_LIOR)
            # filter based on given boundaries
            filtered_df = commons.utils.filter_df_based_on_tuple_list(filtered_df,
                                                                      boundaries_data[chromosome[
                                                                                      3:]])  # todo format chromosome name
            # filter out low and high coverage cpg
            filtered_df = commons.utils.remove_extreme_cpgs_by_coverage(
                filtered_df, top_low_level_to_remove=perc_of_cpg_to_remove_based_on_coverage)
            filtered_dict[patient][chromosome] = filtered_df
    return filtered_dict


def get_filtered_patient_df_dict(methylation_folder, windows_file):
    """
    Loads the data and filters it.
    :param methylation_folder: Path to the folder with the methylation files
    :param windows_file: The windows data
    :return: A dictionary were the keys are the patient and the values a dictionary of chromosome and the methylation df
    """
    all_files = files_tools.get_files_to_work(methylation_folder, pattern=os.path.join("CRC09", "*.pkl.zip"))
    boundaries_data = files_tools.load_compressed_pickle(windows_file)

    patients_dict = get_patient_df_dict(all_files)
    patients_dict = filter_patient_df_dict(patients_dict, boundaries_data)
    return patients_dict


def create_methylation_info_dfs(cpg_context_dict, nc_files, patients_dict, sublineage):
    """
    Creates df of the mean and median methylation for all patients.
    :param cpg_context_dict: Dictionary, for each chromosome holds a df with context information
    :param nc_files: The nc files path
    :param patients_dict: Dictionary with the methylation df for every patient, chromosome
    :param sublineage: Dictionary mapping cell names to their sublineage
    :return: Lists of all the patients mean and median methylation DFs
    """
    context_info, cpg75flank, nc_meth, strong_indices, weak_indices = find_cpg_indices_for_feature(cpg_context_dict,
                                                                                                   nc_files)

    patient_cell_names = set([])
    all_patients_mean = []
    all_patients_median = []
    all_patients_total_mean = []
    for patient in tqdm(patients_dict, desc="get data"):
        patient_cell_names = get_patient_cells(cpg75flank, nc_meth, patient, patient_cell_names, patients_dict)

        mean_methylation_df, pattern_mean_methylation_df_dict, strong_mean_methylation_df, weak_mean_methylation_df, \
        total_mean, strong_mean, weak_mean = \
            get_mean_methylation_dfs(context_info, patient, patient_cell_names, patients_dict, strong_indices,
                                     weak_indices)

        patient_df_mean, patient_df_median, patient_df_all_mean = create_mean_median_data(mean_methylation_df,
                                                                     pattern_mean_methylation_df_dict,
                                                                     strong_mean_methylation_df,
                                                                     weak_mean_methylation_df, patient, sublineage,
                                                                     total_mean, strong_mean, weak_mean)

        all_patients_median.append(patient_df_median)
        all_patients_mean.append(patient_df_mean)
        all_patients_total_mean.append(patient_df_all_mean)

    return all_patients_mean, all_patients_median, all_patients_total_mean


def create_mean_median_data(mean_methylation_df, pattern_mean_methylation_df_dict, strong_mean_methylation_df,
                            weak_mean_methylation_df, patient, sublineage, total_mean, strong_mean, weak_mean):
    """
    From the DFs (all, pattern, strong, weak) containing the mean of each cell per chromosome create two DFs, one for
    the mean and the other for the median for every cell across all chromosomes with the different features.
    :param mean_methylation_df: df with the mean methylation for every cell and chromosome for all CpGs
    :param pattern_mean_methylation_df_dict: df with the mean methylation for every cell and chromosome for CpGs that match patterns
    :param strong_mean_methylation_df: df with the mean methylation for every cell and chromosome for all strong CpGs
    :param weak_mean_methylation_df: df with the mean methylation for every cell and chromosome for all weak CpGs
    :param patient: patient name
    :param sublineage: Dictionary of cell names -> sublineage
    :return: Two DFs, one for the mean and the other for the median for every cell across all chromosomes with the
    different features.
    """
    ### save median data ###
    patient_df_median = pd.DataFrame(index=mean_methylation_df.index)
    patient_df_median.loc[:, 'median'] = mean_methylation_df.median(axis=1)
    patient_df_median.loc[:, 'strong'] = strong_mean_methylation_df.median(axis=1)
    patient_df_median.loc[:, 'weak'] = weak_mean_methylation_df.median(axis=1)
    patient_df_median.loc[:, 'patient'] = patient
    patient_df_median.loc[:, 'region'] = mean_methylation_df.index
    patient_df_median.loc[:, 'region'] = patient_df_median.loc[:, 'region'].apply(get_region)
    patient_df_median.loc[:, 'lesion'] = patient_df_median.loc[:, 'region'].apply(lambda x: x[:2])
    patient_df_median.loc[:, 'sublineage'] = mean_methylation_df.index
    patient_df_median.loc[:, 'sublineage'] = patient_df_median.loc[:, 'sublineage'].apply(get_sublineage,
                                                                                          args=(
                                                                                              sublineage, patient))
    ### save mean data ###
    patient_df_mean = pd.DataFrame(index=mean_methylation_df.index)
    patient_df_mean.loc[:, 'mean'] = mean_methylation_df.mean(axis=1)
    patient_df_mean.loc[:, 'strong'] = strong_mean_methylation_df.mean(axis=1)
    patient_df_mean.loc[:, 'weak'] = weak_mean_methylation_df.mean(axis=1)
    patient_df_mean.loc[:, 'patient'] = patient
    patient_df_mean.loc[:, 'region'] = mean_methylation_df.index
    patient_df_mean.loc[:, 'region'] = patient_df_mean.loc[:, 'region'].apply(get_region)
    patient_df_mean.loc[:, 'lesion'] = patient_df_mean.loc[:, 'region'].apply(lambda x: x[:2])
    patient_df_mean.loc[:, 'sublineage'] = mean_methylation_df.index
    patient_df_mean.loc[:, 'sublineage'] = patient_df_mean.loc[:, 'sublineage'].apply(get_sublineage,
                                                                                      args=(sublineage, patient))
    ### save total mean data ###
    patient_df_all_mean = pd.DataFrame(index=total_mean.index)
    patient_df_all_mean.loc[:, 'mean'] = total_mean
    patient_df_all_mean.loc[:, 'strong'] = strong_mean
    patient_df_all_mean.loc[:, 'weak'] = weak_mean
    patient_df_all_mean.loc[:, 'patient'] = patient
    patient_df_all_mean.loc[:, 'region'] = total_mean.index
    patient_df_all_mean.loc[:, 'region'] = patient_df_mean.loc[:, 'region'].apply(get_region)
    patient_df_all_mean.loc[:, 'lesion'] = patient_df_mean.loc[:, 'region'].apply(lambda x: x[:2])
    patient_df_all_mean.loc[:, 'sublineage'] = total_mean.index
    patient_df_all_mean.loc[:, 'sublineage'] = patient_df_mean.loc[:, 'sublineage'].apply(get_sublineage,
                                                                                      args=(sublineage, patient))
    ### add the data for the pattern ###
    for pattern in PATTERNS_TO_CALC:
        patient_df_median.loc[:, pattern] = pattern_mean_methylation_df_dict[pattern].median(axis=1)
        patient_df_mean.loc[:, pattern] = pattern_mean_methylation_df_dict[pattern].mean(axis=1)
    return patient_df_mean, patient_df_median, patient_df_all_mean


def get_mean_methylation_dfs(context_info, patient, patient_cell_names, patients_dict, strong_indices, weak_indices):
    """
    Calculate the mean methyaltion for every cell and chromosome for CpGs fitting different features
    :param context_info: Df holding context information
    :param patient: Ptient name
    :param patient_cell_names: List of the cell names
    :param patients_dict: Dictionary with the methylation df for every patient, chromosome
    :param strong_indices: Indices of strong CpGs
    :param weak_indices: Indices of weak CpGs
    :return: DFs (all, pattern, strong, weak) containing the mean of each cell per chromosome
    """
    mean_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    strong_mean_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    weak_mean_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    sum_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    num_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    strong_sum_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    strong_num_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    weak_sum_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    weak_num_methylation_df = pd.DataFrame(index=list(patient_cell_names))
    pattern_mean_methylation_df_dict = {}
    for pattern in PATTERNS_TO_CALC:
        pattern_mean_methylation_df_dict[pattern] = pd.DataFrame(index=list(patient_cell_names))
    for chromosome in tqdm(patients_dict[patient], desc="chromosome data calc %s" % patient):
        patient_chromosome_df = patients_dict[patient][chromosome]

        mean_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.mean(patient_chromosome_df, axis=1)
        sum_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.sum(patient_chromosome_df, axis=1)
        num_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.sum(np.isfinite(patient_chromosome_df),
                                                                                 axis=1)

        strong_mean_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.mean(
            patient_chromosome_df.loc[:, patient_chromosome_df.columns & strong_indices[chromosome]], axis=1)
        strong_sum_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.sum(
            patient_chromosome_df.loc[:, patient_chromosome_df.columns & strong_indices[chromosome]], axis=1)
        strong_num_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.sum(np.isfinite(
            patient_chromosome_df.loc[:, patient_chromosome_df.columns & strong_indices[chromosome]]), axis=1)

        weak_mean_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.mean(
            patient_chromosome_df.loc[:, patient_chromosome_df.columns & weak_indices[chromosome]], axis=1)
        weak_sum_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.sum(
            patient_chromosome_df.loc[:, patient_chromosome_df.columns & weak_indices[chromosome]], axis=1)
        weak_num_methylation_df.loc[patient_chromosome_df.index, chromosome] = np.sum(np.isfinite(
            patient_chromosome_df.loc[:, patient_chromosome_df.columns & weak_indices[chromosome]]), axis=1)

        for pattern in PATTERNS_TO_CALC:
            pattern_mean_methylation_df_dict[pattern].loc[
                patient_chromosome_df.index, chromosome] = get_all_pattern_cells_from_df(patient_chromosome_df,
                                                                                         context_info[chromosome],
                                                                                         pattern)
    mean_methylation_df.reset_index().melt(id_vars='index', var_name='chr', value_name='tot_avg')
    total_mean = np.sum(sum_methylation_df, axis=1) / np.sum(num_methylation_df, axis=1)
    strong_mean = np.sum(strong_sum_methylation_df, axis=1) / np.sum(strong_num_methylation_df, axis=1)
    weak_mean = np.sum(weak_sum_methylation_df, axis=1) / np.sum(weak_num_methylation_df, axis=1)
    return mean_methylation_df, pattern_mean_methylation_df_dict, strong_mean_methylation_df, weak_mean_methylation_df, \
           total_mean, strong_mean, weak_mean


def get_patient_cells(cpg75flank, nc_meth, patient, patient_cell_names, patients_dict):
    """
    Find cell names
    :param cpg75flank:
    :param nc_meth:
    :param patient:
    :param patient_cell_names:
    :param patients_dict:
    :return: Set of all the patients cell names
    """
    for chromosome in patients_dict[patient]:
        patient_chromosome_df = patients_dict[patient][chromosome]
        solo_nc_patient_chromosome_df = patient_chromosome_df.loc[:,
                                        patient_chromosome_df.columns & cpg75flank[chromosome] & nc_meth[chromosome]]
        patient_cell_names |= set(solo_nc_patient_chromosome_df.index)
    return patient_cell_names


def find_cpg_indices_for_feature(cpg_context_dict, nc_files):
    """
    Find indices of all CpGs of different features - strong, weak, solo, methylated and the context info
    :param cpg_context_dict: Dictionary, for each chromosome holds a df with context information
    :param nc_files: The nc files path
    :return: Indices of all CpGs - strong, weak, solo, methylated and the context info
    """
    strong = {}
    weak = {}
    cpg75flank = {}
    nc_meth = {}
    context_info = {}
    for chr in range(1, 23):  # todo Dror / Lior - format chromosome name in the CONTEXT_MAP_FILTERED_NO_BL_CPGI
        chr_info = cpg_context_dict['chr%d' % chr]
        strong['chr%d' % chr] = chr_info[chr_info[:, -1] == 1][:, 0]
        weak['chr%d' % chr] = chr_info[chr_info[:, -2] == 1][:, 0]
        cpg75flank["chr%d" % chr] = chr_info[chr_info[:, -7] == SOLO][:, 0]
        print(chr)
        nc_meth_avg = get_nc_avg(chr, chr_info[:, 0], nc_files)
        nc_meth["chr%d" % chr] = nc_meth_avg[nc_meth_avg > 0.5].index
        context_info['chr%d' % chr] = pd.Series(chr_info[:, -3], index=chr_info[:, 0]).apply(
            format_cpg_context_map.convert_context_int_to_str)
    return context_info, cpg75flank, nc_meth, strong, weak


def save_to_file(all_patients_mean, all_patients_median, output_folder, all_patients_total_mean):
    """
    Saves both the mean DF and the median df to a pickle and a csv
    :param output_folder:
    :param all_patients_mean: List of all the patients mean methylation DFs
    :param all_patients_median: List of all the patients median methylation DFs
    :param output_folder: Path of the output folder
    """
    base_path = os.path.join(output_folder, "avg_data_all_NCGN_solo_nc_")
    median_path = base_path + "median.pickle.zlib"
    mean_path = base_path + "mean.pickle.zlib"

    all_patients_median_df = pd.concat(all_patients_median)
    all_patients_mean_df = pd.concat(all_patients_mean)
    all_patients_total_mean_df = pd.concat(all_patients_total_mean)

    # all_patients_median_df.to_pickle(median_path)
    # all_patients_mean_df.to_pickle(mean_path)
    # all_patients_median_df.to_csv(base_path + "median.csv")
    all_patients_total_mean_df.to_csv(base_path + "total_mean_CRC09.csv")


def main():
    args = parse_input()

    sublineage = files_tools.load_compressed_pickle(consts.CONVERT_SUBLINEAGE_LIOR)
    cpg_context_dict = files_tools.get_cpg_context_map(load_with_path=consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI_LIOR)

    patients_dict = get_filtered_patient_df_dict(args.methylation_folder, args.windows_file)
    all_patients_mean, all_patients_median, all_patients_total_mean = create_methylation_info_dfs(cpg_context_dict, consts.NC_LIOR,
                                                                         patients_dict, sublineage)
    save_to_file(all_patients_mean, all_patients_median, args.output_folder, all_patients_total_mean)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
