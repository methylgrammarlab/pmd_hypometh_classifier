"""
Functions to handle PMD data
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import consts
from commons import files_tools, data_tools
from covariance import covariance_to_bedgraph

TOP_LOW_PERCENTAGE_TO_REMOVE = 5
PMD_LABEL = "commonPMD"

pmd_dict_mu = {}


def create_pmd_dict(bedgraph_path, output_path):
    """
    Convert PMD_coordinates_hg19 bedgraph file to a dictionary
    :param bedgraph_path: The bedgraph file to convert
    :param output_path: The output path for the dictionary
    :return:
    """
    all_data = pd.read_csv(bedgraph_path, sep="\t")._values
    chr_dict = {}

    chromosomes = np.unique(all_data[:, 0])
    for chromosome in chromosomes:
        chr_lines = all_data[all_data[:, 0] == chromosome]
        pmd_lines = chr_lines[chr_lines[:, -1] == PMD_LABEL]
        chr_dict[chromosome] = []

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
                chr_dict[chromosome].append((global_start, global_end))
                start_new = True

    files_tools.save_pickle(os.path.join(output_path), "pmd_dict.pickle", chr_dict)


def get_pmd_dict(pmd_file_path=consts.PMD_FILE):
    """
    Get the pmd dict
    :param pmd_file_path: different path for the pmd dict
    :return: A dictionary, each chromosome is the key and the values is a list of tuples each represent the
    start and end of pmd
    """
    if pmd_dict_mu != {}:
        return pmd_dict_mu

    pmd_dict = files_tools.load_pickle(pmd_file_path)
    return pmd_dict


def get_pmd_context_map(context_map_file_path=consts.CONTEXT_MAP_FILTERED_LOCAL_DROR,
                        pmd_file_path=consts.PMD_FILE):
    """
    Get the context map information only for PMDs
    :param context_map_file_path: The path for the context map file
    :param pmd_file_path: THe path for the pmd file
    :return: The contaxt map information only for PMDs
    """
    pmd_context_map = {}
    cpg_context_map = files_tools.get_cpg_context_map(load_with_path=context_map_file_path)
    pmd_dict = get_pmd_dict(pmd_file_path=pmd_file_path)

    for chromosome in cpg_context_map:
        chr_context_map = cpg_context_map[chromosome]
        prev_mask = None
        pmd_list = pmd_dict[chromosome]
        for pmd_tuple in pmd_list:
            start, end = pmd_tuple
            pmd_mask = (chr_context_map[:, 0] >= start) & (chr_context_map[:, 0] <= end)
            prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

        pmd_context_map[chromosome] = chr_context_map[prev_mask]

    return pmd_context_map


def filtered_out_non_pmd(df, chromosome, pmd_file=consts.PMD_FILE, add_pmd_index=False):
    """
    Filter out CpG which doesn't include in a PMD
    :param df: The df to work with
    :param chromosome: The chromosome of this df
    :param add_pmd_index: Should we add a pmd index to each cpg
    :param pmd_file: The pmd file
    :return: Same df but filtered based on the PMD
    """
    pmd_dict = get_pmd_dict(pmd_file)

    # Try getting the chromosome from the pmd dict, support chr16 or just 16
    if chromosome in pmd_dict:
        pmd_list = pmd_dict[chromosome]
    elif "chr%s" % chromosome in pmd_dict:
        pmd_list = pmd_dict["chr%s" % chromosome]
    else:
        raise Exception("Chromosome name is invalid ")

    if add_pmd_index:
        df["pmd_index"] = np.nan

    prev_mask = None
    i = 0

    # Trying to support two ways of working with df - one is where the columns is the CpG and the other is
    # that the index is the CpG location
    for pmd_tuple in pmd_list:
        i += 1
        start, end = pmd_tuple
        try:
            pmd_mask = (df.columns >= start) & (df.columns <= end)
        except Exception:
            pmd_mask = (df.index >= start) & (df.index <= end)

        prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

        if add_pmd_index:
            df.loc[pmd_mask, "pmd_index"] = i

    try:
        data = df.loc[:, prev_mask]
    except Exception:
        data = df.loc[prev_mask, :]

    return data


def convert_bedgraph_to_df_with_pmd_filter(bedgraph_path, chromosome, add_pmd_index=False,
                                           pmd_file=consts.PMD_FILE):
    """
    Convert bedgraph to df with a pmd filter
    :param bedgraph_path: The bedgraph file path
    :param chromosome: The chromosome of this bedgraph
    :param add_pmd_index: Should we add a pmd index to each cpg
    :param pmd_file: The path for the pmd file dictionary
    :return: A df of the bedgraph
    """
    covariance_df = files_tools.load_bedgraph(bedgraph_path)
    return filtered_out_non_pmd(df=covariance_df, chromosome=chromosome, pmd_file=pmd_file,
                                add_pmd_index=add_pmd_index)


def get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict, global_windows_data,
                                                    top_low_level_to_remove=TOP_LOW_PERCENTAGE_TO_REMOVE,
                                                    add_pmd_index=False, pmd_file=consts.PMD_FILE):
    # TODO: need to doc this and maybe split , this is really messy
    patients_dict = {}
    for patient in all_files_dict:
        patients_dict[patient] = {}
        for t in all_files_dict[patient]:
            chromosome, file_path = t
            dff = pd.read_pickle(file_path)

            try:
                chromosome = str(chromosome)
                windows_data = global_windows_data[chromosome]
            except Exception:
                chromosome = int(chromosome)
                windows_data = global_windows_data[chromosome]

            try:
                chromosome = str(chromosome)
                covariance_pmd_df = filtered_out_non_pmd(dff, chromosome, add_pmd_index=add_pmd_index,
                                                         pmd_file=pmd_file)
            except:
                chromosome = int(chromosome)
                covariance_pmd_df = filtered_out_non_pmd(dff, chromosome, add_pmd_index=add_pmd_index,
                                                         pmd_file=pmd_file)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.columns >= start) & (covariance_pmd_df.columns <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            _, df = covariance_to_bedgraph.get_region_df(covariance_pmd_df.loc[:, prev_mask],
                                                         sublineage_cells=[],
                                                         sublineage_name=covariance_to_bedgraph.ALL_CANCER)

            if top_low_level_to_remove != 0:
                df = data_tools.remove_extreme_cpgs_by_coverage(df, top_low_level_to_remove)

            patients_dict[patient][chromosome] = df

    return patients_dict


def add_pmd_index_to_df(df, tuple_list):
    """
    Get a pmd index based on the start and end in the tuple
    :param df: A df
    :param tuple_list: A tuple list
    :return: The df with the pmd index
    """
    df["pmd_index"] = np.nan
    i = 0
    for pmd_tuple in tuple_list:
        i += 1
        start, end = pmd_tuple
        pmd_mask = (df.index >= start) & (df.index <= end)
        df.loc[pmd_mask, "pmd_index"] = i

    return df


def get_pmd_index(df, chromosome, pmd_indexes_dict):
    """
    Get the pmd indexes of each cpg in a data frame
    :param df: The dataframe
    :param chromosome: The chromosome
    :param pmd_indexes_dict: A dictionary with with chr as keys and list of tuples representing pmd as values
    :return: The indexes of each cpg in the df
    """
    try:
        chromosome = str(chromosome)
        windows_data = pmd_indexes_dict[chromosome]
    except Exception:
        chromosome = int(chromosome)
        windows_data = pmd_indexes_dict[chromosome]

    return add_pmd_index_to_df(df, windows_data)["pmd_index"]
