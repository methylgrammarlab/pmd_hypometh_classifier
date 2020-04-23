import glob
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd

CHR_NAME = re.compile("(chr\d+)")
PMD_LABEL = "commonPMD"
pmd_dict_mu = {}

sys.path.append(os.getcwd())
from commons import consts
from commons import files_tools
from covariance import covariance_to_bedgraph

TOP_LOW_PERCENTAGE_TO_REMOVE = 5


def create_pmd_dict(bed_file, output_path):
    all_data = pd.read_csv(bed_file, sep="\t")._values
    chr_dict = {}

    chromose = np.unique(all_data[:, 0])
    for chro in chromose:
        chr_lines = all_data[all_data[:, 0] == chro]
        pmd_lines = chr_lines[chr_lines[:, -1] == PMD_LABEL]
        chr_dict[chro] = []

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
                chr_dict[chro].append((global_start, global_end))
                start_new = True

    with open(os.path.join(output_path, "pmd_dict.pickle"), "wb") as pmd_dict_file:
        pickle.dump(chr_dict, pmd_dict_file)


def split_files(input_path, output_path, dict_path, patients):
    with open(dict_path, "rb") as pmd_dict_f:
        pmd_dict = pickle.load(pmd_dict_f)

    for patient in patients:
        files = glob.glob(os.path.join(input_path, patients, "*_all_data.csv.gzip"))
        for f in files:
            chr_name = CHR_NAME.findall(f)[0]
            csv_data = pd.read_csv(f, compression="gzip")
            pmd_list = pmd_dict[chr_name]
            for pmd_tuple in pmd_list:
                start, end = pmd_tuple
                pmd_mask = (csv_data['locations'] > start) & (csv_data['locations'] < end)
                pmd = csv_data.loc[pmd_mask, :]
                output_data = os.path.join(output_path, patient, "%s_pmd_%s_%s.csv.gzip" % (chr_name,
                                                                                            start, end))
                if not os.path.exists(os.path.dirname(output_data)):
                    os.mkdir(os.path.dirname(output_data))

                pmd.to_csv(output_data, compression="gzip")


def read_pmd_dict(file_path=consts.PMD_FILE):
    if pmd_dict_mu != {}:
        return pmd_dict_mu

    with open(file_path, "rb") as pmd_dict_f:
        pmd_dict = pickle.load(pmd_dict_f)

        return pmd_dict


def get_pmd_context_map():
    pmd_context_map = {}
    cpg_context_map = files_tools.get_cpg_context_map(load_with_path=consts.CONTEXT_MAP_FILTERED_LOCAL_DROR)
    pmd_dict = read_pmd_dict(consts.PMD_FILE_LOCAL_DROR)

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


def get_pmd_df(df, chromosome, add_pmd_index=False):
    pmd_dict = read_pmd_dict(consts.PMD_FILE_LOCAL_DROR)

    if chromosome in pmd_dict:
        pmd_list = pmd_dict[chromosome]
    elif "chr%s" % chromosome in pmd_dict:
        pmd_list = pmd_dict["chr%s" % chromosome]
    else:
        raise Exception("Chromosome name is invalid ")

    prev_mask = None
    i = 0

    for pmd_tuple in pmd_list:
        i += 1
        start, end = pmd_tuple
        try:
            pmd_mask = (df.columns >= start) & (df.columns <= end)
        except Exception:
            pmd_mask = (df.index >= start) & (df.index <= end)

        prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

    try:
        data = df.loc[:, prev_mask]
    except Exception:
        data = df.loc[prev_mask, :]

    return data


def get_covariance_pmd_df(bedgraph_path, chromosome, add_pmd_index=False):
    covariance_df = files_tools.load_badgraph_to_df(bedgraph_path)
    pmd_dict = read_pmd_dict(consts.PMD_FILE_LOCAL_DROR)

    if chromosome in pmd_dict:
        pmd_list = pmd_dict[chromosome]
    elif "chr%s" % chromosome in pmd_dict:
        pmd_list = pmd_dict["chr%s" % chromosome]
    else:
        raise Exception("Chromosome name is invalid ")

    if add_pmd_index:
        covariance_df["pmd_index"] = np.nan

    prev_mask = None
    i = 0
    for pmd_tuple in pmd_list:
        i += 1
        start, end = pmd_tuple
        pmd_mask = np.logical_and(covariance_df.index >= start, covariance_df.index <= end)
        prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

        if add_pmd_index:
            covariance_df.loc[pmd_mask, "pmd_index"] = i

    return covariance_df.loc[prev_mask, :]


def get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict, global_windows_data,
                                                    top_low_level_to_remove=TOP_LOW_PERCENTAGE_TO_REMOVE,
                                                    add_pmd_index=False):
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
                covariance_pmd_df = get_pmd_df(dff, chromosome, add_pmd_index)
            except:
                chromosome = int(chromosome)
                covariance_pmd_df = get_pmd_df(dff, chromosome, add_pmd_index)

            prev_mask = None
            for pmd_tuple in windows_data:
                start, end = pmd_tuple
                pmd_mask = (covariance_pmd_df.columns >= start) & (covariance_pmd_df.columns <= end)
                prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

            _, df = covariance_to_bedgraph.get_region_df(covariance_pmd_df.loc[:, prev_mask],
                                                         sublineage_cells=[],
                                                         sublineage_name=covariance_to_bedgraph.ALL_CANCER)

            ##
            # Get the coverage of cpg
            ##
            cpg_coverage = np.sum(~pd.isnull(df), axis=0)
            cpg_s = cpg_coverage.shape[0]
            n_to_remove = int(cpg_s * top_low_level_to_remove / 100)
            cpg_to_keep = cpg_coverage.index[n_to_remove:-n_to_remove]
            df = df[cpg_to_keep]  # this remove the 5% lower and top

            patients_dict[patient][chromosome] = df

    return patients_dict


def get_pmd_index(df, chromosome, global_windows_data):
    try:
        chromosome = str(chromosome)
        windows_data = global_windows_data[chromosome]
    except Exception:
        chromosome = int(chromosome)
        windows_data = global_windows_data[chromosome]

    df["pmd_index"] = np.nan
    i = 0
    for pmd_tuple in windows_data:
        i += 1
        start, end = pmd_tuple
        pmd_mask = (df.index >= start) & (df.index <= end)
        df.loc[pmd_mask, "pmd_index"] = i

    return df["pmd_index"]
