import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm

from variance.get_valid_cpg import get_nc_avg_all

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from covariance import covariance_to_bedgraph
import commons.consts as consts

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
PMD_DICT_FILE_FORMAT = "%s_pmd_dicts.pkl.zip"
SEQ_SIZE = 150
TOP_LOW_PERCENTAGE_TO_REMOVE = 5
SOLO = 0

STRONG_WEAK = {"strong": "SCGS", "weak": "WCGW"}

LESION_ABBREVIATION = {'ML': 'Liver Metastasis', "PT": "Primary Tumor", "LN": "Lymph Node\nMetastasis",
                       "MP": "Post-treatment\nLiver Metastasis", "NC": "Normal Cell", "MO": "Omental Metastasis"}


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


def get_all_pattern_cells_from_df(df, context_info, pattern):
    pattern_cells = context_info[context_info.str.contains(pattern)]
    return np.mean(df.loc[:, df.columns & pattern_cells.index], axis=1)


def get_total_meth_avg(all_files_dict, nc_meth):
    patients_dict = {}
    for patient in all_files_dict:
        patients_dict[patient] = {}
        for t in all_files_dict[patient]:
            chromosome, file_path = t
            df = pd.read_pickle(file_path)
            df = df.loc[:, df.columns & nc_meth[chromosome]]
            patients_dict[patient][chromosome] = np.nanmean(df, axis=1, dtype=np.float64)
    return patients_dict


def main():
    args = parse_input()

    nc_files = consts.NC_LIOR

    methylation_folder = args.methylation_folder
    all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                global_windows_data)

    sublineage = files_tools.load_compressed_pickle(consts.CONVERT_SUBLINEAGE_LIOR)

    # # Collect global sequence information from the cpg dict
    cpg_dict = files_tools.get_cpg_context_map(load_with_path=consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI_LIOR)
    strong = {}
    weak = {}
    cpg75flank_solo = {}
    nc_meth = {}
    for chr in range(1, 23):
        chr_info = cpg_dict['chr%d' % chr]
        weak['%d' % chr] = chr_info[chr_info[:, -2] == 1][:, 0]
        strong['%d' % chr] = chr_info[chr_info[:, -1] == 1][:, 0]
        cpg75flank_solo["%d" % chr] = chr_info[chr_info[:, -7] == SOLO][:, 0]
        nc_meth_avg = get_nc_avg_all(chr, nc_files)
        nc_meth["%d" % chr] = nc_meth_avg[nc_meth_avg > 0.5].index

    total_meth_dict = get_total_meth_avg(all_files_dict, nc_meth)

    pa_cells = set([])
    all_patients_data = []
    for patient in tqdm(patients_dict, desc="get data"):
        for chromosome in patients_dict[patient]:
            df = patients_dict[patient][chromosome]
            df = df.loc[:, df.columns & nc_meth[chromosome]]
            # df = df.loc[:, df.columns & cpg75flank_solo[chromosome] & nc_meth[chromosome]]
            pa_cells |= set(df.index)

        for chromosome in tqdm(patients_dict[patient], desc="chromosome data calc %s" % patient):
            pa_df_meth = pd.DataFrame(index=list(pa_cells))
            df = patients_dict[patient][chromosome]

            pa_df_meth.loc[df.index, 'patient'] = patient
            pa_df_meth.loc[df.index, 'chr'] = chromosome
            pa_df_meth.loc[:, 'region'] = pa_df_meth.index
            pa_df_meth.loc[:, 'region'] = pa_df_meth.loc[:, 'region'].apply(get_region)
            pa_df_meth.loc[:, 'lesion'] = pa_df_meth.loc[:, 'region'].apply(lambda x: x[:2])
            pa_df_meth.loc[:, 'sublineage'] = pa_df_meth.index
            pa_df_meth.loc[:, 'sublineage'] = pa_df_meth.loc[:, 'sublineage'].apply(get_sublineage,
                                                                                    args=(sublineage, patient))
            pa_df_meth.loc[df.index, 'total'] = total_meth_dict[patient][chromosome]
            pa_df_meth.loc[df.index, 'PMD'] = np.nanmean(df, axis=1, dtype=np.float64)
            pa_df_meth.loc[df.index, 'PMD_solo'] = np.nanmean(df.loc[:, df.columns & cpg75flank_solo[chromosome]], axis=1, dtype=np.float64)
            pa_df_meth.loc[df.index, 'PMD_strong'] = np.nanmean(df.loc[:, df.columns & strong[chromosome]], axis=1, dtype=np.float64)
            pa_df_meth.loc[df.index, 'PMD_weak'] = np.nanmean(df.loc[:, df.columns & weak[chromosome]], axis=1, dtype=np.float64)
            pa_df_meth.loc[df.index, 'PMD_solo_strong'] = np.nanmean(
                df.loc[:, df.columns & strong[chromosome] & cpg75flank_solo[chromosome]], axis=1, dtype=np.float64)
            pa_df_meth.loc[df.index, 'PMD_solo_weak'] = np.nanmean(
                df.loc[:, df.columns & weak[chromosome] & cpg75flank_solo[chromosome]], axis=1, dtype=np.float64)
            all_patients_data.append(pa_df_meth)

    print("finished processing")
    ### save to file ###
    all_patients_data_df = pd.concat(all_patients_data)
    print("concat done")
    path = "all_patients_chr_data.pickle.zlib"
    csv_path = "all_patients_chr_data.csv"
    all_patients_data_df.to_pickle(path)
    print("saved to pickle")
    all_patients_data_df.to_csv(csv_path)
    print("saved to csv")


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()


    # # df = pd.read_pickle(r"C:\Users\liorf\Documents\filtered_by_bl_and_cpgi\CRC01\all_cpg_ratios_CRC01_chr16.dummy.pkl.zip")
    # # print("loaded")
    # # df = df.T
    # # df = df.reset_index()
    # #
    # # # import feather
    # # # import pyarrow
    # path = 'solo_nc/all_cpg_ratios_CRC01_chr16.feather'
    # # # df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["a", "b", "c"])
    # # df.to_feather(path)
    # dff = pd.read_feather(path)
    # print(dff)
    # # pyarrow.feather.write_dataframe(df, path)
    # # feather.write_dataframe(df, path)
    # # feather.read_dataframe(path)
    #
    # return
    #
    # df.to_csv("solo_nc/all_cpg_ratios_CRC01_chr16.csv", index_label='cpg_loc')
    # return