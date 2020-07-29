"""
Take all the bedgraphs created from the scWGBS covariance and create a dataset which combine all the
patients and chromosome to one, next we will use this to create the data for a NN
Note: we decided to work with variance and not covariance so this wasn't used
"""

import argparse
import glob
import os
import sys
import warnings

import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools, consts, sequence_utils
from covariance import covariance_to_bedgraph

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedgraph_files', help='Path to bedgraph files', required=True)
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--nc_files', help='Path to nc files', required=True)
    parser.add_argument('--cells_to_use', help='Path to cells to use file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def get_patients_cpgs_index(patients_information):
    """
    Get a list of all the indexes from the different patient matching from the data frames
    :param patients_information: A list with all the patient, each value here is a tuple with a df in the
    second value
    :return: A set of all the indexes available from all the patients
    :rtype: set
    """
    cpg_locations = None
    for df_tuple in patients_information:
        _, cov_df = df_tuple
        if cpg_locations is not None:
            cpg_locations = set(cov_df.index.values) | cpg_locations
        else:
            cpg_locations = set(cov_df.index.values)
    return cpg_locations


def get_covariance_per_chromosome_filtered(all_bedgraphs, boundaries_dict):
    """
    Arrange the df from all the bedgraph after filtering only for PMD between the boundaries
    :param all_bedgraphs: List of all bedgraphs files
    :param boundaries_dict: Information about boundaries for each chromosome
    :return: A dictionary of chromosome:(patient, df after filtered)
    """
    cov_dict = {}
    for chromosome in all_bedgraphs:
        cov_dict[chromosome] = []
        chr_window_data = boundaries_dict[chromosome]
        for bedgraph in all_bedgraphs[chromosome]:
            patient, _ = consts.PATIENT_CHR_NAME_RE.findall(bedgraph)[0]
            covariance_pmd_df = handle_pmds.convert_bedgraph_to_df_with_pmd_filter(bedgraph_path=bedgraph,
                                                                                   chromosome=chromosome,
                                                                                   add_pmd_index=True)
            df_in_boundaries = handle_pmds.filter_df_based_on_tuple_list(df=covariance_pmd_df,
                                                                         boundaries_list=chr_window_data)
            cov_dict[chromosome].append((patient, df_in_boundaries))

    return cov_dict


def cpg_meth_in_cells(patient, chromosome, cpgs_indexes, meth_files_path, sublineage_name):
    """
    Get the methylation ratio of a patient in a chromosome for specific set of cells\sublineage_name
    :param patient: The patient we are looking at
    :param chromosome: The chromosome we are looking at
    :param cpgs_indexes: The of cpgs
    :param meth_files_path: The path for the patient files
    :param sublineage_name: The sublineage name or cell type name
    :return: A df with the indexes matching those cells type
    """
    methylation_file_path = os.path.join(meth_files_path, patient, consts.SCWGBS_FILE_FORMAT %
                                         (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    _, df = covariance_to_bedgraph.get_region_df(df, sublineage_cells=[], sublineage_name=sublineage_name)

    return df.loc[cpgs_indexes]


def get_nc_avg(chromosome, cpg_indexes, nc_meth_files):
    """
    Get the average values of methylation in original cells based on a file already created
    :param chromosome: The chromsome we are looking at
    :param cpg_indexes: The cpg indexes
    :param nc_meth_files: The nc files path
    :return: The average methylation in those cells
    """
    nc_file = glob.glob(os.path.join(nc_meth_files, "*%s.dummy*" % chromosome))[0]
    df = pd.read_pickle(nc_file)
    mdf = df.mean(axis=1)
    return mdf.loc[cpg_indexes]


def main():
    args = parse_input()

    bedgraphs_path = files_tools.get_files_to_work(args.bedgraph_files, pattern=BEDGRPH_FILE_FORMAT)
    all_bedgraphs = files_tools.convert_paths_list_to_chromosome_based_dict(bedgraphs_path)
    windows_dict = files_tools.load_compressed_pickle(args.windows_file)

    cov_dict = get_covariance_per_chromosome_filtered(all_bedgraphs, windows_dict)

    chromosomes_list = []
    for chromosome in cov_dict:
        cpg_locations = get_patients_cpgs_index(cov_dict[chromosome])
        indexes = list(cpg_locations)
        indexes.sort()

        # Create df with location, chromosome and sequence
        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = indexes
        chromosome_df["chromosome"] = chromosome
        chromosome_df["sequence"] = sequence_utils.get_sequences_for_cpgs(indexes, chromosome)

        chromosome_df = chromosome_df.set_index("location")

        nc_meth_avg = get_nc_avg(chromosome, chromosome_df.index, args.nc_files)
        chromosome_df.loc[nc_meth_avg.index, "nc_avg"] = nc_meth_avg

        for df_tuple in cov_dict[chromosome]:
            patient, patient_df = df_tuple
            patient_num = patient[-2:]
            cpgs_indexes = patient_df.index

            chromosome_df.loc[patient_df.index, "cov%s" % patient_num] = patient_df["coverage"]
            chromosome_df.loc[patient_df.index, "pmd_index"] = patient_df["pmd_index"]

            cancer_meth = cpg_meth_in_cells(patient, chromosome, cpgs_indexes, args.methylation_folder,
                                            sublineage_name=covariance_to_bedgraph.ALL_CANCER)
            chromosome_df.loc[cancer_meth.index, "meth%s" % patient_num] = cancer_meth.mean()

            original_meth = cpg_meth_in_cells(patient, chromosome, cpgs_indexes, args.methylation_folder,
                                              sublineage_name=covariance_to_bedgraph.ONLY_NC)
            chromosome_df.loc[original_meth.index, "orig_meth%s" % patient_num] = original_meth.mean()

            cancer_variance = cpg_meth_in_cells(patient, chromosome, cpgs_indexes, args.methylation_folder,
                                                sublineage_name=covariance_to_bedgraph.ALL_CANCER)
            chromosome_df.loc[cancer_variance.index, "var%s" % patient[-2:]] = cancer_variance.var()

        chromosomes_list.append(chromosome_df.reset_index())

    all_chromosome_df = pd.concat(chromosomes_list)
    all_chromosome_df.to_pickle(os.path.join(args.output_folder, "valid_cpg.pkl"))


if __name__ == '__main__':
    main()
