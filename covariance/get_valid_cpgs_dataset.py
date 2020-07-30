"""
Take all the bedgraphs created from the scWGBS covariance and create a dataset which combine all the
patients and chromosome to one, next we will use this to create the data for a NN
Note: we decided to work with variance and not covariance so this wasn't used
"""

import argparse
import os
import sys
import warnings

import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools, consts, sequence_tools, utils

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedgraph_files', help='Path to bedgraph files', required=True)
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--nc_files', help='Path to nc files', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def get_bedgraph_to_df_per_chromosome_filtered(all_bedgraphs, boundaries_dict):
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
            df_in_boundaries = utils.filter_df_based_on_tuple_list(df=covariance_pmd_df,
                                                                   boundaries_list=chr_window_data)
            cov_dict[chromosome].append((patient, df_in_boundaries))

    return cov_dict


def main():
    args = parse_input()

    bedgraphs_path = files_tools.get_files_to_work(args.bedgraph_files, pattern=BEDGRPH_FILE_FORMAT)
    all_bedgraphs = files_tools.convert_paths_list_to_chromosome_based_dict(bedgraphs_path)
    boundaries_dict = files_tools.load_compressed_pickle(args.windows_file)

    cov_dict = get_bedgraph_to_df_per_chromosome_filtered(all_bedgraphs, boundaries_dict)

    chromosomes_list = []
    for chromosome in cov_dict:
        cpg_locations = utils.get_all_indexes_from_patients_list(cov_dict[chromosome])

        # Create df with location, chromosome and sequence
        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = cpg_locations
        chromosome_df["chromosome"] = chromosome
        chromosome_df["sequence"] = sequence_tools.get_sequences_for_cpgs(cpg_locations, chromosome)

        chromosome_df = chromosome_df.set_index("location")

        nc_meth_avg = utils.get_nc_avg(chromosome, chromosome_df.index, args.nc_files)
        chromosome_df.loc[nc_meth_avg.index, "nc_avg"] = nc_meth_avg

        for df_tuple in cov_dict[chromosome]:
            patient, patient_df = df_tuple
            patient_num = patient[-2:]
            cpgs_indexes = patient_df.index

            # This is covariance although we wrote covariance
            chromosome_df.loc[cpgs_indexes, "cov%s" % patient_num] = patient_df["coverage"]
            chromosome_df.loc[cpgs_indexes, "pmd_index"] = patient_df["pmd_index"]

            cancer_meth = utils.cpg_meth_in_cells(patient, chromosome, cpgs_indexes,
                                                  args.methylation_folder,
                                                  sublineage_name=utils.NOT_NC)
            chromosome_df.loc[cancer_meth.index, "meth%s" % patient_num] = cancer_meth.mean()
            chromosome_df.loc[cancer_meth.index, "var%s" % patient[-2:]] = cancer_meth.var()

            original_meth = utils.cpg_meth_in_cells(patient, chromosome, cpgs_indexes,
                                                    args.methylation_folder,
                                                    sublineage_name=utils.ONLY_NC)
            chromosome_df.loc[original_meth.index, "orig_meth%s" % patient_num] = original_meth.mean()


        chromosomes_list.append(chromosome_df.reset_index())

    all_chromosome_df = pd.concat(chromosomes_list)
    all_chromosome_df.to_pickle(os.path.join(args.output_folder, "valid_cpg.pkl"))


if __name__ == '__main__':
    main()
