import argparse
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

from commons import files_tools, consts
from format_files import handle_pmds
from variance import bulk_get_valid_cpgs_dataset


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input', help='path for the input file or folder', required=True)
    parser.add_argument('--orig_meth_cells', help='path for the file with original methylation cell '
                                                  'information ', required=False)
    parser.add_argument('--parse_format', help='do we parse scwgb or bulk data', required=True, type=str,
                        choices=consts.PARSE_FORMAT_OPTIONS)
    args = parser.parse_args()
    return args


def plot_scwgbs_average_methylation(input_path, output_folder):
    """
    Plot the scwgbs average methyatlionn histogram
    :param input_path: Input folder with the average methylation information
    :param output_folder: Output folder for the images
    """
    pmd_list = []
    non_pmd_list = []
    for f in files_tools.get_files_to_work(input_path, "*.dummy.pkl.zip"):
        all_data = pd.read_pickle(f).mean(axis=1)
        all_data = all_data[~pd.isnull(all_data)]

        chromosome = consts.CHR_FULL_NAME_RE.findall(os.path.basename(f))[0]
        pmd_df = handle_pmds.filtered_out_non_pmd(df=all_data, chromosome=chromosome)
        pmd_list.append(pmd_df)
        non_pmd_list.append(all_data.drop(index=pmd_df.index))

    pmd_df = pd.concat(pmd_list)
    _ = pmd_df.hist(bins=20)
    plt.title("Histogram of average methylation in normal cells in PMD scWGBS")
    plt.xlabel("Avg methylation")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_folder, "average_methylation_normal_cells_pmd_scwgbs.png"))

    non_pmd_df = pd.concat(non_pmd_list)
    _ = non_pmd_df.hist(bins=20)
    plt.title("Histogram of average methylation in normal cells in Non-PMD scWGBS")
    plt.xlabel("Avg methylation")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_folder, "average_methylation_normal_cells_non_pmd_scwgbs.png"))


def plot_bulk_average_methylation(input_path, output_folder, orig_meth_cells):
    all_files = files_tools.get_files_to_work(input_path, consts.BULK_FILE_FORMAT % "*")
    orig_meth_cells = bulk_get_valid_cpgs_dataset.get_cells_to_calculate_original_meth(orig_meth_cells)

    pmd_list = []
    non_pmd_list = []
    for chromosome_file in all_files:
        chromosome = consts.CHR_FULL_NAME_RE.findall(os.path.basename(chromosome_file))[0]
        all_data = pd.read_pickle(chromosome_file)
        pmd_df = handle_pmds.filtered_out_non_pmd(all_data, chromosome, pmd_file=consts.PMD_FILE_LOCAL_DROR)
        non_pmd_df = all_data.drop(index=pmd_df.index)

        orig_meth_pmd = bulk_get_valid_cpgs_dataset.get_cpgs_orig_methylated(pmd_df, orig_meth_cells)
        orig_meth_non_pmd = bulk_get_valid_cpgs_dataset.get_cpgs_orig_methylated(non_pmd_df, orig_meth_cells)

        pmd_list.append(orig_meth_pmd)
        non_pmd_list.append(orig_meth_non_pmd)

    pmd_df = pd.concat(pmd_list)
    _ = pmd_df.hist(bins=20)
    plt.title("Histogram of average methylation in normal cells in PMD Zhou")
    plt.xlabel("Avg methylation")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_folder, "average_methylation_normal_cells_pmd_zhou.png"))

    non_pmd_df = pd.concat(non_pmd_list)
    _ = non_pmd_df.hist(bins=20)
    plt.title("Histogram of average methylation in normal cells in Non-PMD Zhou")
    plt.xlabel("Avg methylation")
    plt.ylabel("Count")
    plt.savefig(os.path.join(output_folder, "average_methylation_normal_cells_non_pmd_zhou.png"))


def main():
    args = parse_input()
    input_path, output_folder, parse_format = args.input, args.output_folder, args.parse_format

    if parse_format == consts.SCWGBS:
        plot_scwgbs_average_methylation(input_path, output_folder)
    else:
        plot_bulk_average_methylation(input_path, output_folder, args.orig_meth_cells)


if __name__ == '__main__':
    main()
