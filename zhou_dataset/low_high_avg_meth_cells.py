import argparse
import glob
import os
import pickle
import re
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('seaborn')

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from format_files import handle_pmds
import tqdm

MEAN = 1
VAR = 2
METHYLATION_FILE_FORMAT = "all_cpg_ratios_chr*_hg19.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+_(chr\d+)_hg19.dummy.pkl.zip")
INT_RE = re.compile("chr(\d+)")

CELL_NAME_COLUMN = "sample"
NC_CELLS_COLUMN = "nc_cells"
VARIANCE_CELLS_COLUMN = "variance_cells"

SEQ_SIZE = 150


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--cell_info', help='Path to list containing the nc cells and variance cells',
                        required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_pmd_df(chromosome_file):
    """
    Get the data frame with information about the chromosome methylation only for pmd
    :param chromosome_file: The file with the information of the cpg in a specific chromosome
    :return: The data frame with information only for PMD
    :rtype: pd.DataFrame
    """
    chromosome = CPG_FORMAT_FILE_RE.findall(os.path.basename(chromosome_file))[0]
    df = pd.read_pickle(chromosome_file)
    return chromosome, handle_pmds.get_pmd_df(df, chromosome)


def get_cpgs_orig_methylated(df, cells_to_use):
    """
    Get a list of cpg which where methylated to begin with using a specific set of cells
    :param df: The df to work on
    :type df: pd.DataFrame
    :param cells_to_use: A list of cells to use, should match the rows names of the df
    :return: The indexes of CpG to use
    """
    df_orig = df.filter(items=cells_to_use, axis=0)
    mean_values = df_orig.mean()
    return mean_values >= 0.6


def plot_graph(avg, output):
    plt.bar([i for i in range(len(avg))], avg)
    plt.xlabel("Cells")
    plt.ylabel("Average methylation")
    plt.title("Average Methylation of cells")
    plt.savefig(os.path.join(output, "average_methylation_of_cells.png"))
    plt.close()
    # plt.show()


def main():
    args = parse_input()

    all_files = glob.glob(os.path.join(args.methylation_folder, METHYLATION_FILE_FORMAT))
    cells_info_data = pd.read_csv(args.cell_info)
    nc_cells = list(cells_info_data[cells_info_data[NC_CELLS_COLUMN] == 1][CELL_NAME_COLUMN])
    nc_cells = [i.strip() for i in nc_cells]
    variance_cells = list(cells_info_data[cells_info_data[VARIANCE_CELLS_COLUMN] == 1][
                              CELL_NAME_COLUMN])
    variance_cells = [i.strip() for i in variance_cells]

    rows_dict = {}
    for chromosome_file in tqdm.tqdm(all_files):
        chromosome, pmd_df = get_pmd_df(chromosome_file)
        pmd_df = pmd_df.fillna(0)  # Replace the nan with 0 for next calculation steps
        cpgs_methylated = get_cpgs_orig_methylated(pmd_df, nc_cells)
        df = pmd_df.filter(items=variance_cells, axis=0)
        df = df.loc[:, cpgs_methylated]  # Leave only methylated cells

        df = handle_pmds.remove_low_high_coverage(df)  # remove the top and low 5 percentage
        for row in df.index.values:
            if row not in rows_dict:
                rows_dict[row] = []
            rows_dict[row].append(df.loc[row].values)

    rows_values = [(row, np.mean(rows_dict[row])) for row in rows_dict]
    rows_values.sort(key=lambda x: x[1])
    # plot_graph([i[1] for i in rows_values], args.output_folder)
    with open(os.path.join(args.output_folder, "average_cells.pkl"), "wb") as f:
        pickle.dump(rows_values, f)


def plot_from_pickle():
    args = parse_input()
    with open(args.methylation_folder, "rb") as f:
        rows_values = pickle.load(f)

    plot_graph([i[1] for i in rows_values], args.output_folder)


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    # main()
    plot_from_pickle()
