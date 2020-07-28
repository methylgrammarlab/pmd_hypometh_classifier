import argparse
import glob
import os
import pickle
import re
import sys

import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from commons import consts
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
TOP_BOTTOM = 15

try:
    import pyfaidx

    PYFAIDX = True

except:
    PYFAIDX = False

if PYFAIDX:
    genome = pyfaidx.Fasta(consts.GENOME_FILE_LOCAL_DROR, sequence_always_upper=True, as_raw=True)


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--cell_info', help='Path to list containing the nc cells and variance cells',
                        required=True)
    parser.add_argument('--cells_avg_meth', help='Path to list containing the cells to use avg methylation',
                        required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_seq_for_cpg(chr_num, i, seq_size):
    chr_info = genome[chr_num]
    # assert chr_info[i] == "G"
    seq_size = int((seq_size - 2) / 2)
    return chr_info[i - seq_size - 1:i + seq_size + 1]


def get_seq_info(ind, chromosome):
    if not chromosome.isdigit():
        chromosome = INT_RE.findall(chromosome)[0]

    seq = []
    for i in ind:
        seq.append(get_seq_for_cpg(chromosome, i, SEQ_SIZE))

    return seq


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


def get_variance_cells_to_use(file_path):
    with open(file_path, "rb") as f:
        temp_variance_cells = pickle.load(f)

    temp_variance_cells.sort(key=lambda x: x[1])
    low = TOP_BOTTOM * len(temp_variance_cells) / 100 + 1
    top = len(temp_variance_cells) - TOP_BOTTOM * len(temp_variance_cells) / 100 - 1
    variance_cells = temp_variance_cells[:int(low)] + temp_variance_cells[int(top):]
    variance_cells = [i[0] for i in variance_cells]
    return variance_cells


def main():
    args = parse_input()

    all_files = glob.glob(os.path.join(args.methylation_folder, METHYLATION_FILE_FORMAT))
    cells_info_data = pd.read_csv(args.cell_info)
    nc_cells = list(cells_info_data[cells_info_data[NC_CELLS_COLUMN] == 1][CELL_NAME_COLUMN])
    nc_cells = [i.strip() for i in nc_cells]
    variance_cells = get_variance_cells_to_use(args.cells_avg_meth)
    pmd_data = handle_pmds.read_pmd_dict()

    df_list = []
    for chromosome_file in tqdm.tqdm(all_files):
        chromosome, pmd_df = get_pmd_df(chromosome_file)
        cpgs_methylated = get_cpgs_orig_methylated(pmd_df, nc_cells)
        df = pmd_df.filter(items=variance_cells, axis=0)
        df = df.loc[:, cpgs_methylated]  # Leave only methylated cells

        df = handle_pmds.remove_low_high_coverage(df)  # remove the top and low 5 percentage
        # Note: we didn't removed CpG with less then 5 samples in low or high like we did in the scWGBS but
        #  we checked that we don't have CpG like this

        chromosome_df = pd.DataFrame(columns=["chromosome", "location"])
        chromosome_df["location"] = df.columns
        chromosome_df["chromosome"] = chromosome
        chromosome_df = chromosome_df.set_index("location")

        chromosome_df = handle_pmds.get_pmd_index_based_on_tuples_list(chromosome_df, pmd_data[chromosome])
        chromosome_df["meth"] = df.mean()
        chromosome_df["var"] = df.var()

        if PYFAIDX:
            chromosome_df["sequence"] = get_seq_info(df.columns, str(chromosome))
        df_list.append(chromosome_df.reset_index())
        # Now we need to create the methylation and the variance + the sequence itself and the pmd index

    final_df = pd.concat(df_list)
    try:
        final_df.to_pickle(os.path.join(args.output_folder, "valid_cpg_zhou.pkl"))
    except:
        final_df.to_pickle("valid_cpg_zhou.pkl")


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    main()
