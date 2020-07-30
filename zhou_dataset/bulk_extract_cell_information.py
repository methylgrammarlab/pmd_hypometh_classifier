"""
Create a list of cells to use for variance and methylation and variance info for them for future parsing
based on extreme cells
"""

import argparse
import os
import pickle
import sys
import warnings

import numpy as np
import tqdm

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from variance import bulk_get_valid_cpgs_dataset
from commons import files_tools, consts


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--variance_meth_cells', help='Path to file containing the cells to calculate the '
                                                      'average methylation and variance', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--percentage_of_extreme_cells', required=False, type=int, default=15,
                        help='Percentage of extreme cells to take from top and bottom, default is 15% from '
                             'each side')
    args = parser.parse_args()
    return args


def main():
    args = parse_input()
    all_files = files_tools.get_files_to_work(args.methylation_folder, consts.BULK_FILE_FORMAT)

    variance_cells = bulk_get_valid_cpgs_dataset.get_cells_to_calculate_variance(
        args.variance_meth_cells, args.percentage_of_extreme_cells)

    cells_dict = {}
    # Go over all the data frames and get all the cpgs which match a cell and find the average methylation
    # and methylation variance for this cell
    for chromosome_file in tqdm.tqdm(all_files):
        chromosome, pmd_df = bulk_get_valid_cpgs_dataset.get_pmd_df(chromosome_file)
        df = pmd_df.filter(items=variance_cells, axis=0)

        for row in df.index.values:
            if row not in cells_dict:
                cells_dict[row] = []
            cells_dict[row].append(df.loc[row][np.isfinite(df.loc[row])].values.astype(np.float64))

    rows_values = [(row, np.mean(np.concatenate(cells_dict[row])), np.var(np.concatenate(cells_dict[row])))
                   for row in cells_dict]
    rows_values.sort(key=lambda x: x[1])

    with open(os.path.join(args.output_folder, "average_cells.pkl"), "wb") as f:
        pickle.dump(rows_values, f)


if __name__ == '__main__':
    main()
