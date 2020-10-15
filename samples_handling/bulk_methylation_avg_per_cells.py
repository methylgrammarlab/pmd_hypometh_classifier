import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm
import time

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools
import commons.consts as consts
from variance import bulk_get_valid_cpgs_dataset

SOLO = 0


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--orig_meth_cells', help='path for the file with original methylation cell information ',
                        required=False)
    parser.add_argument('--filter_nc', help='filter only nc cpg', required=False, type=bool, default=False)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def main():
    args = parse_input()

    all_files = files_tools.get_files_to_work(args.methylation_folder, consts.BULK_FILE_FORMAT % "*")

    cells_info_data = pd.read_csv(args.orig_meth_cells)
    orig_meth_cells = bulk_get_valid_cpgs_dataset.get_cells_to_calculate_original_meth(args.orig_meth_cells,
                                                                                       cells_info_data=cells_info_data)

    sum_methylation_df = pd.DataFrame()  # index=list(cells_info_data.loc[:, "sample"]))
    num_methylation_df = pd.DataFrame()  # index=list(cells_info_data.loc[:, "sample"]))

    for chromosome_file in tqdm(all_files):
        chromosome, pmd_df = bulk_get_valid_cpgs_dataset.get_pmd_df(chromosome_file)
        orig_meth_values = bulk_get_valid_cpgs_dataset.get_cpgs_orig_methylated(pmd_df, orig_meth_cells)
        nc_meth = orig_meth_values[orig_meth_values > 0.5].index

        df = pmd_df.astype(np.float64)
        if args.filter_nc:
            df = df.loc[:, df.columns & nc_meth]

        sum_methylation_df.loc[:, chromosome] = np.sum(df, axis=1)
        num_methylation_df.loc[:, chromosome] = np.sum(np.isfinite(df), axis=1)

    total_mean = np.sum(sum_methylation_df, axis=1) / np.sum(num_methylation_df, axis=1)
    combined = pd.DataFrame(index=total_mean.index)
    combined.loc[:, "mean"] = total_mean
    combined.loc[:, "nc_cells"] = cells_info_data.set_index("sample").loc[:, "nc_cells"]

    base_path = os.path.join(args.output_folder, "bulk_avg_data_all_")
    if args.filter_nc:
        base_path += "NC_"

    combined.to_csv(base_path + "mean_coverage.csv")


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    t0 = time.time()
    main()
    t1 = time.time()
    print("total time:", t1 - t0)
