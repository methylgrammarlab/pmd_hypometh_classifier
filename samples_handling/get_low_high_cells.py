import argparse
import os
import sys
import warnings

import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_stats', help='Path to methylation stats file', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def get_low_high_cells(methylation_df, patient, percentage):
    """
    For the given patient calculate the percentage high and low cells, excluding normal cells.
    :param methylation_df: df containing the average methylation levels for each cell.
    :param patient: The patient to calculate cells for.
    :param percentage: how many percent of cells
    :return: two lists, names of the top and bottom percentage percent of cells
    """
    df = methylation_df.loc[(methylation_df.patient == patient) & (methylation_df.lesion != "NC"), :]
    df = df.sort_values(by=['mean'])
    bottom = min(len(df) - 1, int((len(df) / 100) * percentage))
    top = int((len(df) / 100) * percentage)
    low_perc = list(df.iloc[:bottom].index)
    high_perc = list(df.iloc[::-1][:top].index)
    return high_perc, low_perc


def main():
    args = parse_input()

    methylation_df = pd.read_csv(args.methylation_stats, index_col=0)

    bins_patients_dict = {}
    percentage = 20

    for patient in methylation_df.patient.unique():
        high_perc, low_perc = get_low_high_cells(methylation_df, patient, percentage)
        bins_patients_dict[patient] = {"low": low_perc, "high": high_perc}

    output_path = os.path.join(args.output_folder, "%sperc_low_high_avg_methylation_cells.pickle.zlib" % percentage)
    files_tools.save_as_compressed_pickle(output_path, bins_patients_dict)


if __name__ == '__main__':
    main()
