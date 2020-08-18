import argparse
import os
import sys
import warnings
from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('seaborn')

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import files_tools
import commons.consts as consts


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_dict', help='Path to methylation files', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def convert_sublineage(sublineage_info, patient, cell):
    """
    :param sublineage_info: Dictionary of cell names -> sublineage
    :param patient: Patient name
    :param cell: The cell name
    :return: The cells sublineage if there is, undefined otherwise
    """
    try:
        return sublineage_info[patient + '_' + cell]
    except:
        return "no_info"


def chosen_cell_stats(sublineage, all_dict, output_path):
    """
    Calculate different stats for the low and high cells
    :param sublineage: Dictionary of cell names -> sublineage
    :param all_dict: Dictionary with lists of low and high cells per patient
    :param output_path: output folder path
    """
    region_high_counters = {
        "CRC01": Counter(LN1=0, LN2=0, LN3=0, ML1=0, ML2=0, ML3=0, ML4=0, MP1=0, MP2=0, MP3=0, MP4=0, MP5=0,
                         PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC04": Counter(LN1=0, LN2=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC10": Counter(LN1=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC11": Counter(LN1=0, PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC13": Counter(LN1=0, LN2=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC14": Counter(LN1=0, LN2=0, PT1=0, PT3=0, PT4=0, PT5=0)}
    region_low_counters = {
        "CRC01": Counter(LN1=0, LN2=0, LN3=0, ML1=0, ML2=0, ML3=0, ML4=0, MP1=0, MP2=0, MP3=0, MP4=0, MP5=0,
                         PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC04": Counter(LN1=0, LN2=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC10": Counter(LN1=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC11": Counter(LN1=0, PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC13": Counter(LN1=0, LN2=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC14": Counter(LN1=0, LN2=0, PT1=0, PT3=0, PT4=0, PT5=0)}

    print("Number of cells in each of the low/high groups when using 20%:")
    for patient in all_dict:
        num_of_cells = len(all_dict[patient]["high"])
        print(patient, ":", num_of_cells)

        region_stats(all_dict, num_of_cells, output_path, patient, region_high_counters, region_low_counters)

        sublineage_stats(all_dict, num_of_cells, output_path, patient, sublineage)


def sublineage_stats(all_dict, num_of_cells, output_path, patient, sublineage):
    """
    Calculate sublineage stats, how many of the low/high cells are in every sublineage
    :param all_dict: Dictionary with lists of low and high cells per patient
    :param num_of_cells: number of low/high cells
    :param output_path: output folder path
    :param patient: Patient name
    :param sublineage: Dictionary of cell names -> sublineage
    """
    high_sampling_region = [convert_sublineage(sublineage, patient, cell) for cell in
                            all_dict[patient]["high"]]
    low_sampling_region = [convert_sublineage(sublineage, patient, cell) for cell in
                           all_dict[patient]["low"]]
    sublineage_high_counter = Counter(high_sampling_region)
    sublineage_low_counters = Counter(low_sampling_region)

    all_cells = set(sublineage_high_counter.keys()) | set(sublineage_low_counters.keys())
    df = pd.DataFrame(index=all_cells).sort_index()
    high_df = pd.DataFrame.from_dict(sublineage_high_counter, orient='index', columns=["high"])
    low_df = pd.DataFrame.from_dict(sublineage_low_counters, orient='index', columns=["low"])
    df.loc[high_df.index, "high"] = high_df.high
    df.loc[low_df.index, "low"] = low_df.low

    df.plot.bar()
    plt.title("%s sublineage (%d cells)" % (patient, num_of_cells))
    plt.savefig(os.path.join(output_path, "sublineage_%s_cells" % patient))
    plt.close()


def region_stats(all_dict, num_of_cells, output_path, patient, region_high_counters, region_low_counters):
    """
    Calculate region stats, how many of the low/high cells are in every region
    :param all_dict: Dictionary with lists of low and high cells per patient
    :param num_of_cells: number of low/high cells
    :param output_path: output folder path
    :param patient: Patient name
    :param region_high_counters: Counters for the high cells
    :param region_low_counters: Counters for the low cells
    """
    high_sampling_region = [cell.split("_")[0] for cell in all_dict[patient]["high"]]
    low_sampling_region = [cell.split("_")[0] for cell in all_dict[patient]["low"]]
    region_high_counters[patient].update(high_sampling_region)
    region_low_counters[patient].update(low_sampling_region)
    df = pd.DataFrame.from_dict(region_high_counters[patient], orient='index', columns=["high"]).join(
        pd.DataFrame.from_dict(region_low_counters[patient], orient='index', columns=["low"]),
        lsuffix='_caller',
        rsuffix='_other')
    df.plot.bar()
    plt.title("%s cells sampling region (%d cells)" % (patient, num_of_cells))
    plt.savefig(os.path.join(output_path, "sampling_region_%s_cells" % patient))
    plt.close()


def main():
    args = parse_input()

    sublineage = files_tools.load_compressed_pickle(consts.CONVERT_SUBLINEAGE_LIOR)
    cell_dict = files_tools.load_compressed_pickle(args.cell_dict)
    output_path = args.output_folder

    chosen_cell_stats(sublineage, cell_dict, output_path)


if __name__ == '__main__':
    main()
