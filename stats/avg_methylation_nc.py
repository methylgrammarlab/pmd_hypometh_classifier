#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python
import argparse
import collections
import copy
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import data_tools, files_tools

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_chr%d.dummy.pkl.zip"
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

MET_AVG_FILE_FORMAT = "nc_average_methylation_chr%d.dummy.pkl.zip"

BEDGRAPH_OUTPUT_FILE_FORMAT = "nc_methylated_coverage_chr%d_threshold_%d.bedgraph"
BEDGRAPH_FILE_FORMAT = "nc_methylated_coverage_chr*_threshold_*.bedgraph"
BEDGRAPH_FILE_FORMAT_RE = re.compile(".+nc_methylated_coverage_chr(\d+)_threshold_\d+.bedgraph")
CSV_FILE = "average_methylation_of_nc_%s.csv"
PATIENTS = ['CRC01', 'CRC13', 'CRC11']
# PATIENTS = ['CRC01', 'CRC13', 'CRC04', 'CRC10', 'CRC11']

MET_THRESHOLD = 0.5


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_folder', help='Path to folder of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--nc_avg',
                        help='Create a csv file with counts of how many times each methylation average appears',
                        required=False)
    parser.add_argument('--methylation_diff', help='', required=False)
    parser.add_argument('--nc_coverage_bedgraph',
                        help='Counts how many of the patients cover each CpG and output to a bedgraph', required=False)
    parser.add_argument('--nc_coverage_inds',
                        help='Creates a pickle file with a dictionary of the methylated indices per chromosome',
                        required=False)
    parser.add_argument('--methylation_coverage',
                        help='Receives the above bedgraph and creates a csv of how many CpGs are methylated for how many patients',
                        required=False)
    parser.add_argument('--methylation_average',
                        help='Creates a DF per chromosome t=with the avergae methylation of each patient',
                        required=False)

    args = parser.parse_args()
    return args


def average_nc_methylation(cpg_format_file):
    """
    Creates a list of all the average values of the normal cells per patient.
    :param cpg_format_file: File with the CpG data
    :return: List containing all the average values
    """
    df = pd.read_pickle(cpg_format_file)
    normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    normal_df = df.loc[normal_cell_ids, :]
    average = np.mean(normal_df, axis=0)
    return average.dropna().values


def avg_main(args, output):
    """
    Create a csv file with counts of how many times each methylation average appears
    :param args: The program arguments
    :param output: The output folder
    """
    cpg_format_file_path = os.path.join(args.cpg_format_folder, CPG_FORMAT_FILE_FORMAT)
    all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    counter = collections.Counter()

    patient = os.path.split(args.cpg_format_folder)[-1]

    # for file in tqdm(all_cpg_format_file_paths, desc='files'):
    for file in tqdm(all_cpg_format_file_paths):
        average = average_nc_methylation(file)
        counter.update(average)
    data_tools.counter_to_csv(counter, os.path.join(output, CSV_FILE % patient))


def create_nc_coverage_bedgraph(patients_list, chr, output, bedgraph=False, dump_indices=False):
    """
    creates different stats of the coverage of normal cells - fins
    If dump_indices saves a dictionary of all the methylated indices - over MET_THRESHOLD% of the cells had at least a
    threshold% methylation average and if bedgraph saves to a bedgraph the number per CpG.
    :param patients_list: All the paths of CpG files per chromosome
    :param chr: The current chromosome
    :param output: The output folder
    :param bedgraph: Save to bedgraph?
    :param dump_indices: Return a list of the indices?
    :return: Only if dump_indices
    """
    all_met = []
    all_nan = []
    threshold = 0.6
    path = os.path.join(output, BEDGRAPH_OUTPUT_FILE_FORMAT % (chr, threshold * 10))

    not_nan = None
    for patient in patients_list:
        df = pd.read_pickle(patient)
        normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        normal_df = df.loc[normal_cell_ids, :]
        average = np.mean(normal_df, axis=0)
        met = average >= threshold
        nan = ~np.isnan(average)
        if not_nan is None:
            not_nan = average.index[np.where(average.notnull())[0]]
        else:
            not_nan &= average.index[np.where(average.notnull())[0]]
        all_met.append(met)
        all_nan.append(nan)

    if len(all_met) == 0:
        return None

    all_met_df = pd.concat(all_met, axis=1).astype(int)
    met_coverage = np.sum(all_met_df, axis=1)

    if bedgraph:
        met_coverage_not_nan = met_coverage.loc[not_nan]
        index = pd.Series(met_coverage_not_nan.index).astype(int)
        chromosome = pd.DataFrame(np.full((index.shape[0],), chr))
        bed = pd.concat([chromosome, index, index + 1, met_coverage_not_nan.reset_index(drop=True)], axis=1)
        bed.to_csv(path, sep='\t', header=False, index=False)

    if dump_indices:
        all_nan_df = pd.concat(all_nan, axis=1).astype(int)
        nan_coverage = np.sum(all_nan_df, axis=1)
        met_nan = pd.concat([met_coverage, nan_coverage], axis=1, names=['met', 'nan'])
        met_nan_ratio = met_coverage / nan_coverage
        methylated_indices = nan_coverage.index[met_nan_ratio >= MET_THRESHOLD]
        return list(methylated_indices)


def methylation_diff(patients_list):
    """
    Finds the percentage per chromosome of the number of CpGs that were methylated (average of over 0.6%) out of the
    total amount of covered CpGs.
    :param patients_list:
    :return:
    """
    all_met_ind = []

    not_nan = None
    for patient in patients_list:
        df = pd.read_pickle(patient)
        normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        normal_df = df.loc[normal_cell_ids, :]
        average = np.mean(normal_df, axis=0)
        met_ind = average.index[np.where(average >= 0.6)[0]]
        if not_nan is None:
            not_nan = average.index[np.where(average.notnull())[0]]
        else:
            not_nan &= average.index[np.where(average.notnull())[0]]
        all_met_ind.append(met_ind)

    if len(all_met_ind) == 0:
        return None

    met_and = copy.deepcopy(all_met_ind[0])
    met_or = copy.deepcopy(all_met_ind[0])
    for i in range(1, len(all_met_ind)):
        met_ind = copy.deepcopy(all_met_ind[i])
        met_and &= met_ind
        met_or |= met_ind

    return len(met_and & not_nan) / len(met_or & not_nan) * 100


def diff_main(args, output):
    """
    Different stats depending on the args
    :param output:
    :return:
    """
    f = open(os.path.join(output, "avg_methylation_60.out"), "w")
    all_cpg_format_file_paths = create_chr_paths_dict(args)

    for chr in tqdm(all_cpg_format_file_paths):
        v = methylation_diff(all_cpg_format_file_paths[chr])
        f.write("chr%s:%s\n" % (chr, v))
        print(v)


def nc_coverage_main(args, output):
    """
    Different stats depending on the args
    :param output:
    :return:
    """
    all_cpg_format_file_paths = create_chr_paths_dict(args)

    methylation_indices_dict = {}

    for chr in tqdm(all_cpg_format_file_paths):
        indices_list = create_nc_coverage_bedgraph(all_cpg_format_file_paths[chr], chr, output,
                                                   bedgraph=args.nc_coverage_bedgraph,
                                                   dump_indices=args.nc_coverage_inds)
        methylation_indices_dict[chr] = indices_list
    path = os.path.join(output, "methylated_indices_threshold_%d.pickle.zlib" % (MET_THRESHOLD * 100))
    files_tools.save_as_compressed_pickle(path, methylation_indices_dict)


def create_chr_paths_dict(args):
    all_cpg_format_file_paths = dict()
    for chr in range(1, 23):
        all_cpg_format_file_paths[chr] = []
    for patient in PATIENTS:
        for chr in range(1, 23):
            cpg_format_file_path = os.path.join(args.cpg_format_folder, patient, CPG_FORMAT_FILE_FORMAT % chr)
            all_cpg_format_file_paths[chr] += glob.glob(cpg_format_file_path)
    return all_cpg_format_file_paths


def met_coverage_main(args, output):
    """
    Counts how many of the patients cover each CpG and output to a bedgraph
    :param args: The program arguments
    :param output: The output folder
    """
    cpg_format_file_path = os.path.join(args.cpg_format_folder, BEDGRAPH_FILE_FORMAT)
    all_methylation_coverage_paths = glob.glob(cpg_format_file_path)

    f = open(os.path.join(output, "methylation_coverage_60.csv"), "w")
    f.write("chr, 1, 2, 3, 4, 5\n")
    for path in all_methylation_coverage_paths:
        chr = BEDGRAPH_FILE_FORMAT_RE.findall(path)[0]
        df = pd.read_csv(path, sep='\t')
        agree_percent = []
        num_of_met = len(np.where(df.iloc[:, -1] > 0)[0])
        for i in range(5):  # num of patients
            num_over = len(np.where(df.iloc[:, -1] > i)[0])
            agree_percent.append(num_over / num_of_met * 100)
        f.write("chr%s," % chr)
        f.write(",".join(str(perc) for perc in agree_percent))
        f.write('\n')
    f.close()


def create_patient_series(patients_list):
    """
    Receives a list of all the paths of CpG files per chromosome, creates an average of all the normal cells and return
    as a DF of patient/genome location.
    """
    avg_dict = {}
    for patient in patients_list:
        name = CPG_FORMAT_FILE_RE.findall(patient)[0][0]
        df = pd.read_pickle(patient)
        normal_cell_ids = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
        normal_df = df.loc[normal_cell_ids, :]
        average = np.mean(normal_df, axis=0)
        avg_dict[name] = average
    return pd.DataFrame(avg_dict)


def avg_df_main(args, output):
    """
    Creates a DF per chromosome t=with the avergae methylation of each patient.
    :param args: The program arguments
    :param output: The output folder
    """
    all_cpg_format_file_paths = create_chr_paths_dict(args)

    for chr in tqdm(all_cpg_format_file_paths):
        path = os.path.join(output, MET_AVG_FILE_FORMAT % chr)
        chr_avgs = create_patient_series(all_cpg_format_file_paths[chr])
        chr_avgs.to_pickle(path)


def main():
    args = parse_input()
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    if args.nc_avg:
        avg_main(args, output)
    elif args.methylation_diff:
        diff_main(args, output)
    elif args.methylation_coverage:
        met_coverage_main(args, output)
    elif args.nc_coverage_inds or args.nc_coverage_bedgraph:
        nc_coverage_main(args, output)
    elif args.methylation_average:
        avg_df_main(args, output)


if __name__ == '__main__':
    main()
