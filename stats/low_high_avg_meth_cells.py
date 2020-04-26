import argparse
import glob
import os
import re
import sys
import warnings
from tqdm import tqdm
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('seaborn')

warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from format_files import handle_pmds
from commons import files_tools
from covariance import covariance_to_bedgraph

CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

BEDGRPH_FILE_FORMAT = os.path.join("*", "*.bedgraph")
BEDGRPAH_FORMAT_FILE_RE = re.compile(".*(CRC\d+)_chr(\d+).*")

METHYLATION_FILE_FORMAT = "all_cpg_ratios_%s_chr%s.dummy.pkl.zip"
SEQ_SIZE = 150
TOP_LOW_PERCENTAGE_TO_REMOVE = 5


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--methylation_folder', help='Path to methylation files', required=True)
    parser.add_argument('--windows_file', help='Path to files with windows we want to take', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def get_bedgraph_files(files):
    if os.path.isdir(files):
        file_path = os.path.join(files, BEDGRPH_FILE_FORMAT)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [files]

    return all_file_paths


def get_patient_dict(all_file_paths):
    d = {}
    for file_path in all_file_paths:
        try:
            patient, chromosome = BEDGRPAH_FORMAT_FILE_RE.findall(file_path)[0]
        except:
            continue

        if patient not in d:
            d[patient] = []

        d[patient].append((chromosome, file_path))

    return d


def get_cancer_methylation_of_patient(patient, chromosome, indexes, methylation_folder):
    methylation_file_path = os.path.join(methylation_folder, patient,
                                         METHYLATION_FILE_FORMAT % (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    _, df = covariance_to_bedgraph.get_region_df(df, sublineage_cells=[],
                                                 sublineage_name=covariance_to_bedgraph.ONLY_PT)
    mdf = df.mean()
    return mdf.loc[indexes].values


def main():
    args = parse_input()
    coverage_graph = True
    get_high_low = False
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    methylation_folder = args.methylation_folder
    patients = ["CRC01", "CRC04", "CRC10", "CRC11", "CRC13", "CRC14"]
    path = os.path.join(methylation_folder, "%s", "all*.pkl.zip")
    all_files_dict = get_patient_dict(
        [path for sublist in [glob.glob(path % patient) for patient in patients] for path in sublist])
    # all_files_dict = get_patient_dict(glob.glob(os.path.join(methylation_folder, "*", "*.pkl.zip")))
    global_windows_data = files_tools.load_compressed_pickle(args.windows_file)
    patients_dict = handle_pmds.get_cancer_pmd_df_with_windows_after_cov_filter(all_files_dict,
                                                                                global_windows_data)
    bins_patients_dict = {}
    for patient in tqdm(patients_dict):
        all_means = []
        tot_num_of_cpg = []
        general_num_of_cpg = 0
        for chr in patients_dict[patient]:
            df = patients_dict[patient][chr]
            mean = df.mean(axis=1)
            num_of_cpg = (~df.isnull()).sum(axis=1)
            all_means.append(mean.mul(num_of_cpg, axis=0))
            tot_num_of_cpg.append(num_of_cpg)
            general_num_of_cpg += df.shape[1]
        sum_of_means = pd.concat(all_means, axis=1).sum(axis=1)
        div_num = pd.concat(tot_num_of_cpg, axis=1).sum(axis=1)
        coverage = div_num / general_num_of_cpg
        avg = sum_of_means.div(div_num, axis=0)

        if (coverage_graph):
            save_coverage_graphs(avg, coverage, output, patient)

        if (get_high_low):
            cell_mean = avg.sort_values()
            bottom = min(len(cell_mean) - 1, int((len(cell_mean) / 100) * 20))
            top = int((len(cell_mean) / 100) * 20)
            low_15 = list(cell_mean.iloc[:bottom].index)
            high_15 = list(cell_mean.iloc[-top - 1:-1].index)
            bins_patients_dict[patient] = {"low": low_15, "high": high_15}

    if get_high_low:
        output_path = os.path.join(output, "20perc_low_high_avg_methylation_cells.pickle.zlib")
        files_tools.save_as_compressed_pickle(output_path, bins_patients_dict)


def save_coverage_graphs(avg, coverage, output, patient):
    plt.scatter(coverage, avg)
    plt.xlabel("coverage (percentage)")
    plt.ylabel("average methylation")
    plt.title("%s Average vs Covariance" % patient)
    plt.savefig(os.path.join(output, "%s_coverage.png" % patient))
    plt.close()


def convert_sublineage(sublineage_info, patient, cell):
    try:
        return sublineage_info[patient + '_' + cell]
    except:
        return "no_info"


def chosen_cell_stats(path):
    sublineage = files_tools.load_compressed_pickle(
        R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\stats\top_bottom\convert_sublineage.pickle.zlib")
    region_high_counters = {
        "CRC01": Counter(LN1=0, LN2=0, LN3=0, ML1=0, ML2=0, ML3=0, ML4=0, MP1=0, MP2=0, MP3=0, MP4=0, MP5=0, PT1=0,
                         PT2=0, PT3=0, PT4=0), "CRC04": Counter(LN1=0, LN2=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC10": Counter(LN1=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0), "CRC11": Counter(LN1=0, PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC13": Counter(LN1=0, LN2=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC14": Counter(LN1=0, LN2=0, PT1=0, PT3=0, PT4=0, PT5=0)}
    region_low_counters = {
        "CRC01": Counter(LN1=0, LN2=0, LN3=0, ML1=0, ML2=0, ML3=0, ML4=0, MP1=0, MP2=0, MP3=0, MP4=0, MP5=0, PT1=0,
                         PT2=0, PT3=0, PT4=0), "CRC04": Counter(LN1=0, LN2=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC10": Counter(LN1=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0), "CRC11": Counter(LN1=0, PT1=0, PT2=0, PT3=0, PT4=0),
        "CRC13": Counter(LN1=0, LN2=0, LN3=0, PT1=0, PT2=0, PT3=0, PT4=0, PT5=0),
        "CRC14": Counter(LN1=0, LN2=0, PT1=0, PT3=0, PT4=0, PT5=0)}
    all_dict = files_tools.load_compressed_pickle(path)
    print("Number of cells in each of the low/high groups when using 20%:")
    for patient in all_dict:
        num_of_cells = len(all_dict[patient]["high"])
        print(patient, ":", num_of_cells)

        high_sampling_region = [cell.split("_")[0] for cell in all_dict[patient]["high"]]
        low_sampling_region = [cell.split("_")[0] for cell in all_dict[patient]["low"]]
        region_high_counters[patient].update(high_sampling_region)
        region_low_counters[patient].update(low_sampling_region)
        df = pd.DataFrame.from_dict(region_high_counters[patient], orient='index', columns=["high"]).join(
            pd.DataFrame.from_dict(region_low_counters[patient], orient='index', columns=["low"]), lsuffix='_caller',
            rsuffix='_other')
        df.plot.bar()
        plt.title("%s cells sampling region (%d cells)" % (patient, num_of_cells))
        plt.savefig("top_bottom\sampling_region_%s_cells" % patient)
        plt.close()

        high_sampling_region = [convert_sublineage(sublineage, patient, cell) for cell in all_dict[patient]["high"]]
        low_sampling_region = [convert_sublineage(sublineage, patient, cell) for cell in all_dict[patient]["low"]]
        sublineage_high_counter = Counter(high_sampling_region)
        sublineage_low_counters = Counter(low_sampling_region)
        all_cells = set(sublineage_high_counter.keys()) | set(sublineage_low_counters.keys())
        df = pd.DataFrame(index=all_cells).sort_index()
        high_df = pd.DataFrame.from_dict(sublineage_high_counter, orient='index', columns=["high"])
        low_df =pd.DataFrame.from_dict(sublineage_low_counters, orient='index', columns=["low"])
        df.loc[high_df.index, "high"] = high_df.high
        df.loc[low_df.index, "low"] = low_df.low
        df.plot.bar()
        plt.title("%s sublineage (%d cells)" % (patient, num_of_cells))
        plt.savefig("top_bottom\sublineage_%s_cells" % patient)
        plt.close()
        pass


if __name__ == '__main__':
    # slurm_tools.init_slurm(main)
    # main()
    chosen_cell_stats(
        R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\stats\top_bottom\20perc_low_high_avg_methylation_cells.pickle.zlib")

# To run -  python3 low_high_avg_meth_cells.py --methylation_folder /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/cpg_format/filtered_by_bl_and_cpgi/ --windows_file /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/covariance/cancer/5000_with_750_minimum_pairs/boundries/window_boundries.dummy.pkl.zip --output_folder /vol/sci/bio/data/benjamin.berman/bermanb/projects/scTrio-seq-reanalysis/liordror/stats/low_high_avg_meth_cells
