import seaborn as sns; sns.set(color_codes=True)
import argparse
import glob
import os
import pickle
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import commons.slurm_tools

# sys.path.append(os.path.dirname(os.getcwd()))

CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"  # TODO remove the mini
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")

def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--pmd_boundaries', help='Path to location of the file with the PMD boundries', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    parser.add_argument('--chr', help='Chromosome, all if not provided. e.g. chr16', required=False)
    args = parser.parse_args()
    return args


def create_cov_dend(file, pmd_boundries_list):
    df = pd.read_pickle(file)
    for boundaries in pmd_boundries_list:
        sliced = df.iloc[:, :200]
        cov = sliced.cov()
        ax = sns.heatmap(cov)
        plt.savefig('hmmm')


        # pmd = sns.clustermap(sliced, method=pd.DataFrame.cov)
        # pass


def main():
    args = parse_input()

    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])

    chr = args.chr

    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        if chr:
            cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % chr)

        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)
    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    # TODO ask dror what format we want to save the PMD boundries as...
    pmd_boundries_dict = {"chr16": [(6500000, 6500100)]}

    # for file in tqdm(all_cpg_format_file_paths, desc='files'):
    for file in all_cpg_format_file_paths:
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file)[0]
        create_cov_dend(file, pmd_boundries_dict[chromosome])


if __name__ == '__main__':
    # commons.slurm_tools.init_slurm(main)
    main()