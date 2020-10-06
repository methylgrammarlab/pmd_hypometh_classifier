import argparse
import os
import sys

import pandas as pd

from commons import files_tools, sequence_tools, consts


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder', help='Path to the csv files', required=False, default=None)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    args = parser.parse_args()
    return args


def extract_location_sequence(output_folder):
    """
    Extract the location and sequence of each chromosome except blacksites and CpGI and save it to the output
    :return:
    """
    context_map = files_tools.get_cpg_context_map(
        load_with_path=consts. / vol / sci / bio / data / benjamin.berman / bermanb / projects / genomic - data - misc / cgcontext / cgcontext.hg19,
        only_locations=True)

    for chromosome in context_map:
        sequences = sequence_tools.get_sequences_for_cpgs(indexes_list=context_map[chromosome],
                                                          chromosome=chromosome)
        df = pd.DataFrame()
        df["location"] = context_map[chromosome]
        df["sequences"] = sequences
        df.to_csv(path_or_buf=os.path.join(output_folder, "%s.csv" % chromosome))


def create_bedgraph(input_folder, output_folder):
    chromosome_files = files_tools.get_files_to_work(input_folder, "*.csv")
    df_list = []
    for csv_path in chromosome_files:
        chromosome = consts.CHR_FULL_NAME_RE.findall(os.path.basename(csv_path))[0]
        df = pd.read_csv(csv_path)
        df["chromosome"] = chromosome
        df_list.append(df)

    final_df = pd.concat(df_list)
    final_df["end"] = final_df["location"] + 1
    final_df.to_csv(path_or_buf=os.path.join(output_folder, "prediction.bedgraph"), header=False,
                    index=False, columns=["chromosome", "location", "end", "prediction"], sep=" ")

    final_df_gt09 = final_df[final_df["prediction"] >= 0.9]
    final_df_lt01 = final_df[final_df["prediction"] <= 0.1]
    final_df_lt01.to_csv(path_or_buf=os.path.join(output_folder, "prediction_lt01.bedgraph"), header=False,
                         index=False, columns=["chromosome", "location", "end", "prediction"], sep=" ")
    final_df_gt09.to_csv(path_or_buf=os.path.join(output_folder, "prediction_gt09.bedgraph"), header=False,
                         index=False, columns=["chromosome", "location", "end", "prediction"], sep=" ")


def main():
    args = parse_input()
    input_folder, output_folder = args.input_folder, args.output_folder

    if args.input_folder:
        create_bedgraph(input_folder, output_folder)
    else:
        extract_location_sequence(output_folder)


if __name__ == '__main__':
    main()
