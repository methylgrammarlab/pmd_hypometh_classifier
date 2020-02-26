import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from pandas_profiling import ProfileReport
from tqdm import tqdm

HIGH_T = 0.5
LOW_T = 0.25
MID = 0.5

sys.path.append(os.getcwd())
from format_files.format_cpg_context_map import NUMBER_OF_ORPH_PER_INDEX
from format_files import format_cpg_context_map, handle_pmds

ORPH_COLS = ["num_cpg_in_%s" % i for i in NUMBER_OF_ORPH_PER_INDEX]
CPG_FORMAT_FILE_RE = re.compile(".+(CRC\d+)_(chr\d+).dummy.pkl.zip")
CPG_FORMAT_FILE_FORMAT = "all_cpg_ratios_*_%s.dummy.pkl.zip"


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpg_format_files', help='Path to folder or file of parsed scWGBS', required=True)
    parser.add_argument('--output_folder', help='Path of the output folder', required=False)
    args = parser.parse_args()
    return args


def format_args():
    """
    Format the args for this script
    :return: The path of the files and the output directory
    """
    args = parse_input()
    output = args.output_folder
    if not output:
        output = os.path.dirname(sys.argv[0])
    if os.path.isdir(args.cpg_format_files):
        cpg_format_file_path = os.path.join(args.cpg_format_files, CPG_FORMAT_FILE_FORMAT % '*')
        all_cpg_format_file_paths = glob.glob(cpg_format_file_path)

    else:
        all_cpg_format_file_paths = [args.cpg_format_files]

    return all_cpg_format_file_paths, output


def collect_data_global(df, chr_info, patient, chromosome):
    context_3 = np.array(format_cpg_context_map.get_context_as_str_for_chr(chr_info))
    context_2 = np.array(format_cpg_context_map.get_context_as_str_for_chr_2(chr_info))
    context_1 = np.array(format_cpg_context_map.get_context_as_str_for_chr_1(chr_info))
    orph_info = chr_info[:, 1:14]
    weak = format_cpg_context_map.get_weak_column(chr_info)
    strong = format_cpg_context_map.get_strong_column(chr_info)
    weak_or_strong = np.logical_xor(weak, strong)
    solo_col = format_cpg_context_map.get_orph_35_column(chr_info)

    nc_index = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    pt_index = [cell_id for cell_id in df.index if cell_id.startswith('PT')]

    nc_values = df.loc[nc_index, :].mean(axis=0, skipna=True)
    pt_values = df.loc[pt_index, :].mean(axis=0, skipna=True)

    pt_values[pt_values >= HIGH_T] = 1
    pt_values[pt_values <= LOW_T] = 0
    nc_values[nc_values <= LOW_T] = 0
    nc_values[nc_values >= HIGH_T] = 1

    print("nc meth: %s, perc: %s , nc not-meth:%s\n" %
          (np.sum(nc_values == 1), np.sum(nc_values == 1) / nc_values.size, np.sum(nc_values == 0)))

    ### No change ###
    kept_meth = np.logical_and(pt_values == 1, nc_values == 1)
    print("No changes:%s" % (np.sum(kept_meth) / kept_meth.size))
    f_context_3 = context_3[np.where(kept_meth)]
    f_context_2 = context_2[np.where(kept_meth)]
    f_context_1 = context_1[np.where(kept_meth)]
    f_orph_info = orph_info[np.where(kept_meth)]
    f_weak = weak[np.where(kept_meth)]
    f_strong = strong[np.where(kept_meth)]
    f_weak_or_strong = weak_or_strong[np.where(kept_meth)]

    final_table = np.hstack((f_orph_info,
                             f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
                             f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))

    end_df = pd.DataFrame(final_table, columns=ORPH_COLS +
                                               ["weak", "strong", "weak_or_strong", "context_1", "context_2",
                                                "context_3"])
    profile = ProfileReport(end_df, title='No methylation change report',
                            html={'style': {'full_width': True}})
    profile.to_file(output_file=f"profile_no_met_change_{patient}_{chromosome}.html")

    ### Hypo ###
    hypo_meth = np.logical_and(pt_values == 0, nc_values == 1)
    print("hypo meth: %s" % (np.sum(hypo_meth) / hypo_meth.size))
    f_context_3 = context_3[np.where(hypo_meth)]
    f_context_2 = context_2[np.where(hypo_meth)]
    f_context_1 = context_1[np.where(hypo_meth)]
    f_orph_info = orph_info[np.where(hypo_meth)]
    f_weak = weak[np.where(hypo_meth)]
    f_strong = strong[np.where(hypo_meth)]
    f_weak_or_strong = weak_or_strong[np.where(hypo_meth)]

    final_table = np.hstack((f_orph_info,
                             f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
                             f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))

    end_df = pd.DataFrame(final_table,
                          columns=ORPH_COLS +
                                  ["weak", "strong", "weak_or_strong", "context_1", "context_2", "context_3"])
    profile = ProfileReport(end_df, title='hypo methylation change report',
                            html={'style': {'full_width': True}})
    profile.to_file(output_file=f"profile_hypo_{patient}_{chromosome}.html")

    ### Hyper ###
    # hyper_meth = np.logical_and(pt_values ==1 , nc_values ==0 )
    # print("hyper met: %s" % (np.sum(hyper_meth) / hypo_meth.size))
    # f_context_3 = context_3[np.where(hyper_meth)]
    # f_context_2 = context_2[np.where(hyper_meth)]
    # f_context_1 = context_1[np.where(hyper_meth)]
    # f_orph_info = orph_info[np.where(hyper_meth)]
    # f_weak = weak[np.where(hyper_meth)]
    # f_strong = strong[np.where(hyper_meth)]
    # f_weak_or_strong = weak_or_strong[np.where(hyper_meth)]
    #
    # final_table = np.hstack((f_orph_info,
    #                          f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
    #                          f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))
    #
    # end_df = pd.DataFrame(final_table,
    #                       columns=ORPH_COLS +
    #                               ["weak", "strong", "weak_or_strong", "context_1", "context_2", "context_3"])
    # profile = ProfileReport(end_df, title='hyper methylation change report',
    #                         html={'style': {'full_width': True}})
    # profile.to_file(output_file=f"profile_hyper_{patient}_{chromosome}.html")


def collect_data_not_solo(df, chr_info, patient, chromosome):
    context_3 = np.array(format_cpg_context_map.get_context_as_str_for_chr(chr_info))
    context_2 = np.array(format_cpg_context_map.get_context_as_str_for_chr_2(chr_info))
    context_1 = np.array(format_cpg_context_map.get_context_as_str_for_chr_1(chr_info))
    orph_info = chr_info[:, 1:14]
    weak = format_cpg_context_map.get_weak_column(chr_info)
    strong = format_cpg_context_map.get_strong_column(chr_info)
    weak_or_strong = np.logical_xor(weak, strong)
    solo_col = format_cpg_context_map.get_orph_35_column(chr_info)
    solo_col[solo_col > 0] = 1

    nc_index = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    pt_index = [cell_id for cell_id in df.index if cell_id.startswith('PT')]

    nc_values = df.loc[nc_index, :].mean(axis=0, skipna=True)
    pt_values = df.loc[pt_index, :].mean(axis=0, skipna=True)

    pt_values[pt_values >= HIGH_T] = 1
    pt_values[pt_values <= LOW_T] = 0
    nc_values[nc_values <= LOW_T] = 0
    nc_values[nc_values >= HIGH_T] = 1

    ### No change ###
    kept_meth = np.logical_and(pt_values == 1, nc_values == 1)
    kept_meth = np.logical_and(solo_col, kept_meth)
    print("No changes: %s" % (np.sum(kept_meth) / np.sum(solo_col)))
    f_context_3 = context_3[np.where(kept_meth)]
    f_context_2 = context_2[np.where(kept_meth)]
    f_context_1 = context_1[np.where(kept_meth)]
    f_orph_info = orph_info[np.where(kept_meth)]
    f_weak = weak[np.where(kept_meth)]
    f_strong = strong[np.where(kept_meth)]
    f_weak_or_strong = weak_or_strong[np.where(kept_meth)]

    final_table = np.hstack((f_orph_info,
                             f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
                             f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))

    end_df = pd.DataFrame(final_table, columns=ORPH_COLS +
                                               ["weak", "strong", "weak_or_strong", "context_1", "context_2",
                                                "context_3"])
    profile = ProfileReport(end_df, title='No methylation change report - not solo',
                            html={'style': {'full_width': True}})
    profile.to_file(output_file=f"profile_no_met_change_{patient}_{chromosome}_not_solo.html")

    ### Hypo ###
    hypo_meth = np.logical_and(pt_values == 0, nc_values == 1)
    hypo_meth = np.logical_and(hypo_meth, solo_col)
    print("hypo meth: %s" % (np.sum(hypo_meth) / np.sum(solo_col)))
    f_context_3 = context_3[np.where(hypo_meth)]
    f_context_2 = context_2[np.where(hypo_meth)]
    f_context_1 = context_1[np.where(hypo_meth)]
    f_orph_info = orph_info[np.where(hypo_meth)]
    f_weak = weak[np.where(hypo_meth)]
    f_strong = strong[np.where(hypo_meth)]
    f_weak_or_strong = weak_or_strong[np.where(hypo_meth)]

    final_table = np.hstack((f_orph_info,
                             f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
                             f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))

    end_df = pd.DataFrame(final_table,
                          columns=ORPH_COLS +
                                  ["weak", "strong", "weak_or_strong", "context_1", "context_2", "context_3"])
    profile = ProfileReport(end_df, title='hypo methylation change report - not solo',
                            html={'style': {'full_width': True}})
    profile.to_file(output_file=f"profile_hypo_{patient}_{chromosome}_not_solo.html")


def collect_data_solo(df, chr_info, patient, chromosome):
    context_3 = np.array(format_cpg_context_map.get_context_as_str_for_chr(chr_info))
    context_2 = np.array(format_cpg_context_map.get_context_as_str_for_chr_2(chr_info))
    context_1 = np.array(format_cpg_context_map.get_context_as_str_for_chr_1(chr_info))
    orph_info = chr_info[:, 1:14]
    weak = format_cpg_context_map.get_weak_column(chr_info)
    strong = format_cpg_context_map.get_strong_column(chr_info)
    weak_or_strong = np.logical_xor(weak, strong)
    solo_col = format_cpg_context_map.get_orph_35_column(chr_info)
    solo_col[solo_col > 0] = 1
    solo_col = np.logical_not(solo_col)

    nc_index = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    pt_index = [cell_id for cell_id in df.index if cell_id.startswith('PT')]

    nc_values = df.loc[nc_index, :].mean(axis=0, skipna=True)
    pt_values = df.loc[pt_index, :].mean(axis=0, skipna=True)

    pt_values[pt_values >= HIGH_T] = 1
    pt_values[pt_values <= LOW_T] = 0
    nc_values[nc_values <= LOW_T] = 0
    nc_values[nc_values >= HIGH_T] = 1

    ### No change ###
    kept_meth = np.logical_and(pt_values == 1, nc_values == 1)
    kept_meth = np.logical_and(solo_col, kept_meth)
    print("No changes: %s" % (np.sum(kept_meth) / np.sum(solo_col)))
    f_context_3 = context_3[np.where(kept_meth)]
    f_context_2 = context_2[np.where(kept_meth)]
    f_context_1 = context_1[np.where(kept_meth)]
    f_orph_info = orph_info[np.where(kept_meth)]
    f_weak = weak[np.where(kept_meth)]
    f_strong = strong[np.where(kept_meth)]
    f_weak_or_strong = weak_or_strong[np.where(kept_meth)]

    final_table = np.hstack((f_orph_info,
                             f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
                             f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))

    end_df = pd.DataFrame(final_table, columns=ORPH_COLS +
                                               ["weak", "strong", "weak_or_strong", "context_1", "context_2",
                                                "context_3"])
    profile = ProfileReport(end_df, title='No methylation change report - solo',
                            html={'style': {'full_width': True}})
    profile.to_file(output_file=f"profile_no_met_change_{patient}_{chromosome}_solo.html")

    ### Hypo ###
    hypo_meth = np.logical_and(pt_values == 0, nc_values == 1)
    hypo_meth = np.logical_and(hypo_meth, solo_col)
    print("hypo meth: %s" % (np.sum(hypo_meth) / np.sum(solo_col)))
    f_context_3 = context_3[np.where(hypo_meth)]
    f_context_2 = context_2[np.where(hypo_meth)]
    f_context_1 = context_1[np.where(hypo_meth)]
    f_orph_info = orph_info[np.where(hypo_meth)]
    f_weak = weak[np.where(hypo_meth)]
    f_strong = strong[np.where(hypo_meth)]
    f_weak_or_strong = weak_or_strong[np.where(hypo_meth)]

    final_table = np.hstack((f_orph_info,
                             f_weak[:, None], f_strong[:, None], f_weak_or_strong[:, None],
                             f_context_1[:, None], f_context_2[:, None], f_context_3[:, None]))

    end_df = pd.DataFrame(final_table,
                          columns=ORPH_COLS +
                                  ["weak", "strong", "weak_or_strong", "context_1", "context_2", "context_3"])
    profile = ProfileReport(end_df, title='hypo methylation change report - solo',
                            html={'style': {'full_width': True}})
    profile.to_file(output_file=f"profile_hypo_{patient}_{chromosome}_solo.html")


def collect_data(df, chr_info, patient, chromosome):
    # Seq information
    orph_info = chr_info[:, 1:14]
    context_1 = np.array(format_cpg_context_map.get_context_as_str_for_chr_1(chr_info))
    context_2 = np.array(format_cpg_context_map.get_context_as_str_for_chr_2(chr_info))
    context_3 = np.array(format_cpg_context_map.get_context_as_str_for_chr(chr_info))

    # weak\strong
    weak = format_cpg_context_map.get_weak_column(chr_info)
    strong = format_cpg_context_map.get_strong_column(chr_info)
    weak_or_strong = np.logical_or(weak, strong)
    not_strong_or_weak = np.logical_not(weak_or_strong)

    # solo\ not solo
    solo_col = format_cpg_context_map.get_orph_35_column(chr_info)
    solo_col[solo_col > 0] = 1
    is_solo = solo_col.astype(np.bool)
    not_solo = np.logical_not(is_solo).astype(np.bool)

    # Format methylation level values
    nc_index = [cell_id for cell_id in df.index if cell_id.startswith('NC')]
    pt_index = [cell_id for cell_id in df.index if cell_id.startswith('PT')]

    nc_values = df.loc[nc_index, :].mean(axis=0, skipna=True)
    pt_values = df.loc[pt_index, :].mean(axis=0, skipna=True)

    pt_values[pt_values >= HIGH_T] = 1
    pt_values[pt_values <= LOW_T] = 0
    nc_values[nc_values <= LOW_T] = 0
    nc_values[nc_values >= HIGH_T] = 1

    info_file = open("information_data.txt", "w")
    valid_indexes = np.logical_and(np.logical_or(pt_values == 1, pt_values == 0),
                                   np.logical_or(nc_values == 1, nc_values == 0))

    info_file.write("Perc valid locations: %s\n" % (np.sum(valid_indexes) / valid_indexes.size))
    info_file.write("Perc with meth orig: %s\n" % (np.sum(nc_values == 1) / nc_values.size))

    kept_meth = np.logical_and(pt_values == 1, nc_values == 1)
    kept_meth = np.logical_and(kept_meth, valid_indexes)

    info_file.write("Perc kept: %s\n" % (np.sum(kept_meth) / kept_meth.size))

    lost_meth = np.logical_and(pt_values == 0, nc_values == 1)
    lost_meth = np.logical_and(lost_meth, valid_indexes)

    info_file.write("Perc lost: %s\n" % (np.sum(lost_meth) / kept_meth.size))

    valid_solo = np.logical_and(is_solo, kept_meth)
    info_file.write("Perc of solo kept meth: %s\n" % (np.sum(valid_solo) / np.sum(valid_solo)))

    valid_not_solo = np.logical_and(not_solo, valid_indexes)
    info_file.write("Perc of not solo kept meth: %s\n" % (np.sum(valid_not_solo) / np.sum(valid_not_solo)))

    kept_strong = np.logical_and(kept_meth, strong)
    kept_weak = np.logical_and(kept_meth, weak)
    kept_nothing = np.logical_and(kept_meth, not_strong_or_weak)
    info_file.write("Perc of strong kept meth: %s\n" % (np.sum(kept_strong) / np.sum(strong)))
    info_file.write("Perc of weak kept meth: %s\n" % (np.sum(kept_weak) / np.sum(weak)))
    info_file.write("Perc of else kept meth: %s\n" % (np.sum(kept_nothing) / np.sum(not_strong_or_weak)))

    strong_solo = np.logical_and(solo_col, strong)
    strong_solo = np.logical_and(strong_solo, valid_indexes)
    kept_strong_solo = np.logical_and(strong_solo, kept_meth)
    info_file.write(
        "Perc of strong solo kept from all strong solo: %s\n" % (np.sum(kept_strong_solo) / np.sum(
            strong_solo)))

    weak_solo = np.logical_and(solo_col, weak)
    weak_solo = np.logical_and(weak_solo, valid_indexes)
    kept_weak_solo = np.logical_and(weak_solo, kept_meth)
    info_file.write("Perc of weak solo kept from all weak solo: %s\n" % (np.sum(kept_weak_solo) / np.sum(
        weak_solo)))

    strong_not_solo = np.logical_and(not_solo, strong)
    strong_not_solo = np.logical_and(strong_not_solo, valid_indexes)
    kept_strong_not_solo = np.logical_and(strong_not_solo, kept_meth)
    info_file.write(
        "Perc of strong not solo kept from all strong not solo: %s\n" %
        (np.sum(kept_strong_not_solo) / np.sum(strong_not_solo)))

    weak_not_solo = np.logical_and(not_solo, weak)
    weak_not_solo = np.logical_and(weak_not_solo, valid_indexes)
    kept_weak_not_solo = np.logical_and(weak_not_solo, kept_meth)
    info_file.write("Perc of weak not solo kept from all weak not solo: %s\n" %
                    (np.sum(kept_weak_not_solo) / np.sum(weak_not_solo)))

    for context in np.unique(context_1):
        info_file.write("\ncontext:%s\n" % context)
        context_index = np.logical_and(context_1 == context, valid_indexes)
        context_kept = np.logical_and(context_index, kept_meth)
        context_kept_in_solo = np.logical_and(context_kept, is_solo)
        context_in_solo = np.logical_and(context_index, is_solo)
        context_kept_in_not_solo = np.logical_and(context_kept, not_solo)
        context_in_not_solo = np.logical_and(context_index, not_solo)
        context_in_else = np.logical_and(context_index, not_strong_or_weak)
        context_kept_in_else = np.logical_and(context_kept, not_strong_or_weak)

        info_file.write("Context %s kept methylation: %s of the times\n"
                        % (context, np.sum(context_kept) / np.sum(context_index)))
        info_file.write("Context %s kept methylation in solo: %s of the times\n"
                        % (context, np.sum(context_kept_in_solo) / np.sum(context_in_solo)))
        info_file.write("Context %s kept methylation in not solo: %s of the times\n"
                        % (context, np.sum(context_kept_in_not_solo) / np.sum(context_in_not_solo)))
        info_file.write("Context %s kept methylation in else: %s of the times\n"
                        % (context, np.sum(context_kept_in_else) / np.sum(context_in_else)))

    info_file.close()


def main():
    input_files, output_dir = format_args()
    pmd_context_map = handle_pmds.get_pmd_context_map()

    for file_path in tqdm(input_files):
        patient, chromosome = CPG_FORMAT_FILE_RE.findall(file_path)[0]

        df = pd.read_pickle(file_path)
        pmd_df = handle_pmds.get_pmd_df(df, chromosome)
        # collect_data_global(pmd_df, pmd_context_map[chromosome], patient, chromosome)
        # collect_data_solo(pmd_df, pmd_context_map[chromosome], patient, chromosome)
        # collect_data_not_solo(pmd_df, pmd_context_map[chromosome], patient, chromosome)
        # collect_data(pmd_df, pmd_context_map[chromosome], patient, chromosome)


if __name__ == '__main__':
    main()
