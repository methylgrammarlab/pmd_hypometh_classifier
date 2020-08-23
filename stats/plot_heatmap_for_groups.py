import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use('seaborn')
sns.set_style('whitegrid')

PATIENTS = ["CRC01", "CRC11", "CRC13", "CRC10"]


# PATIENTS = ["CRC13"]


def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', help='Path of the output folder', required=False,
                        default=os.path.dirname(sys.argv[0]))
    parser.add_argument('--input_file', help='path for the input file', required=True)

    args = parser.parse_args()
    return args


def plot_heatmap(df, title, labelx, labely, output_path):
    """
    Plot a heatmap figure
    :param df: The dataframe
    :param title: The title of the figure
    :param labelx: The x label
    :param labely: The y label
    :param output_path: The output path to save the image
    """
    plt.subplots_adjust(top=0.9)
    min, max = np.nanmin(df.replace(0, np.nan).values), np.nanmax(df.replace(0, np.nan).values)
    midpoint = (max - min) / 2
    sns.heatmap(df, cmap='coolwarm', center=midpoint, vmin=min, vmax=max)
    # sns.heatmap(df, cbar=1, cmap="YlGnBu")
    plt.title(title % (labely, labelx), fontsize=20)
    plt.xlabel("%s" % labelx, fontsize=16)
    plt.ylabel("%s" % labely, fontsize=16)

    plt.savefig(output_path)
    plt.close()


def create_heatmap_matrix(df, x_label, y_label, x_step, y_step):
    """
    Create the heatmap matrixs
    :param df: The dataframe with all the information about the x label and y label
    :param x_label:The label to be in the x axis
    :param y_label:The label to be in the y axis
    :param x_step: The step of each tick in the x
    :param y_step: The step of each tick in the y
    """
    x = np.arange(0, 1, x_step)
    y = np.arange(df[y_label].min(), df[y_label].max(), y_step)
    print("Y length is %s" % len(y))
    data_matrix = np.zeros(shape=(len(y), len(x)))
    size_matrix = np.zeros(shape=(len(y), len(x)))

    row = len(y) - 1
    for i in y:
        column = 0
        for j in x:
            indexes = np.logical_and(np.logical_and(df[x_label] >= j, df[x_label] < j + x_step),
                                     np.logical_and(df[y_label] >= i, df[y_label] < i + y_step))
            temp_df = df[indexes]
            num_of_cpg = temp_df.shape[0]
            perc_of_solo_wcgw = np.sum(temp_df["solo-WCGW"]) / num_of_cpg * 100 if num_of_cpg > 500 \
                else np.nan

            # df.loc[indexes, "soloWCGW_perc"] = perc_of_solo_wcgw
            data_matrix[row, column] = perc_of_solo_wcgw
            size_matrix[row, column] = num_of_cpg
            column += 1

        row -= 1
    x = [float("%.2f" % i) for i in x]
    yr = [float("%.3f" % i) for i in y]

    yr.reverse()
    df_data_matrix = pd.DataFrame(data_matrix, index=yr, columns=x)
    df_size_matrix = pd.DataFrame(size_matrix, index=yr, columns=x)
    return df_data_matrix, df_size_matrix, df

def plot_heatmap_to_scwgbs():
    """
    Plot heatmap for scWGBS (Bian) - we decided not to usedidn't used
    """
    args = parse_input()
    input_path = args.input_file

    full_df = pd.read_csv(input_path)
    full_df["ccpg"] = full_df["sequence"].str.count("CG")
    df = full_df[full_df["ccpg"] < 4]
    del full_df

    df = df[df["orig_meth_avg"] >= 0.7]
    df["small_seq"] = df["sequence"].str[73:77]
    solo_wcgw = np.logical_and(df["ccpg"] == 1, df["small_seq"].str.contains("[AT]CG[AT]", regex=True))
    df["solo-WCGW"] = solo_wcgw
    df["meth_diff"] = df["orig_meth"] - df["meth"]

    patients_df = []
    for patient in PATIENTS:
        patient_df = df[df["patient"] == patient]

        x_label = "meth"
        x_step = 0.01

        y_label = "var"
        y_step = (df[y_label].max() - df[y_label].min()) / 100

        df_data_matrix, df_size_matrix, patient_df = create_heatmap_matrix(patient_df, x_label, y_label,
                                                                           x_step, y_step)

        plot_heatmap(df_data_matrix, title="%s vs %s for %%CpG1-WCGW" + " %s" % patient, labelx=x_label,
                     labely=y_label, output_path=os.path.join(args.output_folder, "%s_vs_%s_perc_cpg_%s.png"
                                                              % (x_label, y_label, patient)))

        plot_heatmap(df_size_matrix, title="%s vs %s total CpGs" + " %s" % patient, labelx=x_label,
                     labely=y_label, output_path=os.path.join(args.output_folder,
                                                              "%s_vs_%s_total_cpg_%s.png" %
                                                              (x_label, y_label, patient)))

        # Save data
        df_data_matrix.to_pickle(os.path.join(args.output_folder, "%s_%s_%s_matrix.pkl" % (x_label, y_label,
                                                                                           patient)))
        df_size_matrix.to_pickle(os.path.join(args.output_folder, "%s_%s_%s_size.pkl" % (x_label, y_label,
                                                                                         patient)))

        patients_df.append(patient_df)


def plot_heatmap_to_bulk():
    """
    Plot the heatmap for the bulk(Zhou_ dataset
    """
    args = parse_input()
    input_path = args.input_file
    full_df = pd.read_pickle(input_path)

    x_label = "meth"
    x_step = 0.01

    y_label = "coveriance"
    y_step = 0.001

    full_df["ccpg"] = full_df["sequence"].str.count("CG")
    df = full_df[full_df["ccpg"] < 4]
    df = df[df["orig_meth"] >= 0.7]
    df["small_seq"] = df["sequence"].str[73:77]
    solo_wcgw = np.logical_and(df["ccpg"] == 1, df["small_seq"].str.contains("[AT]CG[AT]", regex=True))
    df["solo-WCGW"] = solo_wcgw
    df["meth_diff"] = df["orig_meth"] - df["meth"]

    df_data_matrix, df_size_matrix, df = create_heatmap_matrix(df, x_label, y_label, x_step, y_step)

    plot_heatmap(df_data_matrix, title="%s vs %s for %%CpG1-WCGW", labelx=x_label, labely=y_label,
                 output_path=os.path.join(args.output_folder, "%s_vs_%s_perc_cpg.png" % (x_label, y_label)))

    plot_heatmap(df_data_matrix, title="%s vs %s total CpGs", labelx=x_label, labely=y_label,
                 output_path=os.path.join(args.output_folder, "%s_vs_%s_total_cpg.png" % (x_label, y_label)))

    # Save data
    df_data_matrix.to_pickle(os.path.join(args.output_folder, "%s_%s_matrix.pkl" % (x_label, y_label)))
    df_size_matrix.to_pickle(os.path.join(args.output_folder, "%s_%s_size.pkl" % (x_label, y_label)))

if __name__ == '__main__':
    plot_heatmap_to_bulk()
    # plot_heatmap_to_scwgbs()
