# Trying to compare servral databases to see if they act the same

import matplotlib.pyplot as plt
import seaborn as sns

from classifier import convert_sc_dataset_to_nn
from format_files.handle_pmds import filtered_out_non_pmd

plt.style.use('seaborn')

import pandas as pd

ben_path = r"H:\Study\university\Computational-Biology\Year " \
           r"3\Projects\proj_scwgbs\resource\other\methbase_Berman_2012_Human_ColonCancer_chr1.bed"
coad_path = r"H:\Study\university\Computational-Biology\Year " \
            r"3\Projects\proj_scwgbs\resource\other\TCGA_COAD_A00R_chr1.bed"
our_path_crc1_chr1 = r"H:\Study\university\Computational-Biology\Year " \
                     r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC01\all_cpg_ratios_CRC01_chr1.dummy.pkl.zip"
our_path_crc11_chr1 = r"H:\Study\university\Computational-Biology\Year " \
                      r"3\Projects\proj_scwgbs\resource\cpg_format\filtered_by_bl_and_cpgi\CRC11" \
                      r"\all_cpg_ratios_CRC11_chr1.dummy.pkl.zip"
our_combine_path = None


def get_seq_info(ind, chromosome):
    seq = []
    for i in ind:
        seq.append(convert_sc_dataset_to_nn.get_seq_for_cpg(chromosome, i, 150))

    return seq


def crc11_vs_coad():
    coad_data = pd.read_csv(coad_path, names=["chr", "start", "end", "meth", "else"], delimiter="\t",
                            usecols=["start", "meth"], index_col="start")
    crc11_data = pd.read_pickle(our_path_crc11_chr1)

    region_cell_ids = [cell_id for cell_id in crc11_data.index if not cell_id.startswith('NC')]
    crc11_data = crc11_data.loc[region_cell_ids, :].mean()
    crc11_data = crc11_data[~pd.isnull(crc11_data)]

    crc11_data.index = crc11_data.index - 1
    crc11_data = crc11_data.to_frame(name="meth")

    # Only PMD
    coad_data = filtered_out_non_pmd(coad_data, "1")
    crc11_data = filtered_out_non_pmd(crc11_data, "1")

    coad_index = set(coad_data.index.values)
    crc11_index = set(crc11_data.index)

    crc11_and_coad_index = list(coad_index & crc11_index)
    crc11_and_coad_index.sort()
    #
    # sns.distplot(coad_data, hist=False, kde=True, kde_kws={'linewidth': 2}, label="coad")
    # sns.distplot(crc11_data, hist=False, kde=True, kde_kws={'linewidth': 2}, label="crc11")
    # plt.title("Methylation dist across chr1 for crc11 and coad")
    # plt.savefig("meth_dist_chr1_crc11_coad.png")
    #
    # plt.close()

    # for scatter plot
    crc11_coad = pd.DataFrame()
    crc11_coad["start"] = crc11_and_coad_index
    crc11_coad = crc11_coad.set_index("start")
    crc11_coad["coad"] = coad_data.loc[crc11_coad.index].values
    crc11_coad["crc11_meth"] = crc11_data.loc[crc11_coad.index].values
    crc11_coad["sequence"] = get_seq_info(list(crc11_coad.index), 1)
    crc11_coad["small_seq"] = crc11_coad["sequence"].str[73:77]

    density_rows = crc11_coad[crc11_coad["sequence"].str.count("CG") == 3]

    sns_plot = sns.jointplot(x="coad", y="crc11_meth", data=density_rows, kind="kde")
    plt.title("crc11 vs COAD cpg==3 methylation level")
    sns_plot.savefig("scatter_meth_crc11_coad_dens_cpg==3.png")

    plt.close()


def crc1_vs_berman():
    ben_data = pd.read_csv(ben_path, names=["chr", "start", "end", "meth", "else"], delimiter="\t",
                           usecols=["start", "meth"], index_col="start")
    crc1_data = pd.read_pickle(our_path_crc1_chr1)

    region_cell_ids = [cell_id for cell_id in crc1_data.index if not cell_id.startswith('NC')]
    crc1_data = crc1_data.loc[region_cell_ids, :].mean()
    crc1_data = crc1_data[~pd.isnull(crc1_data)]

    crc1_data.index = crc1_data.index - 1
    crc1_data = crc1_data.to_frame(name="meth")

    # Only PMD
    ben_data = filtered_out_non_pmd(ben_data, "1")
    crc1_data = filtered_out_non_pmd(crc1_data, "1")

    ben_index = set(ben_data.index.values)
    crc1_index = set(crc1_data.index)

    crc1_and_ben_index = list(ben_index & crc1_index)
    crc1_and_ben_index.sort()
    #
    # sns.distplot(ben_data, hist=False, kde=True, kde_kws={'linewidth': 2}, label="berman")
    # sns.distplot(crc1_data, hist=False, kde=True, kde_kws={'linewidth': 2}, label="crc1")
    # plt.title("Methylation dist across crc1 and berman- chr 1")
    # plt.savefig("meth_dist_chr1_crc11_bermna.png")
    #
    # plt.close()

    # for scatter plot
    crc1_berman = pd.DataFrame()
    crc1_berman["start"] = crc1_and_ben_index
    crc1_berman = crc1_berman.set_index("start")
    crc1_berman["berman"] = ben_data.loc[crc1_berman.index].values
    crc1_berman["crc01_meth"] = crc1_data.loc[crc1_berman.index].values
    crc1_berman["sequence"] = get_seq_info(list(crc1_berman.index), 1)
    crc1_berman["small_seq"] = crc1_berman["sequence"].str[73:77]
    density_rows = crc1_berman[crc1_berman["sequence"].str.count("CG") == 1]

    sns_plot = sns.jointplot(x="berman", y="crc01_meth", data=density_rows, kind="kde")
    plt.title("crc01 vs berman cpg==1 methylation level")

    sns_plot.savefig("scatter_meth_crc01_berman_dens_cpg=1.png")
    plt.close()


if __name__ == '__main__':
    crc1_vs_berman()
    crc11_vs_coad()
