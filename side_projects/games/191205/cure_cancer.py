import csv
import glob
import os
import pickle
import sys
import zlib

import numpy as np
import pandas as pd

# ALL_CPG_PATH = r":/vol/sci/bio/data/benjamin.berman/bermanb/projects/genomic-data-misc/cgcontext/cgcontext" \
#               r".hg19/chr16.seq_context.13flanks.txt"
FILE_SUFFIX = "*CRC11*.singleC.cpg.txt"
ALL_CPG_PATH = r"external/chr16.seq_context.13flanks.txt"
PMD_DICT = {"chr16": [(6500000, 6700000), (19000000, 21000000), (8000000, 8200000)]}
pos_dict = {}
POSITION_INDEX = 0


def read_file(file_path):
    with open(file_path) as f:
        csv_file = csv.DictReader(f, delimiter="\t")
        for _ in csv_file.reader:
            break

        for line in csv_file.reader:
            chr = line[0]
            if chr not in PMD_DICT.keys():
                continue

            pos = int(line[1])
            others = [int(i) for i in line[4:6]]
            if line[3] == "-":
                pos -= 1

            if pos not in pos_dict:
                pos_dict[pos] = np.array(others)

            else:
                pos_dict[pos] += np.array(others)

        data_array = np.array([np.insert(pos_dict[pos], 0, pos) for pos in pos_dict]).astype(np.float)

    return data_array


def get_pmd_table(cells, start, end, chr):
    chr_index = np.where(np.logical_and(chr > start, chr < end))
    pmd_full_cpg = chr[chr_index]

    all_pmds = np.full((len(cells) + 1, pmd_full_cpg.size), None, dtype=np.float)
    all_pmds[0, :] = pmd_full_cpg

    i = 1
    for cell in cells:
        cell_indexes = np.where(
            np.logical_and(cell[:, POSITION_INDEX] > start, cell[:, POSITION_INDEX] < end))
        pmd = cell[cell_indexes]
        match = np.isin(pmd_full_cpg, pmd)
        ratio = pmd[:, 2] / pmd[:, 1]
        all_pmds[i, match] = ratio
        i += 1

    return all_pmds


def get_chr_cpg():
    with open(ALL_CPG_PATH, "r") as cpg_file:
        csv_file = csv.DictReader(cpg_file, delimiter="\t")
        for _ in csv_file.reader:
            break

        cpg_positions = [line[2] for line in csv_file.reader]

    return np.array(cpg_positions, dtype=np.int)


def save_data(formatted_data):
    with open("3_data.pickle.zlib", "wb") as z:
        z.write(zlib.compress(pickle.dumps(formatted_data, pickle.HIGHEST_PROTOCOL), 9))


def main():
    no_data = False
    if no_data:
        all_files = glob.glob(os.path.join(sys.argv[1], FILE_SUFFIX))
        formatted_data = [read_file(file_path) for file_path in all_files]
        save_data(formatted_data)
        return
    else:
        with open("3_data.pickle.zlib", "rb") as z:
            data = z.read()
            decompressed_data = zlib.decompress(data)
            formatted_data = pickle.loads(decompressed_data)

            del decompressed_data
            del data

    chr_cpg = get_chr_cpg()
    for pmd_boundaries in PMD_DICT["chr16"]:
        pmd_matrix = get_pmd_table(formatted_data, pmd_boundaries[0], pmd_boundaries[1], chr_cpg)

        # none_cols = np.any(np.isnan(pmd_matrix), axis=0)
        #
        df = pd.DataFrame(data=pmd_matrix[1:, :], columns=pmd_matrix[0, :].astype(np.int))
        cov = df.cov()
        show(cov, pmd_matrix[0, :])

        # Y = pdist(pmd_matrix[1:], pd.DataFrame.corr)
        pass


def show(cor, labels):
    labels = labels.astype(np.str)
    from matplotlib import pyplot as plt
    from matplotlib import cm as cm

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cmap = cm.get_cmap('jet', 30)
    cax = ax1.imshow(cor, interpolation="nearest", cmap=cmap)
    ax1.grid(True)
    plt.title('SOmething')
    ax1.set_xticklabels(labels, fontsize=6)
    ax1.set_yticklabels(labels, fontsize=6)
    fig.colorbar(cax, ticks=[-.5, -.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5])
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    plt.show()


if __name__ == '__main__':
    main()
