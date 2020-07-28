import csv

import numpy as np
from sklearn.metrics import r2_score

FILE_NAME = r"external/GSM3153400_scTrioSeq2Met_CRC11_NC_418.singleC.cpg.txt"

START_INDEX = 6500000
END_INDEX = 6700000
POSITION_INDEX = 0
STRAND_INDEX = "2"

chr16 = []
pos_dict = {}


def main():
    with open(FILE_NAME) as f:
        csv_file = csv.DictReader(f, delimiter="\t")

        for line in csv_file.reader:
            if line[0] == "chr16":
                pos = int(line[1])
                others = [int(i) for i in line[4:-3]]
                if line[3] == "-":
                    pos -= 1

                if pos not in pos_dict:
                    pos_dict[pos] = np.array(others)

                else:
                    pos_dict[pos] += np.array(others)

    data_array = np.array([np.insert(pos_dict[pos], 0, pos) for pos in pos_dict])
    indexes = np.where(
        np.logical_and(data_array[:, POSITION_INDEX] > START_INDEX,
                       data_array[:, POSITION_INDEX] < END_INDEX))
    pmd = data_array[indexes]
    ratio = pmd[:, 2] / pmd[:, 1]
    # pmd_with_ratio = np.column_stack((pmd, ratio))
    # pmd_with_ratio[:,0] = pmd_with_ratio[:,0].astype(np.int)

    r2 = []
    for index in range(ratio.size):
        print(index)
        ratio_copy = np.copy(ratio)
        value = ratio_copy[index]
        ratio_copy_deleted = np.delete(ratio_copy, index)
        value_stretched = np.tile(value, ratio_copy_deleted.size)
        score = r2_score(ratio_copy_deleted, value_stretched)
        r2.append(score)


if __name__ == '__main__':
    main()
