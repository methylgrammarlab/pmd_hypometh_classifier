#!/cs/usr/liorf/PycharmProjects/proj_scwgbs/venv/bin python
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))

patient = R"C:\Users\liorf\OneDrive\Documents\University\year 3\Project\proj_scwgbs\stats\CRC11"
format = 'NC_pairwise_coverage_*_*_numNC_*_window_500.dummy.pkl.zip'

def main():
    dict = {}
    all_files = glob.glob(os.path.join(patient, format))
    num_of_cells = 15
    sums = np.zeros((num_of_cells, num_of_cells))
    for path in all_files:
        df = pd.read_pickle(path)
        ones_df = np.where(df > 0, 1, 0)
        ones_df[:, -1] = df.iloc[:, -1]
        for i in range(num_of_cells):
            sums[i, :] = np.sum(ones_df[np.where(ones_df[:, -1] == i)], axis=0)[:-1]

        pass
    index = ['%d_nc' % i for i in range(15)]
    ax = pd.DataFrame(sums, index=index).plot.bar()
    plt.show()





if __name__ == '__main__':
    main()