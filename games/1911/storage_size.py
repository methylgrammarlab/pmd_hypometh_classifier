import csv
import glob
import pickle
import re
import sys
import h5py
import os
import numpy as np
import zlib

REGEX = re.compile(".+(CRC\d+)_(\w+)_(\d+).singleC.cpg")
SUFFIX = "*.singleC.cpg.txt"

chr16 = []
chr_dict = {}
chr_np = {}


def read_file(fin):
    with open(fin) as f:
        csv_file = csv.DictReader(f, delimiter="\t")
        for _ in csv_file.reader:
            break

        for line in csv_file.reader:
            chr = line[0]
            if chr not in chr_dict:
                chr_dict[line[0]] = {}

            pos = int(line[1])
            others = [int(i) for i in line[4:-3]]
            if line[3] == "-":
                pos -= 1

            if pos not in chr_dict[chr]:
                chr_dict[chr][pos] = np.array(others)

            else:
                chr_dict[chr][pos] += np.array(others)

    for chr in chr_dict:
        chr_array = np.array([np.insert(chr_dict[chr][pos], 0, pos) for pos in chr_dict[chr]],dtype=np.int)
        chr_np[chr] = chr_array

    return chr_np


if __name__ == '__main__':

    all_file = glob.glob(os.path.join(sys.argv[1],SUFFIX))
    print(len(all_file))
    patients = {}
    i = 0
    for f in all_file:
        i += 1
        name, cell, num = REGEX.findall(f)[0]
        print(i)
        if name not in patients:
            patients[name] = []

        chrs_dict = read_file(f)
        patients[name].append(chrs_dict)

    with open("one.pickle", "wb") as one:
        pickle.dump(chr_dict, one)

    with open("one.pickle.zlib","wb") as z:
        z.write(zlib.compress(pickle.dumps(chr_dict, pickle.HIGHEST_PROTOCOL),9))
    with h5py.File('one.hdf5', 'w') as f:
        f.create_dataset("one", data=one)




    for p in patients:
        with open("%s.pickle" % p, "w") as pa_f:
            pickle.dumps(patients[p], pa_f)

        with h5py.File('%s.hdf5' % p, 'w') as f:
            f.create_dataset("%s" % p, data=patients[p])

    with open("all.pickle", "wb") as all_f:
        pickle.dump(patients, all_f)

    with h5py.File('all.hdf5', 'w') as f:
        f.create_dataset("all", data=patients)
