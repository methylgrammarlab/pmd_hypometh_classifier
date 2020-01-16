import datetime
import glob
import os
import pickle
import re
import sys
import zlib

import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
import commons.consts as consts

COLUMNS = ["read_number", "counter"]


def save_as_compressed_pickle(output_file, data):
    """
    Save the data as pickle after compressing it with zlib
    :param output_file: The path of the output file
    :type output_file: str
    :param data: The data to save
    """
    with open(output_file, "wb") as patient_file:
        patient_file.write(zlib.compress(pickle.dumps(data, pickle.HIGHEST_PROTOCOL), 9))


def load_compressed_pickle(file_path):
    """
    Load a file, decompress it and load the data as if it was pickled
    The way to read data which was saved with the `save_as_compressed_pickle` function
    :param file_path: The path of the file to upload
    :type file_path: str
    :return: The data
    """
    with open(file_path, "rb") as data_file:
        data = data_file.read()
        decompressed_data = zlib.decompress(data)
        formatted_data = pickle.loads(decompressed_data)

        return formatted_data


def load_compressed_pickle_not_zlib(file_path):
    """
    Load a file, decompress it and load the data as if it was pickled
    The way to read data which was saved with the `save_as_compressed_pickle` function
    :param file_path: The path of the file to upload
    :type file_path: str
    :return: The data
    """
    with open(file_path, "rb") as data_file:
        data = data_file.read()
        formatted_data = pickle.loads(data)

        return formatted_data


def init_slurm(func):
    """
    Some simple prints when running code in slurm
    :param func:
    :return:
    """
    print("Run the script: %s" % " ".join(sys.argv))
    start_time = datetime.datetime.now()
    print("Start time: %s" % start_time)
    try:
        func()
    except Exception as ex:
        print("Run ended due to exception: %s" % ex)
    else:
        print("Run ended with no exceptions")
    end_time = datetime.datetime.now()
    delta = end_time - start_time
    print("End time:%s\nTotal time is:%s seconds" % (end_time, delta.total_seconds()))


def counter_to_csv(counter, output_path):
    """
    Save a counter obj to csv
    :param counter:
    :param output_path:
    :return:
    """
    counter_df = pd.DataFrame.from_dict(counter, orient='index').reset_index()
    counter_df.columns = COLUMNS
    counter_df.to_csv(output_path)


def get_all_cpg_locations_across_chr(full_name=False, full_data=False):
    """
    Get a dictionary of chr name (number) and a list of all the locations of CpG
    """
    all_files = glob.glob(os.path.join(consts.ALL_SEQ_PATH, "*.pickle.zlib"))
    chr_dict = {}
    for f in all_files:
        chr_name = re.findall("\d+", os.path.basename(f))[0]
        data = load_compressed_pickle(f)
        if full_name and not full_data:
            chr_dict["chr%s" % chr_name] = data[:, 0]
        elif not full_data and not full_data:
            chr_dict[chr_name] = data[:, 0]
        elif full_data and full_data:
            chr_dict["chr%s" % chr_name] = data
        else:
            chr_dict[chr_name] = data

    return chr_dict
