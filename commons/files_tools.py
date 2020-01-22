import glob
import os
import pickle
import re
import sys
import zlib

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


def get_cpg_context_map(drop_chr_prefix=False, get_full_mapping=False):
    """
    Get a dictionary of chr name (number) and a list of all the locations of CpG
    @:param get_full_mapping: If true only give the locations
    """
    all_files = glob.glob(os.path.join(consts.ALL_SEQ_PATH, "*.pickle.zlib"))
    chr_dict = {}
    for f in all_files:
        chr_name = re.findall("\d+", os.path.basename(f))[0]
        data = load_compressed_pickle(f)
        if drop_chr_prefix and not get_full_mapping:
            chr_dict["chr%s" % chr_name] = data[:, 0]
        elif not get_full_mapping and not get_full_mapping:
            chr_dict[chr_name] = data[:, 0]
        elif get_full_mapping and get_full_mapping:
            chr_dict["chr%s" % chr_name] = data
        else:
            chr_dict[chr_name] = data

    return chr_dict
