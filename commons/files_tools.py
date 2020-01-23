import glob
import os
import pickle
import re
import sys
import zlib

sys.path.append(os.path.dirname(os.getcwd()))
import commons.consts as consts

COLUMNS = ["read_number", "counter"]

chr_dict = {}


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


def get_cpg_context_map(only_locations=False, load_with_path=consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI):
    """
    Get a dictionary of chr name (number) and a list of all the locations of CpG
    @:param only_locations: If true only give the locations of cpg
    """
    if chr_dict != {}  and consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI == load_with_path:
        return chr_dict

    all_files = glob.glob(os.path.join(load_with_path, "*.pickle.zlib"))
    for f in all_files:
        chr_name = re.findall("chr\d+", os.path.basename(f))[0]
        data = load_compressed_pickle(f)
        chr_dict["%s" % chr_name] = data[:, 0] if only_locations else data

    return chr_dict
