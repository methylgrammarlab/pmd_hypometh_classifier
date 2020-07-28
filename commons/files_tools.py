"""
A set of common useful functions to do data manipulations
"""

import glob
import os
import pickle
import re
import sys
import zlib

import pandas as pd

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())
from commons import consts

context_map = {}


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


def load_pickle(file_path):
    """
    Load a pickle file based on a path
    :param file_path: The path of the file to upload
    :type file_path: str
    :return: The data
    """
    with open(file_path, "rb") as data_file:
        return pickle.load(data_file)


def save_pickle(file_path, data):
    """
    Save pickle data
    :param file_path: The path to save
    :param data: The data to save
    """
    with open(file_path, "wb") as f:
        pickle.dump(data, f)

def get_cpg_context_map(only_locations=False, load_with_path=consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI):
    """
    Get a dictionary of chr name (chr16 style) and the indexs of all cpgs based on the context map file
    If only_location=False we only get a np.array of the locations, if it's true we get the entire
    informaiton we have about each location: CpGX density for several lengths, weak or strong,
    context flanking
    :param only_locations: If true only give the locations of cpg
    :param load_with_path: The path to the context map file
    :rtype: dict
    """
    if context_map != {} and consts.CONTEXT_MAP_FILTERED_NO_BL_CPGI == load_with_path:
        return context_map

    all_files = glob.glob(os.path.join(load_with_path, "*.pickle.zlib"))
    for f in all_files:
        chr_name = re.findall("chr\d+", os.path.basename(f))[0]
        data = load_compressed_pickle(f)
        context_map["%s" % chr_name] = data[:, 0] if only_locations else data

    return context_map


def load_bedgraph(bedgraph_path):
    """
    Load a bedgraph file to a data frame
    :param bedgraph_path: The path for the bedgraph file
    :return: A data frame of representing the bedgraph
    :rtype: pd.DataFrame
    """
    return pd.read_csv(bedgraph_path, sep="\t", names=["chr", "start", "end", "coverage"],
                       usecols=["start", "coverage"], index_col="start")
