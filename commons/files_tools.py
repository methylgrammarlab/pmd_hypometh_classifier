"""
A set of common useful functions to do data manipulations
"""

import glob
import os
import pickle
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
        chr_name = consts.CHR_FULL_NAME_RE.findall(os.path.basename(f))[0]
        data = load_compressed_pickle(f)
        context_map[chr_name] = data[:, 0] if only_locations else data

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


def convert_paths_list_to_chromosome_based_dict(files):
    """
    Convert a lis of paths for files to a dictionary were key is the chromosome and the value is a list of
    files matching this chromosome
    :param files: Files list
    :return: A dictionary
    :rtype: dict
    """
    files_dict = {}
    for file_path in files:
        try:
            patient, chromosome = consts.PATIENT_CHR_NAME_RE.findall(file_path)[0]
            if chromosome not in files_dict:
                files_dict[chromosome] = []

            files_dict[chromosome].append(file_path)

        except:
            print("Error: can't use file %s because it doesn't match the format" % file_path)

    return files_dict


def get_files_to_work(folder_file_path, pattern):
    """
    Get a list of files matching the required pattern
    :param folder_file_path: The path for the file or folder
    :param pattern: A suffix for folders, will be used for the glob
    :return: A list with all the paths matching this format
    :rtype: list{str}
    """
    if os.path.isdir(folder_file_path) or "*" in folder_file_path:
        file_path = os.path.join(folder_file_path, pattern)
        all_file_paths = glob.glob(file_path)

    else:
        all_file_paths = [folder_file_path]

    return all_file_paths
