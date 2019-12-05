import os
import pickle
import zlib

from format_files.format_sc_data import FILE_DETAILS_RE


def save_as_compressed_pickle(file_path, data):
    with open(file_path, "wb") as patient_file:
        patient_file.write(zlib.compress(pickle.dumps(data, pickle.HIGHEST_PROTOCOL), 9))


def load_compressed_pickle(file_path):
    with open(file_path, "rb") as data_file:
        data = data_file.read()
        decompressed_data = zlib.decompress(data)
        formatted_data = pickle.loads(decompressed_data)

        return formatted_data
