import pickle
import zlib


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
