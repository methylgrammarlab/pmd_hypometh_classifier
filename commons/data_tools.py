import statistics

import pandas as pd

from commons.files_tools import COLUMNS


def counter_to_csv(counter, output_path):
    """
    Save a counter obj to csv
    :param counter:
    :param output_path:
    """
    counter_df = pd.DataFrame.from_dict(counter, orient='index').reset_index()
    counter_df.columns = COLUMNS
    counter_df.to_csv(output_path)


def mean_of_counter_obj(counter):
    """
    Calculate the weighted mean of a counter obj
    :param counter: A counter obj
    """
    sum_of_numbers = 0
    for i in counter.items():
        sum_of_numbers += i[0] * i[1]

    count = sum(counter.values())
    mean = sum_of_numbers / count
    return mean


def median_of_counter_obj(counter):
    """
    Calculate the weighted median of a counter obj
    :param counter: A counter obj
    """
    l = []
    for item in counter.items():
        l += [item[0]] * int(item[1])

    return statistics.median(l)


def extend_lists(value_list, each_value_size):
    """
    Extend a list of values to full a specific size requirements
    Example:
        value_list = [1,2,3], each_value_size = [1,3,5]
        return:  [1,2,2,2,3,3,3,3,3]

    :param value_list: The value to populate the new list
    :param each_value_size: The number of times each value should appear
    :return: The new list
    """
    l = []
    for i in range(len(value_list)):
        l += [value_list[i]] * each_value_size[i]

    return l


def numpy_counter_to_dict(counter):
    """
    Convert a counter fill with numpy value to a dict of str
    :param counter: The counter obj
    :return: A dict of string representing the counter
    """
    return {str(item[0]): item[1] for item in counter.items()}