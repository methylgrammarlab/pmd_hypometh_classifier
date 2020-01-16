import statistics

import pandas as pd

from commons.files_tools import COLUMNS


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


def mean_of_counter_obj(counter):
    sum_of_numbers = 0
    for i in counter.items():
        sum_of_numbers += i[0] * i[1]

    count = sum(counter.values())
    mean = sum_of_numbers / count
    return mean


def median_of_counter_obj(counter):
    l = []
    for item in counter.items():
        l += [item[1]] * int(item[0])

    return statistics.median(l)


def extend_lists(value_list, each_value_size):
    l = []
    for i in range(len(value_list)):
        l += [value_list[i]] * each_value_size[i]

    return l


def counter_to_dict(counter):
    return {str(item[0]): item[1] for item in counter.items()}