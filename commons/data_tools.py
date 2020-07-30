"""
A set of common useful functions to do data manipulations
"""
import statistics

import matplotlib.pyplot as plt
import pandas as pd


def counter_to_csv(counter, output_path, columns=["read_number", "counter"]):
    """
    Save a counter obj to csv
    :param columns: The columns of the csv
    :param counter: The counter obj
    :param output_path: The output path for the csv
    """
    counter_df = pd.DataFrame.from_dict(counter, orient='index').reset_index()
    counter_df.columns = columns
    counter_df.to_csv(output_path, index=False)


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
    items = []
    for item in counter.items():
        items += [item[0]] * int(item[1])

    return statistics.median(items)


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
    items = []
    for i in range(len(value_list)):
        items += [value_list[i]] * each_value_size[i]

    return items


def numpy_counter_to_dict(counter):
    """
    Convert a counter fill with numpy value to a dict of str
    :param counter: The counter obj
    :return: A dict of string representing the counter
    """
    return {str(item[0]): item[1] for item in counter.items()}


def dict_to_histogram(counter_or_dict, color='#0E74E3', xlabel=None, ylabel=None, title=None, save_path=None):
    """
    Convert a dictionary or a counter into an histogram and show or save it
    :param counter_or_dict: The counter\dict obj
    :param color: The color of the hist columns
    :param xlabel: The xlabel or None
    :param ylabel: The ylabel or None
    :param title: The title or None
    :param save_path: The path if we want to save the histogram or None if only to show
    """
    labels, values = counter_or_dict.keys(), counter_or_dict.values()
    plt.bar(x=labels, height=values, tick_label=labels, color=color)

    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)

    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

    plt.close()


