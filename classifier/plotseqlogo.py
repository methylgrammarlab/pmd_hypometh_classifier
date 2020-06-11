#######################################################################################
### Code adoped from  https://github.com/c/concise (c) 2016, Žiga Avsec Gagneur lab ### #######################################################################################

import math
import re
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from concise.preprocessing.sequence import DNA, RNA, AMINO_ACIDS
from concise.utils.letters import all_letters
from descartes.patch import PolygonPatch
from shapely import affinity
from shapely.wkt import loads as load_wkt


def standardize_polygons_str(data_str):
    """Given a POLYGON string, standardize the coordinates to a 1x1 grid.
    Input : data_str (taken from above)
    Output: tuple of polygon objects
    """
    # find all of the polygons in the letter (for instance an A
    # needs to be constructed from 2 polygons)
    path_strs = re.findall("\(\(([^\)]+?)\)\)", data_str.strip())

    # convert the data into a numpy array
    polygons_data = []
    for path_str in path_strs:
        data = np.array([
            tuple(map(float, x.split())) for x in path_str.strip().split(",")])
        polygons_data.append(data)

    # standardize the coordinates
    min_coords = np.vstack(data.min(0) for data in polygons_data).min(0)
    max_coords = np.vstack(data.max(0) for data in polygons_data).max(0)
    for data in polygons_data:
        data[:, ] -= min_coords
        data[:, ] /= (max_coords - min_coords)

    polygons = []
    for data in polygons_data:
        polygons.append(load_wkt(
            "POLYGON((%s))" % ",".join(" ".join(map(str, x)) for x in data)))

    return tuple(polygons)


def update_all_letters():
    weak_strong = {
        'AW': """POLYGON ((0 72.9062, 9.625 72.9062, 16.6094 13.7188, 24.9062 52.875, 35.2031 52.875, 43.6094 13.625, 50.5938 72.9062, 60.2031 72.9062, 49.3125 0, 39.9844 0, 30.0781 43.3125, 20.2188 0, 10.8906 0, 0 72.9062))""",
        'CS': """POLYGON ((49.4219 70.4062, 49.4219 60.4062, 44.918 62.9297, 40.4062 64.75, 35.8711 65.8516, 31.2969 66.2188, 27.9883 66.0156, 25.0469 65.4062, 22.4727 64.3906, 20.2656 62.9688, 18.4951 61.2021, 17.2305 59.1523, 16.4717 56.8193, 16.2188 54.2031, 16.3848 51.9253, 16.8828 49.9355, 17.7129 48.2339, 18.875 46.8203, 20.4922 45.603, 22.6875 44.4902, 25.4609 43.4819, 28.8125 42.5781, 33.9844 41.4062, 38.7617 40.04, 42.8438 38.332, 46.2305 36.2822, 48.9219 33.8906, 50.9727 31.0957, 52.4375 27.8359, 53.3164 24.1113, 53.6094 19.9219, 53.2065 15.0273, 51.998 10.75, 51.0917 8.8428, 49.9839 7.0898, 48.6747 5.4912, 47.1641 4.0469, 43.5757 1.6543, 39.2559 -0.0547, 34.2046 -1.0801, 28.4219 -1.4219, 23.2832 -1.1465, 18.1172 -0.3203, 12.9277 1.0566, 7.7188 2.9844, 7.7188 13.4844, 13.1777 10.3867, 18.3359 8.25, 23.3613 7.0078, 28.4219 6.5938, 31.9487 6.8027, 35.0605 7.4297, 37.7573 8.4746, 40.0391 9.9375, 41.854 11.7754, 43.1504 13.9453, 43.9282 16.4473, 44.1875 19.2812, 44.0093 21.854, 43.4746 24.1035, 42.5835 26.0298, 41.3359 27.6328, 39.645 28.9917, 37.4238 30.1855, 34.6724 31.2144, 31.3906 32.0781, 26.125 33.2969, 21.3945 34.5918, 17.3594 36.1797, 14.0195 38.0605, 11.375 40.2344, 9.3652 42.7529, 7.9297 45.668, 7.0684 48.9795, 6.7812 52.6875, 7.1919 57.3276, 8.4238 61.4824, 10.4771 65.1519, 13.3516 68.3359, 16.9067 70.9097, 21.002 72.748, 25.6372 73.8511, 30.8125 74.2188, 35.123 73.9805, 39.6484 73.2656, 44.4082 72.0742, 49.4219 70.4062))""",
        'GS': """POLYGON ((49.4219 70.4062, 49.4219 60.4062, 44.918 62.9297, 40.4062 64.75, 35.8711 65.8516, 31.2969 66.2188, 27.9883 66.0156, 25.0469 65.4062, 22.4727 64.3906, 20.2656 62.9688, 18.4951 61.2021, 17.2305 59.1523, 16.4717 56.8193, 16.2188 54.2031, 16.3848 51.9253, 16.8828 49.9355, 17.7129 48.2339, 18.875 46.8203, 20.4922 45.603, 22.6875 44.4902, 25.4609 43.4819, 28.8125 42.5781, 33.9844 41.4062, 38.7617 40.04, 42.8438 38.332, 46.2305 36.2822, 48.9219 33.8906, 50.9727 31.0957, 52.4375 27.8359, 53.3164 24.1113, 53.6094 19.9219, 53.2065 15.0273, 51.998 10.75, 51.0917 8.8428, 49.9839 7.0898, 48.6747 5.4912, 47.1641 4.0469, 43.5757 1.6543, 39.2559 -0.0547, 34.2046 -1.0801, 28.4219 -1.4219, 23.2832 -1.1465, 18.1172 -0.3203, 12.9277 1.0566, 7.7188 2.9844, 7.7188 13.4844, 13.1777 10.3867, 18.3359 8.25, 23.3613 7.0078, 28.4219 6.5938, 31.9487 6.8027, 35.0605 7.4297, 37.7573 8.4746, 40.0391 9.9375, 41.854 11.7754, 43.1504 13.9453, 43.9282 16.4473, 44.1875 19.2812, 44.0093 21.854, 43.4746 24.1035, 42.5835 26.0298, 41.3359 27.6328, 39.645 28.9917, 37.4238 30.1855, 34.6724 31.2144, 31.3906 32.0781, 26.125 33.2969, 21.3945 34.5918, 17.3594 36.1797, 14.0195 38.0605, 11.375 40.2344, 9.3652 42.7529, 7.9297 45.668, 7.0684 48.9795, 6.7812 52.6875, 7.1919 57.3276, 8.4238 61.4824, 10.4771 65.1519, 13.3516 68.3359, 16.9067 70.9097, 21.002 72.748, 25.6372 73.8511, 30.8125 74.2188, 35.123 73.9805, 39.6484 73.2656, 44.4082 72.0742, 49.4219 70.4062))""",
        'TW': """POLYGON ((0 72.9062, 9.625 72.9062, 16.6094 13.7188, 24.9062 52.875, 35.2031 52.875, 43.6094 13.625, 50.5938 72.9062, 60.2031 72.9062, 49.3125 0, 39.9844 0, 30.0781 43.3125, 20.2188 0, 10.8906 0, 0 72.9062))""",
        }
    all_letters.update(weak_strong)


# ----------------------
update_all_letters()
letter_polygons = {k: standardize_polygons_str(v) for k, v in all_letters.items()}

VOCABS = {"DNA": OrderedDict([("A", "green"),
                              ("C", "blue"),
                              ("G", "orange"),
                              ("T", "red")]),
          "DNAWS": OrderedDict([("AW", "red"),
                                ("CS", "green"),
                                ("GS", "green"),
                                ("TW", "red")]),
          "RNA": OrderedDict([("A", "green"),
                              ("C", "blue"),
                              ("G", "orange"),
                              ("U", "red")]),
          "AA": OrderedDict([('A', '#CCFF00'),
                             ('B', "orange"),
                             ('C', '#FFFF00'),
                             ('D', '#FF0000'),
                             ('E', '#FF0066'),
                             ('F', '#00FF66'),
                             ('G', '#FF9900'),
                             ('H', '#0066FF'),
                             ('I', '#66FF00'),
                             ('K', '#6600FF'),
                             ('L', '#33FF00'),
                             ('M', '#00FF00'),
                             ('N', '#CC00FF'),
                             ('P', '#FFCC00'),
                             ('Q', '#FF00CC'),
                             ('R', '#0000FF'),
                             ('S', '#FF3300'),
                             ('T', '#FF6600'),
                             ('V', '#99FF00'),
                             ('W', '#00CCFF'),
                             ('Y', '#00FFCC'),
                             ('Z', 'blue')]),
          "RNAStruct": OrderedDict([("P", "red"),
                                    ("H", "green"),
                                    ("I", "blue"),
                                    ("M", "orange"),
                                    ("E", "violet")]),
          }
# make sure things are in order
VOCABS["AA"] = OrderedDict((k, VOCABS["AA"][k]) for k in AMINO_ACIDS)
VOCABS["DNA"] = OrderedDict((k, VOCABS["DNA"][k]) for k in DNA)
VOCABS["RNA"] = OrderedDict((k, VOCABS["RNA"][k]) for k in RNA)


def add_letter_to_axis(ax, let, col, x, y, height):
    """Add 'let' with position x,y and height height to matplotlib axis 'ax'.
    """
    if len(let) == 2:
        colors = [col, "white"]
    elif len(let) == 1:
        colors = [col]
    else:
        raise ValueError("3 or more Polygons are not supported")

    for polygon, color in zip(let, colors):
        new_polygon = affinity.scale(
            polygon, yfact=height, origin=(0, 0, 0))
        new_polygon = affinity.translate(
            new_polygon, xoff=x, yoff=y)
        patch = PolygonPatch(new_polygon, edgecolor=color, facecolor=color)
        # patch = PolygonPatch(new_polygon, edgecolor=None, facecolor=color)
        ax.add_patch(patch)
    return


def seqlogo(letter_heights, vocab="DNA", ax=None, yl=None):
    """Make a logo plot
    # Arguments
        letter_heights: "motif length" x "vocabulary size" numpy array
    Can also contain negative values.
        vocab: str, Vocabulary name. Can be: DNA, RNA, AA, RNAStruct.
        ax: matplotlib axis
    """
    ax = ax or plt.gca()

    assert letter_heights.shape[1] == len(VOCABS[vocab])
    x_range = [1, letter_heights.shape[0]]
    pos_heights = np.copy(letter_heights)
    pos_heights[letter_heights < 0] = 0
    neg_heights = np.copy(letter_heights)
    neg_heights[letter_heights > 0] = 0

    for x_pos, heights in enumerate(letter_heights):
        letters_and_heights = sorted(zip(heights, list(VOCABS[vocab].keys())))
        y_pos_pos = 0.0
        y_neg_pos = 0.0
        for height, letter in letters_and_heights:
            color = VOCABS[vocab][letter]
            polygons = letter_polygons[letter]
            if height > 0:
                add_letter_to_axis(ax, polygons, color, 0.5 + x_pos, y_pos_pos, height)
                y_pos_pos += height
            elif height == 0 and ((x_pos == len(letter_heights) / 2 and letter == 'G') or (
                    x_pos == len(letter_heights) / 2 - 1 and letter == 'C')):
                height = 0.01
                add_letter_to_axis(ax, polygons, "gray", 0.5 + x_pos, y_pos_pos, height)
                y_pos_pos += height
            elif height == 0:
                pass
            else:
                add_letter_to_axis(ax, polygons, color, 0.5 + x_pos, y_neg_pos, height)
                y_neg_pos += height

    # if add_hline:
    #     ax.axhline(color="black", linewidth=1)
    # ax.set_xlim(x_range[0] - 1, x_range[1] + 1)

    if yl:
        ax.set_ylim(-yl, yl)


    ax.grid(False)
    # ax.set_xticks(list(range(*x_range)) + [x_range[-1]])
    ax.set_xticks([], [])
    ax.yaxis.set_tick_params(labelsize=8)
    ax.set_aspect(aspect='auto', adjustable='box')
    ax.autoscale_view()


def seqlogo_fig(letter_heights, vocab="DNA", figsize=(10, 2), ncol=1, plot_name=None, yl=None):
    """
    # Arguments
        plot_name: Title of the plot. Can be a list of names
    """
    fig = plt.figure(figsize=figsize)

    if len(letter_heights.shape) == 3:
        #
        n_plots = letter_heights.shape[2]
        nrow = math.ceil(n_plots / ncol)
        if isinstance(plot_name, list):
            assert len(plot_name) == n_plots
    else:
        n_plots = 1
        nrow = 1
        ncol = 1

    for i in range(n_plots):
        if len(letter_heights.shape) == 3:
            w_cur = letter_heights[:, :, i]
        else:
            w_cur = letter_heights
        ax = plt.subplot(nrow, ncol, i + 1)
        plt.tight_layout(h_pad=0.01)

        # plot the motif
        seqlogo(w_cur, vocab, ax, yl)

        # add the title
        if plot_name is not None:
            if n_plots > 0:
                if isinstance(plot_name, list):
                    pln = plot_name[i]
                else:
                    pln = plot_name + " {0}".format(i)
            else:
                pln = plot_name
            ax.set_title(pln)
    return fig
