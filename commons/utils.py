"""
Set of functions which help do common actions on data frames of our data
"""
import glob
import os

import numpy as np
import pandas as pd

from commons import consts

ALL = 'ALL'
NC_AND_PT = "PT"
ONLY_PT = "ONLY_PT"
NOT_NC = "ALL_CANCER"
LN_AND_PT = "LN_AND_PT"
ONLY_NC = "NC"


def remove_extreme_cpgs_by_coverage(df, top_low_level_to_remove=5):
    """
    Remove CpG from a data frame which are extreme in coverage
    :param df: The data frame
    :type df: pd.DataFrame
    :param top_low_level_to_remove: The percentage of CpC to remove from top and bottom, keep in mind that
    this is percentage and not value base, meaning if every value is 3 that we will only remove
    top_low_level_to_remove*2 % of the values and not all
    :type top_low_level_to_remove: int
    :return: A new df with the cpg we want to keep
    """
    cpg_coverage = np.sum(~pd.isnull(df), axis=0)
    cpg_coverage = cpg_coverage.sort_values()
    cpg_s = cpg_coverage.shape[0]
    n_to_remove = int(cpg_s * top_low_level_to_remove / 100)

    cpg_to_keep = cpg_coverage.index[n_to_remove:-n_to_remove]
    return df[cpg_to_keep]  # this remove the top_low_level_to_remove lower and top


def get_nc_avg(chromosome, cpg_indexes, nc_meth_files):
    """
    Get the average values of methylation in original cells based on a file already created
    :param chromosome: The chromsome we are looking at
    :param cpg_indexes: The cpg indexes
    :param nc_meth_files: The nc files path
    :return: The average methylation in those cells
    """
    nc_file = glob.glob(os.path.join(nc_meth_files, "*%s.dummy*" % chromosome))[0]
    df = pd.read_pickle(nc_file)
    mdf = df.mean(axis=1)
    return mdf.loc[cpg_indexes]


def cpg_meth_in_cells(patient, chromosome, cpgs_indexes, meth_files_path, sublineage_name):
    """
    Get the methylation ratio of a patient in a chromosome for specific set of cells\sublineage_name
    :param patient: The patient we are looking at
    :param chromosome: The chromosome we are looking at
    :param cpgs_indexes: The of cpgs
    :param meth_files_path: The path for the patient files
    :param sublineage_name: The sublineage name or cell type name
    :return: A df with the indexes matching those cells type
    """
    methylation_file_path = os.path.join(meth_files_path, patient, consts.SCWGBS_FILE_FORMAT %
                                         (patient, chromosome))
    df = pd.read_pickle(methylation_file_path)
    df = filter_df_based_on_region_name(df, region_name=sublineage_name)

    return df.loc[cpgs_indexes]


def get_all_indexes_from_patients_list(patients_information):
    """
    Get a list of all the indexes from the different patient matching from the data frames
    :param patients_information: A list with all the patient, each value here is a tuple with a df in the
    second value
    :return: A set of all the indexes available from all the patients
    :rtype: set
    """
    cpg_locations = None
    for df_tuple in patients_information:
        _, df = df_tuple
        if cpg_locations is not None:
            cpg_locations = set(df.index.values) | cpg_locations
        else:
            cpg_locations = set(df.index.values)

    indexes = list(cpg_locations)
    indexes.sort()
    return indexes


def filter_df_based_on_cells(full_df, region_cells):
    """
    Filter a data frame of cells based on list of cells
    :param full_df: The full original df
    :param region_cells: List of cells to use of the region
    :type region_cells: list[str]
    :return: A df filtered dataframe based based on the region name
    :rtype: pd.DataFrame
    """
    region_cell_ids = [cell_id for cell_id in full_df.index if cell_id.startswith('NC')]
    for sublineage_cell in region_cells:
        region_cell_ids.extend(cell_id for cell_id in full_df.index if cell_id.startswith(sublineage_cell))

    return full_df.loc[region_cell_ids, :]


def filter_df_based_on_region_name(full_df, region_name):
    """
    Filter a data frame of cells based on a region name
    :param full_df: The full original df
    :param region_name: The region name to use, based on the consts
    :type region_name: str
    :return: A df filtered dataframe based based on the region name
    :rtype: pd.DataFrame
    """
    # All cells
    if region_name == ALL:
        region_df = full_df

    elif region_name == ONLY_NC:
        region_cell_ids = [cell_id for cell_id in full_df.index if cell_id.startswith('NC')]
        region_df = full_df.loc[region_cell_ids, :]

    elif region_name == LN_AND_PT:
        region_cell_ids = [cell_id for cell_id in full_df.index if (cell_id.startswith('PT') or
                                                                    cell_id.startswith("LN"))]
        region_df = full_df.loc[region_cell_ids, :]

    elif region_name == NOT_NC:
        region_cell_ids = [cell_id for cell_id in full_df.index if not cell_id.startswith('NC')]
        region_df = full_df.loc[region_cell_ids, :]

    elif region_name == NC_AND_PT:
        region_cell_ids = [cell_id for cell_id in full_df.index if cell_id.startswith('NC')]
        region_cell_ids.extend(cell_id for cell_id in full_df.index if cell_id.startswith("PT"))
        region_df = full_df.loc[region_cell_ids, :]

    elif region_name == ONLY_PT:
        region_cell_ids = [cell_id for cell_id in full_df.index if cell_id.startswith('PT')]
        region_df = full_df.loc[region_cell_ids, :]

    else:
        raise NotImplementedError("No region names: %s" % region_name)

    return region_df


def filter_df_based_on_tuple_list(df, boundaries_list, add_index=False, add_index_name="pmd_index"):
    """
    Filter out a df to only contains indexes which included in the tuple list
    :param df: The df to work with
    :param boundaries_list: A tuple list with boundaries, each tuple is (start,end)
    :param add_index: should we add index of each tuple
    :param add_index_name: The name of the index column we should add
    :return: The df filtered out and maybe with index
    """
    prev_mask = None
    i = 0

    if add_index:
        df[add_index_name] = -1

    # Trying to support two ways of working with df - one is where the columns is the CpG and the other is
    # that the index is the CpG location
    for boundary in boundaries_list:
        i += 1
        start, end = boundary
        try:
            pmd_mask = (df.columns >= start) & (df.columns <= end)
        except Exception:
            pmd_mask = (df.index >= start) & (df.index <= end)

        prev_mask = np.logical_or(pmd_mask, prev_mask) if prev_mask is not None else pmd_mask

        if add_index:
            df.loc[pmd_mask, add_index_name] = i

    try:
        filtered_df = df.loc[:, prev_mask]
    except Exception:
        filtered_df = df.loc[prev_mask, :]

    return filtered_df
