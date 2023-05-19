# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:10:07 2021

@author: ShaharGroup-fyu
"""
# I will intergrate the normalization process into the contact map class and make Chunk class a subclass of Contact map
import numpy as np
# fitting function of standard contact probability map.
# The fitting_function or the calculation of the fitting should be computed in a seperate file
def fitting_function(df, a, b):
    """

    :param df: interaction data
    :type df: dataframe
    :return: int, fitting result
    """
    return a * df['distance'] ** b


# Use standard contact probability curve to calculate interaction strength
def normalization(target_map, a1=13.12, b1=-2.32, inter_cutoff=(2, 1, -1, -2), value_list=(2, 1, 0, -1, -2)):
    """

    :param target_map: interaction data
    :type target_map: dataframe
    :param a1:
    :param b1:
    :param inter_cutoff: interaction strength cutoff
    :type inter_cutoff: list
    :param value_list: interaction strength (corresponding to the inter_cutoff)
    :type value_list: list
    :return:
    """
    # Create an empty array for calculating interaction.
    # (Interaction is the binary value representing the interaction type
    # Relative_Strength is the ratio between pairwise contact probability and standard curve)
    print(target_map.columns)
    target_map['gs_standard'] = target_map.apply(fitting_function, axis=1 ,args=(a1, b1))
    target_map['relative_strength'] = np.where(target_map['cont_prob'] == 0, 0,
                                               np.log(target_map['cont_prob'] / target_map['gs_standard']))
    col = 'relative_strength'
    choice_list = value_list
    condition = [target_map[col] >= inter_cutoff[0],
                 (target_map[col] >= inter_cutoff[1]) & (target_map[col] < inter_cutoff[0]),
                 (target_map[col] < inter_cutoff[1]) & (target_map[col] > inter_cutoff[2]),
                 (target_map[col] <= inter_cutoff[2]) & (target_map[col] > inter_cutoff[3]),
                 (target_map[col] <= inter_cutoff[3])
                 ]
    target_map['plot_value'] = np.select(condition, choice_list, default=0)
    return target_map
