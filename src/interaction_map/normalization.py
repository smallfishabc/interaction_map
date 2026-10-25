# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:10:07 2021

@author: ShaharGroup-fyu
"""
# I will intergrate the normalization process into the contact map class and make Chunk class a subclass of Contact map
import numpy as np

from . import profiles


# fitting function of standard contact probability map.
# The fitting_function or the calculation of the fitting should be computed in a seperate file
def fitting_function(df, a, b):
    """

    :param df: interaction data
    :type df: dataframe
    :return: int, fitting result
    """
    return a * df['distance'] ** b


# Reference contact probability for a sequence separation. Same curve as fitting_function,
# but it takes the separation directly so it can be used on an array or a single value.
def reference_probability(separation, a, b):
    """

    :param separation: sequence separation, j - i
    :type separation: int/array
    :return: reference contact probability
    """
    return a * np.asarray(separation, dtype=float) ** b


# Bin the interaction strength into the five levels the map draws.
# Split out of normalization so the same binning can be applied with different cutoffs
# without recomputing the contact map. The defaults are the ones normalization uses.
def classify(relative_strength, inter_cutoff=None, value_list=None):
    """
    :param relative_strength: ln(p_obs / p_ref), a series or array
    :param inter_cutoff: bin edges (strong att, weak att, weak rep, strong rep).
                         Defaults to the all-atom profile's.
    :param value_list: levels assigned to those bins
    :return: array of interaction levels
    """
    if inter_cutoff is None:
        inter_cutoff = profiles.DEFAULT_PROFILE.inter_cutoff
    if value_list is None:
        value_list = profiles.DEFAULT_PROFILE.value_list
    col = relative_strength
    condition = [col >= inter_cutoff[0],
                 (col >= inter_cutoff[1]) & (col < inter_cutoff[0]),
                 (col < inter_cutoff[1]) & (col > inter_cutoff[2]),
                 (col <= inter_cutoff[2]) & (col > inter_cutoff[3]),
                 (col <= inter_cutoff[3])
                 ]
    return np.select(condition, value_list, default=0)


# Use standard contact probability curve to calculate interaction strength
def normalization(target_map, profile=None, a1=None, b1=None, inter_cutoff=None, value_list=None,
                  prob_floor=None):
    """

    :param target_map: interaction data
    :type target_map: dataframe
    :param profile: analysis profile supplying the reference curve and the bin edges
    :type profile: str/profiles.Profile/None
    :param a1: override the profile's reference curve prefactor
    :param b1: override the profile's reference curve exponent
    :param inter_cutoff: interaction strength cutoff
    :type inter_cutoff: list
    :param value_list: interaction strength (corresponding to the inter_cutoff)
    :type value_list: list
    :param prob_floor: contact probabilities at or below this value are forced to level 0
    :type prob_floor: float
    :return:
    """
    profile = profiles.get_profile(profile)
    a1 = profile.reference_a if a1 is None else a1
    b1 = profile.reference_b if b1 is None else b1
    inter_cutoff = profile.inter_cutoff if inter_cutoff is None else inter_cutoff
    value_list = profile.value_list if value_list is None else value_list
    prob_floor = profile.prob_floor if prob_floor is None else prob_floor

    # Create an empty array for calculating interaction.
    # (Interaction is the binary value representing the interaction type
    # Relative_Strength is the ratio between pairwise contact probability and standard curve)
    # Left as the original row-wise apply. A vectorised expression gives bit-identical
    # results here, so there is nothing to gain by changing the form of a calibrated step.
    target_map['gs_standard'] = target_map.apply(fitting_function, axis=1, args=(a1, b1))
    target_map['relative_strength'] = np.where(target_map['cont_prob'] == 0, 0,
                                               np.log(target_map['cont_prob'] / target_map['gs_standard']))
    target_map['plot_value'] = classify(target_map['relative_strength'], inter_cutoff, value_list)
    # The published method applies a noise cutoff to the contact probability as well as to the
    # interaction strength. prob_floor = 0 reproduces the branches this was combined from.
    if prob_floor > 0:
        target_map.loc[target_map['cont_prob'] <= prob_floor, 'plot_value'] = 0
    return target_map
