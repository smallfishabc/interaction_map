import numpy as np
import pandas as pd


# In this module, we define 2 function to calculate the 'true' interaction strength based on normalization
# and the plotting result. Here the strength should correlate with the IDP solution sensitivity.

# This function is designed to directly calculate the overall interaction strength based on data in the interaction
# plot function.

def calculate_overall_strength(pairs, raw_value, interaction, contactmap, inter_type='attr', method='root'):
    full_strength = np.zeros(pairs.shape[0])
    for index, edge in enumerate(pairs):
        distance = edge[1] - edge[0]
        value = raw_value[index]
        interact = interaction[index]
        contact = contactmap[index]
        strength = calculate_individual_strength2(distance, value)
        full_strength[index] = strength
    overall_strength = full_strength.sum()
    return full_strength, overall_strength


# This function is designed to calculate the interaction strength of individual line during the execution of
# interaction_plot function.

def calculate_individual_strength(distance, raw_value):
    strength = np.sqrt(distance)
    return strength

def calculate_individual_strength2(distance, raw_value):
    strength = np.sqrt(distance)*raw_value
    return strength
