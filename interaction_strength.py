import numpy as np


# In this module, we define 2 function to calculate the 'true' interaction strength based on normalization
# and the plotting result. Here the strength should correlate with the IDP solution sensitivity.

# This function is designed to directly calculate the overall interaction strength based on data in the interaction
# plot function.
def calculate_overall_strength(pairs, raw_value, interaction, contactmap, function, factor):
    # Create an empty numpy array for interaction strength calculation
    full_strength = np.zeros(pairs.shape[0])
    for index, edge in enumerate(pairs):
        distance = edge[1] - edge[0]
        if function is calculate_individual_strength_contact_ratio:
            value = raw_value[index]
        else:
            value = contactmap[index]
        interact = interaction[index]
        contact = contactmap[index]
        distance_factor=factor(distance)
        strength = function(distance_factor, value)
        full_strength[index] = strength
    overall_strength = full_strength.sum()
    return full_strength, overall_strength


# Define a function to calculate multiple function and distance_factor
def calculate_overall_strength_multiple(pairs, raw_value, interaction, contactmap):
    factor_list=[distance_factor_none,distance_factor_A,distance_factor_B]
    function_list=[calculate_individual_strength_none,calculate_individual_strength_contact_ratio,calculate_individual_strength_contact_probability]
    list_a=[]
    for i in factor_list:
        for j in function_list:
            full_strength,overall_strength=calculate_overall_strength(pairs, raw_value, interaction, contactmap, j, i)
            list_a.append(overall_strength)
    return list_a

# Test of different possible function and distance_factor

# Long-range interaction may have a stronger influence on the IDP conformational ensemble. Therefore, we introduced
# a distance factor

# Distance factor is proportional to the residue number between the residue pair.
def distance_factor_A(distance):
    factor = distance
    return factor
# Distance factor is proportional to the square root of the residue number between the residue pair.
def distance_factor_B(distance):
    factor = np.sqrt(distance)
    return factor
# No distance factor
def distance_factor_none(distance):
    factor=1
    return factor

# Interaction strength is proportional to the distance factor
def calculate_individual_strength_none(distance_factor, contact_none, standard_nono=0):
    strength = distance_factor
    return strength

# Interaction strength is proportional to the distance factor and the contact probability ratio.
def calculate_individual_strength_contact_ratio(distance_factor, contact_ratio, standard_ratio=0):
    # Standard_ratio should be 1 if we want to subtract ideal polymer contact probability from calculation
    strength = distance_factor * (contact_ratio-standard_ratio)
    return strength

# Interaction strength is proportional the the distance factor and the contact probability value
def calculate_individual_strength_contact_probability(distance_factor, contact_value, standard_value=0):
    # Standard_value should be calculated if we want to subtract ideal polymer contact probability from calculation
    strength = distance_factor * (contact_value - standard_value)
    return strength
