import numpy as np
import normalization as nl


# In this module, we define 2 function to calculate the 'true' interaction strength based on normalization
# and the plotting result. Here the strength should correlate with the IDP solution sensitivity.

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
    factor = 1
    return factor


# Interaction strength is proportional to the distance factor
def calculate_individual_strength_none(distance_factor, contact_none, standard_nono=0):
    strength = distance_factor
    return strength


# Interaction strength is proportional to the distance factor and the contact probability ratio.
def calculate_individual_strength_contact_ratio(distance_factor, contact_ratio, standard_ratio=0):
    # Standard_ratio should be 1 if we want to subtract ideal polymer contact probability from calculation
    strength = distance_factor * (contact_ratio - standard_ratio)
    return strength


# Interaction strength is proportional the the distance factor and the contact probability value
def calculate_individual_strength_contact_probability(distance_factor, contact_value, standard_value=0):
    # Standard_value should be calculated if we want to subtract ideal polymer contact probability from calculation
    strength = distance_factor * (contact_value - standard_value)
    return strength

# This function is designed to directly calculate the overall interaction strength based on data in the interaction
# plot function.
def calculate_overall_strength(interaction, intertype, function=calculate_individual_strength_none, factor=distance_factor_none):
    # Create an empty numpy array for interaction strength calculation
    full_strength = np.zeros(len(interaction))
    # Retrieve contact probability and calculate interaction strength for every residue pairs.
    for i in range(len(interaction)):
        # Calculate the distance between two residue
        distance = edge[1] - edge[0]
        # Determine which function and algorithm to be used based on the input
        if function is calculate_individual_strength_contact_ratio:
            # Contact probability ratio. Standard ratio is calculated for prepare the difference between contact
            # probability and standard curve.
            value = raw_value[index]
            if value >= 0:
                standard_ratio = 1
            elif value < 0:
                standard_ratio = -1
        else:
            # Contact probability value.Standard ratio is calculated for prepare the difference between contact
            # probability and standard curve.
            value = contactmap[index]
            standard_value = nl.fitting_function(distance, a=1.64, b=-1.32)
        # We retrieve the other values just for debugging.
        interact = interaction[index]
        contact = contactmap[index]
        # Based on the distance between residue pair, we will calculate the distance_factor
        distance_factor = factor(distance)
        # Based on the distance factor and the contact probability value, we can calculate the interaction strength
        strength = function(distance_factor, value)
        # Store the interaction strength into corresponding array
        full_strength[index] = strength
    # Overall strength is the sum over of the each individual interaction.
    overall_strength = full_strength.sum()
    return full_strength, overall_strength


# Define a function to calculate multiple function and distance_factor
def calculate_overall_strength_multiple(pairs, raw_value, interaction, contactmap):
    # We defined 3 distance factor based on different assumption
    factor_list = [distance_factor_none, distance_factor_A, distance_factor_B]
    # We defined 3 contact probability factor based on different assumption
    function_list = [calculate_individual_strength_none, calculate_individual_strength_contact_ratio,
                     calculate_individual_strength_contact_probability]
    # Create an empty list to store the result. Will change to 2-d array in next update
    list_a = []
    # Calculate all assumption and append them into a list
    for i in factor_list:
        for j in function_list:
            full_strength, overall_strength = calculate_overall_strength(pairs, raw_value, interaction, contactmap, j,
                                                                         i)
            list_a.append(overall_strength)
    return list_a



