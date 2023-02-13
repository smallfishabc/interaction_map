# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:10:07 2021

@author: ShaharGroup-fyu
"""
# I will intergrate the normalization process into the contact map class and make Chunk class a subclass of Contact map
import numpy as np


# Define a Chunk class for identifying residue chunk and motif.
# We will covert pairwise contact probability from
class Chunk:
    def __init__(self, residue, distance, name, chunk_size, contactmap, pairs):
        # name of the chunk class
        self.name = name
        # which residue is the center residue of the chunk
        self.residue = residue
        # Distance between the Chunk pairs eg: two chunk 1-5 and 10-15 will have a distance of 10
        self.distance = distance
        # How many residues in a chunk
        self.chunk_size = chunk_size
        # Define the center pair
        self.pair = self.search_pairs()
        # Define the center residue on the left
        self.left_residue = self.pair[0]
        # Define the center residue on the right
        self.right_residue = self.pair[1]
        # Define the chunk on the left
        self.left_residue_list = self.generate_residue_list('left')
        # Define the chunk on the right
        self.right_residue_list = self.generate_residue_list('right')
        # We calculate the chunk interaction based on individual pairwise interaction. Therefore we need to search for
        # residue pairs
        # Maybe search the pair between all residue in the chunk
        self.pair_list = self.pair_list()
        # Calculate the interaction in the pair list
        self.interact_list = self.calc_interact(contactmap, pairs)
        # Summarize the interaction in the chunk pair
        self.interact = self.sum_interact()

    # May be modified in nex t update
    def search_pairs(self):
        return [self.residue, self.residue + self.distance]

    # Define the entire chunk based on the center residue
    def generate_residue_list(self, position):
        # Create a empty residue list
        residue_list = []
        # Select the center residue based on the position
        if position == 'left':
            i = self.left_residue
        if position == 'right':
            i = self.right_residue
        # Calculate the half size of the chunk
        size = (self.chunk_size - 1) / 2
        k = 0
        # Append correct residue index to the list
        while k < self.chunk_size:
            residue_list.append(i - size + k)
            k += 1
        return residue_list

    # Will be updated to include all pairwise interaction between the two chunk
    def pair_list(self):
        plist = []
        for index, i in enumerate(self.left_residue_list):
            pair = [i, self.right_residue_list[index]]
            plist.append(pair)
        return plist

    # Retrieve the contact probability from the contact map
    def calc_interact(self, contactmap, pairs):
        interaction_list = []
        temppairs = pairs.tolist()
        for i in self.pair_list:
            indexpairs = temppairs.index(i)
            print(i)
            print(indexpairs)
            interaction_value = contactmap[indexpairs]
            print(interaction_value)
            interaction_list.append(interaction_value)
        return interaction_list

    # Calculate the average contact probability
    def sum_interact(self):
        average = np.log(np.average(self.interact_list))
        return average


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
def normalization(target_map, a1=1.64, b1=-1.32, inter_cutoff=(1.5, 0.5, -1, -2), value_list=(2, 1, 0, -1, -2)):
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
