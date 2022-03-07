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
        if position is 'left':
            i = self.left_residue
        if position is 'right':
            i = self.right_residue
        # Calculate the half size of the chunk
        size = (self.chunk_size - 1) / 2
        k=0
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

# Fitting function of standard contact probability map.
#  A simple exponential function
def fitting_function(x, a, b):
    return a * x ** b

# Use standard contact probability curve to calculate interaction strength (Key algorithm of my script)
def normalization(targetmap, pairs, a1=1.64, b1=-1.32):
    # Create an empty array for calculating interaction.(Interaction is the binary value representing the interaction type
    # raw_value is the ratio between pairwise contact probability and standard curve)
    interaction = np.zeros(pairs.shape[0])
    raw_value = np.zeros(pairs.shape[0])
    # Loop over the entire interaction map
    for index, i in enumerate(targetmap):
        # Calculate the distance of the residue pair.
        distance = (pairs[index][1] - pairs[index][0])
        # Calculate the contact probability of selected distance on the standard curve.
        adjustment = fitting_function(distance, a1, b1)
        # Create a variable to store the ratio
        value = 0
        # If the contact probability is not zero
        if i != 0:
            # Interaction = (Pcontact)/(Pcontact of GS)
            value = np.log(i / adjustment)
        else:
            # Interaction strength = 0
            value = 0
        # Store the ratio into the array
        raw_value[index] = value
        # Based on the ratio, determine the interaction type and store in the array
        if value > 2 and i > 0.001:
            interaction[index] = 2
        elif value > 1 and i > 0.001:
            interaction[index] = 1
        elif value < -2 and i > 0.001:
            interaction[index] = -2
        elif value < -1 and i > 0.001:
            interaction[index] = -1
        # print(pairs[index], value, i, adjustment)
    return interaction, raw_value
