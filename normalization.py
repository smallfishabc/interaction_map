# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:10:07 2021

@author: ShaharGroup-fyu
"""
# I will intergrate the normalization process into the contact map class and make Chunk class a subclass of Contact map
import numpy as np


class Chunk:
    def __init__(self, residue, distance, name, chunk_size, contactmap, pairs):
        self.name = name
        self.residue = residue
        self.distance = distance
        self.chunk_size = chunk_size
        self.pair = self.search_pairs()
        self.left_residue = self.pair[0]
        self.right_residue = self.pair[1]
        self.left_residue_list = self.generate_residue_list('left')
        self.right_residue_list = self.generate_residue_list('right')
        self.pair_list = self.pair_list()
        self.interact_list = self.calc_interact(contactmap, pairs)
        self.interact = self.sum_interact()

    def search_pairs(self):
        return [self.residue, self.residue + self.distance]

    def generate_residue_list(self, position):
        residue_list = []
        if position is 'left':
            i = self.left_residue
        if position is 'right':
            i = self.right_residue
        k = 0
        size = (self.chunk_size - 1) / 2
        while k < self.chunk_size:
            residue_list.append(i - size + k)
            k += 1
        return residue_list

    def pair_list(self):
        plist = []
        for index, i in enumerate(self.left_residue_list):
            pair = [i, self.right_residue_list[index]]
            plist.append(pair)
        return plist

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

    def sum_interact(self):
        average = np.log(np.average(self.interact_list))
        return average


def fitting_function(x, a, b):
    return a * x ** b


def normalization(targetmap, pairs, a1=1.64, b1=-1.32):
    interaction = np.zeros(pairs.shape[0])
    raw_value = np.zeros(pairs.shape[0])
    for index, i in enumerate(targetmap):
        distance = (pairs[index][1] - pairs[index][0])
        adjustment = fitting_function(distance, a1, b1)
        value = 0
        if i != 0:
            value = np.log(i / adjustment)

        else:
            value = 0
        raw_value[index] = value
        if value > 1.5 and i > 0.001:
            interaction[index] = 2
        elif value > 0.5 and i > 0.001:
            interaction[index] = 1
        elif value < -1.5 and i > 0.001:
            interaction[index] = -2
        elif value < -1 and i > 0.001:
            interaction[index] = -1
        print(pairs[index], value, i, adjustment)
    return interaction, raw_value
