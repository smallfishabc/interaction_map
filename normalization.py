# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:10:07 2021

@author: ShaharGroup-fyu
"""
# I will intergrate the normalization process into the contact map class and make Chunk class a subclass of Contact map
import numpy as np




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
def chunk_normalization(targetmap, pairs,a1=1.29,b1=-1.22):
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
        if value > 1 and i > 0.001:
            interaction[index] = 2
        elif value > 0.5 and i > 0.001:
            interaction[index] = 1
        elif value < -1 and i > 0.001:
            interaction[index] = -2
        elif value < -0.5 and i > 0.001:
            interaction[index] = -1
        print(pairs[index], value, i, adjustment)
    return interaction, raw_value