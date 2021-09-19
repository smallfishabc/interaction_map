# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os

import numpy as np

import contactmapgeneration
import contactmapgeneration as cg
import contactmapplot as cm
import normalization as nl
import path
import chunk
import matplotlib.pyplot as plt

def interactionmap_chunk(name,proteinpath, psi, residue):
    print(proteinpath, psi, residue)
    os.chdir(proteinpath)
    [seq, length] = path.readsequence()
    graphg = cm.create_network(seq, length)
    (pos, layout) = cm.create_position(graphg)
    (negacharged, posicharged, aromatic) = cm.seq_color(seq)
    (fig, ax) = cm.create_color_coding(seq, graphg, pos, negacharged, posicharged, aromatic)
    workpath = path.subdir(proteinpath, psi, residue)
    contact=contactmapgeneration.generate_contactmap(name,readfromfile=0)
    location,value=chunk.full_mapping_chunk(name,contact.contact,contact.pair)
    for index,i in enumerate(location):
        interaction,raw_value=nl.chunk_normalization(value[index],i)
        cm.interaction_map_chunk(seq, length, interaction, raw_value, i, layout , ax, 'contact_S_0', value[index])
    #    cm.interaction_map_chunk(seq, length, interaction, raw_value, i, layout , ax, 'contact_S_0', value[index])
    plt.show()

