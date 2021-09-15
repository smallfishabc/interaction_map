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


def interactionmap_pairwise(name,proteinpath, psi, residue):
    print(proteinpath, psi, residue)
    os.chdir(proteinpath)
    [seq, length] = path.readsequence()
    workpath = path.subdir(proteinpath, psi, residue)
    contact=contactmapgeneration.generate_contactmap(name,readfromfile=0)
    location,value =chunk.mapping_chunk(name,7,contact.contact,contact.pair)
    interaction, raw_value = nl.normalization(value,location,a1=1.29,b1=-1.22)
    #contMean, pairs = average_contmean(trajlist)
    #testchunk=chunk.Chunk(15,10,'15_10',5,contact.contact,contact.pair)
    # testchunklist=normalization.chunklist(contMean,pairs,10,5)
    # interaction = normalization.chunk_normalization(testchunklist)
    # interaction, raw_value = nl.normalization(contact.contact, contact.pair)
    # np.savetxt("interaction.csv", interaction, delimiter=",")
    # np.savetxt("raw_value.csv", raw_value, delimiter=",")
    # print(interaction)
    cm.interaction_map(seq, length, interaction, raw_value, location, 'contact_S_0', value)

