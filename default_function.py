# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os
import json
import numpy as np

import contactmapgeneration
import contactmapgeneration as cg
import contactmapplot as cm
import normalization as nl
import readpath
import pandas as pd

def interactionmap_pairwise(name, proteinpath, psi, residue, read_from_file=0):
#    if read_from_file = 0:
#        generate=1
    #read_from_file=0
    print(proteinpath, psi, residue)
    os.chdir(proteinpath)
    [seq, length] = readpath.readsequence()
    workpath = readpath.subdir(proteinpath, psi, residue)
    # trajlist = [0, 1, 2, 3, 4]
    if read_from_file == 0:
        print('start1')
        try:
            contact = contactmapgeneration.generate_contactmap(name)
        except:
            contact = contactmapgeneration.generate_contactmap(name,repeats=3)
        print('end1')
        interaction, raw_value = nl.normalization(contact.contact, contact.pair)
        print('end11')
        np.savetxt("interaction.csv", interaction, delimiter=",")
        np.savetxt("raw_value.csv", raw_value, delimiter=",")
    else:
        contact = contactmapgeneration.generate_contactmap(name, read_from_file)
        interaction = np.loadtxt("interaction.csv", delimiter=",")
        raw_value = np.loadtxt("raw_value.csv", delimiter=",")
    print('start2')
    att1, att2, rep1, rep2 = cm.interaction_map_calc(seq, length, interaction, raw_value, contact.pair, 'contact_S_0',
                                                contact.contact)
    print('end2')
    return att1, att2, rep1, rep2

def interactionmap_pairwise_pre(name, proteinpath, psi, residue, read_from_file):
    print(proteinpath, psi, residue)
    os.chdir(proteinpath)
    [seq, length] = readpath.readsequence()
    workpath = readpath.subdir(proteinpath, psi, residue)
    contact = contactmapgeneration.generate_contactmap(name, read_from_file)
    del contact
    print(contact)

