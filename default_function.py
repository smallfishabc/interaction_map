# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os
import numpy as np
import contactmapgeneration
import contactmapplot as cm
import normalization as nl
import readpath
import contactmapcalculation as cc

# Main function to generate interaction map and calculate the interaction strength
def interactionmap_pairwise(name, proteinpath, psi, residue, read_from_file=0):
    # Print out the target protein and solution condition
    print(proteinpath, psi, residue,read_from_file)
    # Change dir to the target directory to read data and save analysis results
    os.chdir(proteinpath)
    # Read sequence from seq.fasta
    [seq, length] = readpath.readsequence()
    # TBD
    workpath = readpath.subdir(proteinpath, psi, residue)
    # Select desired trajectories from
    # trajlist = [0, 1, 2, 3, 4]
    # If we need to generate interaction map for a new IDP
    if read_from_file == 0:
        # Tracking the running status for debugging
        print('start1')
        # Currently we have two type of simulation. One has 5 repeats. The other one has 3 repeats.
        # This try function is designed for the compatibility.
        try:
            # Generate the contact map array
            contact = contactmapgeneration.generate_contactmap(name)
        except:
            contact = contactmapgeneration.generate_contactmap(name, repeats=3)
        # Tracking the running status for debugging
        print('end1')
        # Calculate interaction type based on ideal polymer model
        interaction, raw_value = nl.normalization(contact.contact, contact.pair)
        print('end11')
        # Save the interaction map's raw value into the csv file for backuo
        np.savetxt("interaction.csv", interaction, delimiter=",")
        np.savetxt("raw_value.csv", raw_value, delimiter=",")
    # If we need to recalculate a existing IDP
    else:
        # load data from saved contact map data
        contact = contactmapgeneration.generate_contactmap(name, read_from_file)
        # load the interaction map data from csv file
        interaction = np.loadtxt("interaction.csv", delimiter=",")
        raw_value = np.loadtxt("raw_value.csv", delimiter=",")
    print('start2')
    # Calculate the interaction strength based on the raw value
    att1, att2, rep1, rep2 = cc.interaction_map_calc(seq, length, interaction, raw_value, contact.pair, 'contact_S_0',
                                                     contact.contact)
    cm.interaction_map(seq, length, interaction, raw_value, contact.pair, 'figname', contact.contact)
    print('end2')
    # Return the overall interaction strength
    return att1, att2, rep1, rep2

# Generating contact map is the speed limit step for this script. This function is designed to generate and store
# the contact map into file for further analysis. (Under development)
def interactionmap_pairwise_pre(name, proteinpath, psi, residue, read_from_file):
    # Check the protein setup
    print(proteinpath, psi, residue)
    # Change directory to protein folder
    os.chdir(proteinpath)
    # Read sequence from the fast file
    [seq, length] = readpath.readsequence()
    # Change directory to actually working directory
    workpath = readpath.subdir(proteinpath, psi, residue)
    # Generate contact map and save in the text file
    contact = contactmapgeneration.generate_contactmap(name, read_from_file)
    # Free up memory space
    del contact
# Interaction map generation main function for single trajectory
def interactionmap_pairwise_single_traj(path,pdb_name,xtc_name,sequence):
    #Change working directory to target directory
    os.chdir(path)
    #Calculate the contact map
    contact = contactmapgeneration.generate_contactmap_single_traj(path,pdb_name,xtc_name)
    #Calculate the interaction strength
    interaction, raw_value = nl.normalization(contact.contact, contact.pair)
    length=len(sequence)
    print(contact)
    figname=pdb_name.split('.')[0]
    cm.interaction_map(sequence, length, interaction, raw_value, contact.pair, figname, contact.contact)