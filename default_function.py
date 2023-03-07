# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os
import numpy as np
import contact_map_generation
import interaction_plot as cm
import normalization as nl
import readpath
# This one is designed for multiple trajctory.
# It will be removed for the final publication
# The input of this funcion is the trajectory home path
# It will return a trajectory full path and a protein sequence
def multi_traj_pre(name, path, psi, residue):
    _traj_path= os.path.join(path,residue,'_'.join(['S',psi]))
    _seq = readpath.readsequence(path)
    return _traj_path, _seq


# Main function to generate interaction map and calculate the interaction strength
# condition can be modified to kwrags later
# How to set up several default parameters in python
def interaction_map_pairwise(name, traj_path, sequence, output_dir, repeat, read_from_file=False):
    """
    :param name: protein name
    :type name: str
    :param traj_path: simulation trajectory file location
    :type name: str
    :param read_from_file: read or recalculate the contact map
    :type read_from_file: bool
    :param condition: additional simulation information
    :type list
    :return:
    """
    # Print out the target protein and solution condition
    print(traj_path, read_from_file)
    # Change dir to the target directory to read data and save analysis results
    os.chdir(traj_path)
    # Read sequence from seq.fasta
    seq=sequence
    length=len(seq)
    if not read_from_file:
        contact = contactmapgeneration.generate_contact(name, repeats=repeat)
        # Tracking the running status for debugging
        # Calculate interaction type based on ideal polymer model
        interaction = nl.normalization(contact.contact)
        # Save the interaction map's raw value into the csv file for backup
        interaction.to_csv(os.path.join(output_dir,"interaction_1202.csv"))
    # If we need to recalculate a existing IDP
    else:
        os.chdir(output_dir)
        # load data from saved contact map data
        contact = contactmapgeneration.generate_contact(name, read_from_file)
        # load the interaction map data from csv file
        interaction = nl.normalization(contact.contact)

    # Calculate the interaction strength based on the raw value
    cm.interaction_map(seq, length, interaction, name)
    # Return the overall interaction strength
    return interaction
'''
# Generating contact map is the speed limit step for this script. This function is designed to generate and store
# the contact map into file for further analysis. (Under development)
def interactionmap_pairwise_pre(name, traj_path, sequence, output_dir, repeat, read_from_file=False):
    # Print out the target protein and solution condition
    print(traj_path, read_from_file)
    # Change dir to the target directory to read data and save analysis results
    os.chdir(traj_path)
    # Read sequence from seq.fasta
    seq=sequence
    length=len(seq)
    # Generate contact map and save in the text file
    contact = contactmapgeneration.generate_contact(name, repeats=repeat)
    # Free up memory space
    del contact
'''
'''
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
'''