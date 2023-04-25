# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os
import contact_map_generation as cg
import interaction_plot as ip
import normalization as nl
import readpath
import pandas as pd

# This one is designed for multiple trajectory.
# It will be removed for the final publication
# The input of this function is the trajectory home path
# It will return a trajectory full path and a protein sequence
def multi_traj_pre(name, path, psi, residue):
    _traj_path = os.path.join(path, residue, '_'.join(['S', psi]))
    _seq = readpath.readsequence(path)
    return _traj_path, _seq


# Main function to generate interaction map and calculate the interaction strength
# condition can be modified to kwargs later
# How to set up several default parameters in python
def interaction_map_pairwise(name, traj_path, sequence, output_dir, pdb_top='__START_0.pdb', xtc_input=5, read_from_file=False):
    """
    :param name: protein name
    :type name: str
    :param traj_path: simulation trajectory file location
    :type name: str
    :param read_from_file: read or recalculate the contact map
    :type read_from_file: bool
    :return:
    """
    # Print out the target protein and solution condition
    print(traj_path, read_from_file)
    # Change dir to the target directory to read data and save analysis results
    os.chdir(traj_path)
    # Read sequence from seq.fasta
    seq = sequence
    length = len(seq)
    if not read_from_file:
        contact = cg.generate_contact(name, pdb_top, xtc_input)
        # Tracking the running status for debugging
        # Calculate interaction type based on ideal polymer model
        interaction = nl.normalization(contact.contact)
        # Save the interaction map's raw value into the csv file for backup
        interaction.to_csv(os.path.join(output_dir, "interaction_0317.csv"))
    # If we need to recalculate an existing IDP
    else:
        os.chdir(output_dir)
        # load data from saved contact map data
        contact = cg.generate_contact(name, read_from_file=True)
        # load the interaction map data from csv file
        interaction = nl.normalization(contact.contact)

    # Calculate the interaction strength based on the raw value
    ip.interaction_map(seq, length, interaction, name)
    # Return the overall interaction strength
    return interaction


