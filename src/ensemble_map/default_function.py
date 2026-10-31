# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os
from . import contact_map_generation as cg
from . import interaction_plot as ip
from . import normalization as nl
from . import profiles
from . import readpath
from . import search_chunk as sc
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
def interaction_map_pairwise(name, traj_path, sequence, output_dir, pdb_top='__START_0.pdb', xtc_input=5,
                             read_from_file=False, profile=None, show=True,
                             interaction_csv='interaction_0317.csv'):
    """
    :param name: protein name
    :type name: str
    :param traj_path: simulation trajectory file location
    :type name: str
    :param read_from_file: read or recalculate the contact map
    :type read_from_file: bool
    :param profile: analysis profile, 'all_atom' or 'cg'
    :param show: display the figure. Turn this off for batch runs.
    :return:
    """
    profile = profiles.get_profile(profile)
    # Print out the target protein and solution condition
    print(traj_path, profile.name, read_from_file)
    # Change dir to the target directory to read data and save analysis results
    os.chdir(traj_path)
    # Read sequence from seq.fasta
    seq = sequence
    length = len(seq)
    if not read_from_file:
        contact = cg.generate_contact(name, pdb_top, xtc_input, profile=profile, sequence=seq,
                                      output_dir=output_dir)
        # Tracking the running status for debugging
        # Calculate interaction type based on ideal polymer model
        interaction = nl.normalization(contact.contact, profile=profile)
        # Save the interaction map's raw value into the csv file for backup
        interaction.to_csv(os.path.join(output_dir, interaction_csv))
    # If we need to recalculate an existing IDP
    else:
        os.chdir(output_dir)
        # load data from saved contact map data
        contact = cg.generate_contact(name, read_from_file=True, profile=profile)
        # load the interaction map data from csv file
        interaction = nl.normalization(contact.contact, profile=profile)

    # Calculate the interaction strength based on the raw value.
    # The figure goes beside the tables. When output_dir is the trajectory folder, which
    # is the default, this is the same place it has always been written.
    figure_name = os.path.join(output_dir, name) if output_dir else name
    ip.interaction_map(seq, length, interaction, figure_name, profile=profile, show=show)
    # Return the overall interaction strength
    return interaction


# Conditional analysis. Split the ensemble on whether one interaction is formed, build a
# contact map for each side, and compare their ensemble properties.
#
# Ported from 1115_new_start. The frames are selected with the same 3x3 window rule as
# before, now expressed as a ChunkContact object so the window and the cutoff are visible
# and can be varied. As in the original, the contact maps and the observables are computed
# on the full trajectory, not the protein slice.
def search_interaction(name, r1, r2, traj_path, sequence, output_dir, pdb_top='__START_0.pdb',
                       xtc_input=5, read_from_file=False, profile=None, show=True,
                       window=1, chunk=True, cutoff=None):
    """
    :param r1: first residue of the target interaction, numbered as in the contact table
    :param r2: second residue of the target interaction
    :param window: half width of the residue window. 1 is the published 3x3 window.
    :param chunk: False narrows the selection to the single pair r1, r2
    :return: (interaction table with the contact, interaction table without it, ensemble summary)
    """
    profile = profiles.get_profile(profile)
    os.chdir(traj_path)
    traj, sliced_traj, frames = cg.load_traj_protein(pdb_top, xtc_input, profile=profile)

    selector = sc.ChunkContact(r1, r2, window=window, cutoff=cutoff, chunk=chunk, profile=profile)
    print('selecting frames with {0}'.format(selector))
    selected_frames = sc.frames_where(traj, selector)
    n_with = int(selected_frames.sum())
    print('  {0} of {1} frames contain the interaction ({2:.1%})'.format(
        n_with, len(selected_frames), n_with / float(len(selected_frames))))
    if n_with == 0 or n_with == len(selected_frames):
        print('  WARNING: one side of the split is empty, the comparison is meaningless')

    contact, contact_non = sc.new_interaction_map(name, traj, selected_frames, profile=profile,
                                                  output_dir=output_dir)
    interaction = nl.normalization(contact.contact, profile=profile)
    interaction_non = nl.normalization(contact_non.contact, profile=profile)

    figure_name = os.path.join(output_dir, name) if output_dir else name
    ip.interaction_map(sequence, len(sequence), interaction, figure_name, profile=profile, show=show)

    ensemble = sc.compare_ensemble(traj, sliced_traj, len(sequence),
                                   contact_select=selected_frames, profile=profile)
    ensemble_pd = pd.DataFrame(ensemble)
    suffix = '_'.join([str(r1), str(r2)]) + '.csv'
    ensemble_pd.to_csv(os.path.join(output_dir, '_'.join(['ensemble', suffix])), index=True)
    interaction.to_csv(os.path.join(output_dir, '_'.join(['with', 'interaction', suffix])), index=False)
    interaction_non.to_csv(os.path.join(output_dir, '_'.join(['without', 'interaction', suffix])), index=False)
    return interaction, interaction_non, ensemble_pd


# Compare ensemble properties with and without a chosen set of interactions.
#
# search_interaction above answers the question for one interaction and also rebuilds the
# contact maps. This one takes any number of interactions, combines them with 'any' or
# 'all', and reports only the ensemble properties, which is usually what the question is.
#
#   mode 'any'   frames where at least one of the interactions is formed
#   mode 'all'   frames where every one of them is formed at the same time
def compare_interactions(name, pairs, traj_path, sequence, output_dir, pdb_top='__START_0.pdb',
                         xtc_input=5, profile=None, mode='any', window=1, chunk=True,
                         cutoff=None, with_maps=False, show=False):
    """
    :param pairs: iterable of (r1, r2), numbered as in the contact table
    :param mode: 'any' or 'all'
    :param with_maps: also rebuild and save a contact map for each side
    :return: (dataframe of the comparison, boolean mask over frames)
    """
    profile = profiles.get_profile(profile)
    os.chdir(traj_path)
    traj, sliced_traj, frames = cg.load_traj_protein(pdb_top, xtc_input, profile=profile)

    selector = sc.interaction_set(pairs, mode=mode, window=window, chunk=chunk,
                                  cutoff=cutoff, profile=profile)
    print('selecting frames: {0}'.format(selector))
    selected = sc.frames_where(traj, selector)
    n_with = int(selected.sum())
    total = len(selected)
    print('  {0} of {1} frames ({2:.1%}) satisfy the {3!r} condition'.format(
        n_with, total, n_with / float(total), mode))
    if n_with == 0 or n_with == total:
        print('  WARNING: one side of the split is empty, the comparison is meaningless')
        return None, selected

    comparison = pd.DataFrame(
        sc.compare_ensemble(traj, sliced_traj, len(sequence), contact_select=selected,
                            profile=profile))
    # Record what was asked for alongside the answer, so the file can be read on its own.
    comparison.loc['mode'] = [mode, mode]
    comparison.loc['interactions'] = ['; '.join('{0}-{1}'.format(a, b) for a, b in pairs)] * 2

    label = selector.label()
    comparison.to_csv(os.path.join(output_dir, 'ensemble_' + label + '.csv'), index=True)

    if with_maps:
        contact, contact_non = sc.new_interaction_map(name, traj, selected, profile=profile,
                                                      output_dir=output_dir)
        interaction = nl.normalization(contact.contact, profile=profile)
        interaction_non = nl.normalization(contact_non.contact, profile=profile)
        interaction.to_csv(os.path.join(output_dir, 'with_interaction_' + label + '.csv'), index=False)
        interaction_non.to_csv(os.path.join(output_dir, 'without_interaction_' + label + '.csv'),
                               index=False)
        ip.interaction_map(sequence, len(sequence), interaction,
                           os.path.join(output_dir, name + '_' + label),
                           profile=profile, show=show)
    return comparison, selected
