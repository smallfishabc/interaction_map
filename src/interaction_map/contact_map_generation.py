# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:23:10 2021

@author: ShaharGroup-fyu
"""
# Optional module

import os

import mdtraj as md
import numpy as np
import pandas as pd

from . import profiles
from . import schema


# One letter code for a residue. mdtraj knows the standard amino acids; coarse-grained
# topologies name their beads with the one letter code already, and anything else, such as
# an ACE or NME cap, comes back as X.
def residue_letter(residue):
    code = residue.code
    if code:
        return code
    if len(residue.name) == 1:
        return residue.name
    return 'X'


# Trajectory file names for a run. An integer means the numbered repeats of the simulation
# database; a string means one file; anything else is taken as an explicit list.
def trajectory_file_list(xtc_input):
    if isinstance(xtc_input, int):
        return ["__traj_{}.xtc".format(i) for i in range(xtc_input)]
    if isinstance(xtc_input, str):
        return [xtc_input]
    return list(xtc_input)


# Frames in each trajectory file, read from the file itself without loading coordinates.
# Returns an empty list if any file cannot be inspected, which simply means no spread is
# reported.
def replica_frame_counts(xtc_input):
    counts = []
    for path in trajectory_file_list(xtc_input):
        try:
            with md.open(path) as handle:
                counts.append(len(handle))
        except Exception:
            return []
    return counts


# We need to first load a trajectory file to use Contact map class
class ContactProbData:
    """
    Containing contact probability and the pairwise information. Will be saved to a dataframe.(Can be saved to csv)
    Traj is the trajectory file pass to the class
    Read from file can save time and memory when dealing with large protein sequences

    Attributes:
        name  The protein and the contact map name
        cutoff The cutoff distance of contact_map unit:nm
        profile The analysis profile (all_atom / cg) deciding how pairs are selected
        contact 4 columns Dataframes with contact probability and residue pair information
    """

    def __init__(self, name, cutoff, traj, read_from_file=0, profile=None, replica_frames=None,
                 output_dir=None):
        # Contact name include different information
        self.name = name
        # Cutoff distance for determining the contact
        self.cutoff = cutoff
        # Analysis profile. Decides pair selection and residue numbering.
        self.profile = profiles.get_profile(profile)
        # Frames contributed by each repeat, in the order they were concatenated. Used to
        # report the spread between repeats; the averages themselves are unaffected.
        self.replica_frames = list(replica_frames) if replica_frames else []
        # Where the contact csv is written and read. None means the working directory,
        # which is where it has always gone.
        self.output_dir = output_dir

        if read_from_file:
            # Read contact map from saved file
            self.read_contact()
        else:
            # Calculate contact probability and residue pairs
            self.contact = self.compute_contact(traj)
            self.update_csv()

    # Define a function to calculate contact probability
    def compute_contact(self, traj):
        # Calculate contact map using the MDtraj library. Cont stored the distance between two residues.
        if self.profile.pair_selection == 'mdtraj_ca':
            # The scheme 'CA' and ignore_nonprotein removed ions and ACE NME caps.
            # Note that contacts='all' already drops every pair with |i - j| < 3.
            [cont, pairs] = md.compute_contacts(traj, contacts='all', scheme='CA', ignore_nonprotein=True)
        elif self.profile.pair_selection == 'all_pairs':
            # Coarse-grained models carry one bead per residue, so atom indices and residue
            # indices coincide and select_pairs can be used to keep every pair, including the
            # short separations that contacts='all' would discard.
            indices = traj.top.select_pairs('all', 'all')
            [cont, pairs] = md.compute_contacts(traj, indices)
        else:
            raise ValueError('Unknown pair_selection {!r}'.format(self.profile.pair_selection))
        # Determine whether the distance is smaller than the cutoff.
        contact = (cont < self.cutoff)
        # Calculate the average contact frequency of each residue pair
        cont_mean = contact.mean(0)
        # Shift the raw mdtraj index so r_1 / r_2 are 1-based sequence positions.
        offset = self.profile.residue_index_offset
        cont_df = pd.DataFrame({'r_1': pairs[:, 0] + offset, 'r_2': pairs[:, 1] + offset, 'cont_prob': cont_mean[:]})
        # int64 explicitly, not 'int': the latter is the C long, so it is 32 bit on Windows
        # and 64 bit elsewhere, and a CSV round trip always comes back 64 bit.
        cont_df['r_1'] = cont_df['r_1'].astype('int64')
        cont_df['r_2'] = cont_df['r_2'].astype('int64')
        cont_df['distance'] = cont_df['r_2'] - cont_df['r_1']
        # Residue identity straight from the topology, so a pair can be read without
        # consulting the sequence file separately.
        residues = list(traj.top.residues)
        cont_df['r_1_aa'] = [residue_letter(residues[i]) for i in pairs[:, 0]]
        cont_df['r_2_aa'] = [residue_letter(residues[i]) for i in pairs[:, 1]]
        cont_df['r_1_chain'] = [residues[i].chain.index for i in pairs[:, 0]]
        cont_df['r_2_chain'] = [residues[i].chain.index for i in pairs[:, 1]]
        # Spread between repeats. cont_prob above is the average over every frame and is
        # left exactly as it was; these columns only describe it.
        cont_df['n_frames'] = contact.shape[0]
        self._add_replica_spread(cont_df, contact)
        # Free up the memory
        del contact
        del pairs
        return cont_df

    # Standard error of the per repeat contact probabilities.
    # The repeats were concatenated at load time, so they are recovered by slicing the
    # frame axis back up with the counts recorded for each file. Nothing here touches
    # cont_prob.
    def _add_replica_spread(self, cont_df, contact):
        counts = self.replica_frames
        if len(counts) < 2 or sum(counts) != contact.shape[0]:
            # One repeat, or frame counts that do not account for the trajectory. No
            # spread can be reported, and saying so is better than inventing one.
            cont_df['n_rep'] = len(counts) if counts else 1
            cont_df['p_sem'] = np.nan
            return
        edges = np.cumsum([0] + list(counts))
        per_replica = np.vstack([contact[a:b].mean(0) for a, b in zip(edges[:-1], edges[1:])])
        cont_df['n_rep'] = len(counts)
        cont_df['p_sem'] = per_replica.std(axis=0, ddof=1) / np.sqrt(len(counts))

    # The file name keeps its date stamp on purpose: it identifies the version of the
    # analysis that produced the file.
    def contact_path(self):
        filename = '_'.join([str(self.name), str(self.cutoff), "contact_df_1201.csv"])
        return os.path.join(self.output_dir, filename) if self.output_dir else filename

    # Define a function to read the average contact probability dataframe from file
    def read_contact(self):
        self.contact = pd.read_csv(self.contact_path(), index_col=0)
        # A file written by an older version, or a truncated one, shows up here rather
        # than as a confusing result several steps later. Warning only.
        schema.warn_if_invalid(self.contact, schema.CONTACT_TABLE)

    # Define a function to save the average contact probability dataframe into file
    def save_contact(self):
        self.contact.to_csv(self.contact_path())

    def update_csv(self):
        self.save_contact()


# Normally, we input repeats number into this function to get the full trajectory file of the contact map
# If switch == 1, we can input single traj file.
# We can also input a list to traj_selection indicating the index of individual repeats to import part of
# the entire simulation.
# We will also create stdoutput function to calculate the standard deviation of the traj by only return a list of traj.
def check_residue_numbering(contact_df, sequence):
    """
    Warn if r_1 / r_2 do not line up with the sequence that will be drawn underneath them.

    The all-atom topologies carry an ACE cap at residue index 0, so mdtraj's own indices
    already count from 1 for the first real residue. An uncapped topology would come out one
    residue low and every arc would be drawn against the wrong residue. Diagnostic only,
    nothing is altered.

    :param contact_df: contact probability dataframe
    :param sequence: protein sequence
    :return: True when the numbering spans 1..len(sequence)
    """
    low = int(contact_df['r_1'].min())
    high = int(contact_df['r_2'].max())
    if low == 1 and high == len(sequence):
        return True
    print('WARNING: residue numbering runs {0}..{1} but the sequence has {2} residues. '
          'Arcs may be drawn against the wrong residues. Check the topology for terminal '
          'caps and the profile residue_index_offset.'.format(low, high, len(sequence)))
    return False


def topology_sequence(traj):
    """
    One letter sequence of the standard residues in a topology, in order.

    Terminal caps and ions have no one letter code and are skipped, so for a capped
    all-atom topology this returns the protein sequence itself. Coarse-grained topologies
    name their beads with one letter codes that mdtraj does not recognise as residues, so
    this comes back empty for them.

    :param traj: mdtraj trajectory or topology
    :return: str
    """
    top = traj.top if hasattr(traj, 'top') else traj
    return ''.join(r.code for r in top.residues if r.code is not None)


def check_sequence_against_topology(traj, sequence):
    """
    Warn when the sequence about to be drawn is not the one that was simulated.

    Returns True when they agree, False when they do not, and None when the topology
    carries no readable residue codes, which is the normal case for coarse-grained models.
    Diagnostic only, nothing is altered.

    :param traj: mdtraj trajectory
    :param sequence: sequence read from seq.fasta or seq.txt
    :return: True / False / None
    """
    from_topology = topology_sequence(traj)
    if not from_topology:
        return None
    if from_topology == sequence:
        return True
    if len(from_topology) != len(sequence):
        print('WARNING: the topology holds {0} residues but the sequence file has {1}. '
              'The map would be drawn against the wrong residues.'.format(
                  len(from_topology), len(sequence)))
        return False
    for position, (a, b) in enumerate(zip(from_topology, sequence), start=1):
        if a != b:
            print('WARNING: sequence file and topology disagree at residue {0}: '
                  'file has {1!r}, topology has {2!r}.'.format(position, b, a))
            return False
    return False


def load_traj_protein(pdb_top='__START_0.pdb', xtc_input='__traj_0.xtc', profile=None):
    '''
    :param pdb_top: pdb file name as topology file
    :type: string
    :param xtc_input: xtc name list o/ repeat number (will auto generate name list)
    :type: list/int
    :param profile: analysis profile deciding whether the topology is sliced to 'protein'
    :return: mdtraj.trajectory, protein-only trajectory, trajectory framenumber
    '''
    profile = profiles.get_profile(profile)
    list_xtc = trajectory_file_list(xtc_input)
    t = load_xtc(list_xtc, pdb_top)
    if profile.slice_protein:
        # Select the protein from the pdb file
        u = t.top.select('protein')
        r = t.atom_slice(u)
    else:
        # Coarse-grained topologies carry no caps or ions, and the beads are not recognised
        # as protein by mdtraj, so the selection is skipped.
        r = t
    # Get frame number of the trajectory file.
    return t, r, r.n_frames


# input a xtc list to load 1 or more xtc files
def load_xtc(xtc, pdb):
    """
    Load xtc/xtc file list with correct pdb topology.
    :param xtc: a list contains xtc file names
    :type xtc: list/string
    :param pdb: pdb file name
    :type pdb: string
    :return: mdtraj.Trajectory
    """
    return md.load(xtc, top=pdb)


# What if we input a single xtc file, this need to be discuss and changed
def generate_contact(protein_name, pdb_top='__START_0.pdb', xtc_input=5, cutoff=None, read_from_file=0,
                     profile=None, sequence=None, output_dir=None):
    """
    :param protein_name: input protein_name/ contact_map_name
    :type protein_name: str
    :param cutoff: cutoff of contact_map unit:nm. Taken from the profile when not given.
    :type cutoff: float
    :param read_from_file: read or recalculate the contact map
    :type read_from_file: bool
    :param profile: analysis profile (all_atom / cg)
    :param sequence: when given, the topology is checked against it and a warning printed
                     on a mismatch. Nothing is altered either way.
    :param output_dir: where the contact csv is written. None means the working directory.
    :return: a Contactmap object with the corresponding protein contact map
    """
    profile = profiles.get_profile(profile)
    cutoff = profile.cutoff if cutoff is None else cutoff

    if read_from_file:
        contact = ContactProbData(protein_name, cutoff, traj=None, read_from_file=1, profile=profile,
                                  output_dir=output_dir)
        return contact

    traj, sliced_traj, frames = load_traj_protein(pdb_top, xtc_input, profile=profile)
    if sequence is not None:
        check_sequence_against_topology(sliced_traj, sequence)
    contact = ContactProbData(protein_name, cutoff, sliced_traj, profile=profile,
                              replica_frames=replica_frame_counts(xtc_input),
                              output_dir=output_dir)
    return contact
