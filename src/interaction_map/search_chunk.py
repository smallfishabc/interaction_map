# This file is designed to search interactions after identified the strongest interaction.
# 1. Search what interaction happens simultaneously with the interaction
# 2. Compare of radius of gyration and helicity distribution with/without the interaction
#
# Ported from 1115_new_start. The calculations are unchanged; what is new is that the
# frame selector is an explicit object rather than a nested comprehension, the profile
# decides the cutoff and the residue numbering, and the comparison records how many frames
# went into each side of the split.
import mdtraj as md
import numpy as np
import pandas as pd
import statistics as st

from . import contact_map_generation as cg
from . import profiles


# A rule for deciding whether a frame contains the interaction of interest.
#
# The published choice is a 3x3 window: the interaction counts as present when ANY pair in
# {r1-1, r1, r1+1} x {r2-1, r2, r2+1} is within the cutoff. That is the single most
# consequential knob in the conditional analysis, so it is spelled out here rather than
# buried in the code that uses it.
class ChunkContact(object):
    """
    Attributes:
        r1, r2  residue positions, numbered the way the contact tables are, counting from 1
        window  how far either side of r1 and r2 to look. 1 gives the 3x3 window.
        cutoff  contact distance, nm
        chunk   False collapses the window to the single pair r1, r2
    """

    def __init__(self, r1, r2, window=1, cutoff=None, chunk=True, profile=None):
        self.profile = profiles.get_profile(profile)
        self.r1 = int(r1)
        self.r2 = int(r2)
        self.window = int(window)
        self.cutoff = self.profile.cutoff if cutoff is None else cutoff
        self.chunk = chunk

    # Residue index pairs as mdtraj numbers them, which is the stored position minus the
    # profile's offset.
    def pairs(self):
        offset = self.profile.residue_index_offset
        a = self.r1 - offset
        b = self.r2 - offset
        if not self.chunk:
            return [[a, b]]
        return [[x, y]
                for x in range(a - self.window, a + self.window + 1)
                for y in range(b - self.window, b + self.window + 1)]

    # True for every frame in which this interaction is formed.
    def mask(self, traj):
        # The scheme has to match the model. 1115 hardcoded 'ca', which is right for the
        # all-atom topologies it was written for but finds nothing in a coarse-grained one,
        # where the beads are named after the residue; mdtraj then warns and drops every
        # pair, leaving an empty array.
        [cont, pairs] = md.compute_contacts(traj, contacts=self.pairs(),
                                            scheme=self.profile.contact_scheme,
                                            ignore_nonprotein=True)
        return np.any(cont < self.cutoff, axis=1)

    def describe(self):
        return {'r1': self.r1, 'r2': self.r2, 'window': self.window if self.chunk else 0,
                'cutoff': self.cutoff, 'chunk': bool(self.chunk),
                'profile': self.profile.name, 'n_pairs': len(self.pairs())}

    def label(self):
        return '{0}_{1}'.format(self.r1, self.r2)

    def __repr__(self):
        return 'ChunkContact({0}, {1}, window={2}, cutoff={3})'.format(
            self.r1, self.r2, self.window, self.cutoff)


# Several interactions treated as one condition.
#
#   mode 'any'  a frame counts when at least one of them is formed
#   mode 'all'  a frame counts only when every one of them is formed
#
# 'any' asks what the ensemble looks like when this group of contacts is in play at all;
# 'all' asks what it looks like when they are formed together, which is the stricter
# question and will usually leave far fewer frames on the 'with' side.
class InteractionSet(object):

    MODES = ('any', 'all')

    def __init__(self, selectors, mode='any'):
        self.selectors = list(selectors)
        if not self.selectors:
            raise ValueError('InteractionSet needs at least one interaction')
        if mode not in self.MODES:
            raise ValueError('mode must be one of {0}, got {1!r}'.format(self.MODES, mode))
        self.mode = mode

    def mask(self, traj):
        combined = None
        for selector in self.selectors:
            current = selector.mask(traj)
            if combined is None:
                combined = current
            elif self.mode == 'any':
                combined = combined | current
            else:
                combined = combined & current
        return combined

    def describe(self):
        return {'mode': self.mode,
                'n_interactions': len(self.selectors),
                'interactions': [s.describe() for s in self.selectors]}

    def label(self):
        return (self.mode + '_') + '__'.join(s.label() for s in self.selectors)

    def __repr__(self):
        return 'InteractionSet({0}, mode={1!r})'.format(
            ', '.join(repr(s) for s in self.selectors), self.mode)


# Build a selector for one or more residue pairs.
def interaction_set(pairs, mode='any', window=1, cutoff=None, chunk=True, profile=None):
    """
    :param pairs: iterable of (r1, r2), numbered as in the contact table
    :param mode: 'any' for frames where at least one is formed, 'all' for frames where
                 every one is formed
    :return: InteractionSet
    """
    selectors = [ChunkContact(r1, r2, window=window, cutoff=cutoff, chunk=chunk, profile=profile)
                 for r1, r2 in pairs]
    return InteractionSet(selectors, mode=mode)


# Get frame numbers with the target interaction
def inter_frames(r1, r2, traj, chunk=True, cutoff=None, window=1, profile=None):
    """
    :return: boolean mask over frames, True where the interaction is present
    """
    selector = ChunkContact(r1, r2, window=window, cutoff=cutoff, chunk=chunk, profile=profile)
    return frames_where(traj, selector)


# Apply a selector to a trajectory. Works for a single interaction or a set of them.
def frames_where(traj, selector):
    return selector.mask(traj)


# Use r[contact_select] to call the filtered trajectory

# Search other interactions with the strongest interaction.
def new_interaction_map(name, traj, contact_select, cutoff=None, profile=None, output_dir=None):
    profile = profiles.get_profile(profile)
    cutoff = profile.cutoff if cutoff is None else cutoff
    contact = cg.ContactProbData(name + 'with', cutoff, traj[contact_select], profile=profile,
                                 output_dir=output_dir)
    contact_non = cg.ContactProbData(name + 'without', cutoff, traj[~contact_select], profile=profile,
                                     output_dir=output_dir)
    return contact, contact_non


# Compare radius of gyration and helicity

# Calculate the Hydrogen-bond number using Mdtraj library
# Here the r is the trajectory that only containing the protein atoms
# j is the number of frames in the trajectory.
def calc_HB(r, j, length):
    # Calculate the Hydogen_bond
    hbo = md.wernet_nilsson(r)
    # This is a couter for the loop over frames
    op = 0
    # A list to store the final result.
    HB = []
    # looping over each frame
    while op < j:
        # hbo[op] is a sub numpy array stored donor and acceptor
        # information of the H-bonds. The length of the array is the number of
        # H-bonds
        if (len(hbo[op])):
            HB.append(len(hbo[op]) / length)
        op += 1
    # Average H-bonds per residue per frame
    return (sum(HB) / j)


# Calculate the helical propensity fo the protein.
# Here the t is the full trajectory file(We can try to use r instead)
def calc_Heli(t, length):
    # Calculate the secondary structure using dssp
    dssp = md.compute_dssp(t)
    # Create an empty array with the row number equals to the number of residues
    # We may not need the number 1 here
    dssp_count = np.zeros((1, t.n_residues))
    # Loop over frames
    for i in range(t.n_frames):
        # Loop over residues
        for j in range(t.n_residues):
            # In dssp, letter 'H' means alpha-helix
            if dssp[i, j] in 'H':
                dssp_count[0, j] += 1
    # Get the average of the dssp_count
    dssp_prob = np.divide(dssp_count, t.n_frames)
    summary = dssp_prob.sum(axis=1)
    return (summary[0] / length)


# calc_beta is added 03/28/2022
# Same algorithm as alpha
def calc_Beta(t, length):
    dssp = md.compute_dssp(t)
    dssp_count = np.zeros((1, t.n_residues))
    for i in range(t.n_frames):
        for j in range(t.n_residues):
            if dssp[i, j] in 'E':
                dssp_count[0, j] += 1
    dssp_prob = np.divide(dssp_count, t.n_frames)
    summary = dssp_prob.sum(axis=1)
    return (summary[0] / length)


# Average radius of gyration of the ensemble
def calc_Rg(r):
    d = md.compute_rg(r)
    return (st.mean(d))


# Average end to end distance of the ensemble
# The defination of the end to end distance is the distance between
# the first alpha carbon and the last alpha carbon of the protein sequence.
#
# Note that the alpha carbon indices come from t while the distances are measured in r.
# That works because the protein atoms sit at the front of the topology and keep their
# indices when it is sliced. Left as it is: it is how every existing number was produced.
def calc_Re(r, t, profile=None):
    profile = profiles.get_profile(profile)
    topology = t.topology
    if profile.has_backbone:
        rpology = topology.select_atom_indices(selection='alpha')
    else:
        # One bead per residue, so there are no alpha carbons to select. The chain ends are
        # the first and last bead.
        rpology = [0, topology.n_atoms - 1]
    d = md.compute_distances(r, [[rpology[0], rpology[-1]]])
    listtemp = []
    for temp in d:
        listtemp.append(float(temp[0]))
    return (st.mean(listtemp))


# Which trajectory object each observable sees is deliberate and inherited:
# Rg and the end to end distance are measured on the protein slice, while the secondary
# structure and the hydrogen bonds are measured on the full trajectory. Changing that
# would change every number this has ever produced.
def calculate_ensemble(traj, sliced_traj, length, profile=None):
    profile = profiles.get_profile(profile)
    values = {'Rg': calc_Rg(sliced_traj),
              'Re': calc_Re(sliced_traj, traj, profile=profile),
              'Frames': traj.n_frames}
    if profile.has_backbone:
        values['Heli'] = calc_Heli(traj, length)
        values['Beta'] = calc_Beta(traj, length)
        values['HB'] = calc_HB(traj, traj.n_frames, length)
    else:
        # DSSP needs a backbone and the hydrogen bond search needs hydrogens. Neither
        # exists in a one bead per residue model, so these are not measurable rather than
        # zero, and saying so is better than reporting a number that means nothing.
        values['Heli'] = float('nan')
        values['Beta'] = float('nan')
        values['HB'] = float('nan')
    return values


def compare_ensemble(traj, sliced_traj, length, contact_select, profile=None):
    """
    Ensemble properties of the frames that contain the interaction and those that do not.

    'Frames' records how many frames went into each side. A split where one side holds a
    handful of frames means something very different from an even one, and the original
    output did not say which it was.
    """
    return {'interaction': calculate_ensemble(traj[contact_select], sliced_traj[contact_select],
                                              length, profile=profile),
            'non_interaction': calculate_ensemble(traj[~contact_select], sliced_traj[~contact_select],
                                                  length, profile=profile)}
