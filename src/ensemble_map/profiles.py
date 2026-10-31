# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 2026

Simulation profiles for the Interaction Map.

Every constant that differs between the all-atom analysis (CAMPARI/ABSINTH) and the
coarse-grained analysis (CALVADOS) is collected here. The two sets of numbers used to
live on two separate branches (0811_final_all_atom_version and 0519_CG), which meant
editing one silently left the other behind. Keeping them side by side makes the
difference explicit and lets a single code path serve both.

@author: ShaharGroup-fyu
"""

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class Profile:
    """
    One complete set of analysis constants.

    Attributes:
        name                Profile identifier, also used in output file names
        cutoff              Contact cutoff distance, nm
        pair_selection      'mdtraj_ca'  -> md.compute_contacts(contacts='all', scheme='CA'),
                                            which already drops pairs with |i-j| < 3
                            'all_pairs'  -> top.select_pairs('all', 'all'), which keeps every
                                            pair. Valid only for one-bead-per-residue models.
        has_backbone        Whether the model carries a real protein backbone. All-atom
                            does, so alpha carbons, DSSP and hydrogen bonds are available.
                            One bead per residue does not, so helicity, beta content and
                            hydrogen bonds are not defined and the end to end distance has
                            to be measured between the first and last bead instead of
                            between alpha carbons.
        contact_scheme      How the distance between two residues is measured.
                            all-atom: 'ca', the alpha carbons.
                            CG: 'closest-heavy', mdtraj's default. CALVADOS beads are named
                            after the residue rather than CA, so asking for 'ca' finds no
                            alpha carbon and silently discards every pair.
        slice_protein       Slice the topology down to 'protein' before computing contacts
        residue_index_offset
                            Added to the raw mdtraj residue index so that r_1 / r_2 are
                            1-based sequence positions, which is what both the plot layer and
                            every existing CSV expect.
                            all-atom: 0. The CAMPARI topologies carry an ACE cap at index 0 and
                            an NME cap at the end. ignore_nonprotein drops the caps from the
                            pair list but not from the numbering, so the first real residue is
                            already index 1 and no shift is needed.
                            CG: 1. CALVADOS topologies have no caps, so the first residue is
                            index 0 and has to be shifted up by one.
                            This is why an UNCAPPED all-atom topology would come out one
                            residue low. check_residue_numbering() in contact_map_generation
                            catches that case.
        reference_a         Reference contact probability curve, p_ref(s) = a * s ** b
        reference_b         Exponent of the reference curve
        inter_cutoff        Interaction strength bin edges (strong att, weak att, weak rep,
                            strong rep), applied to ln(p_obs / p_ref)
        value_list          Interaction levels assigned to those bins
        prob_floor          Contact probabilities at or below this value are forced to level 0.
                            The published method applies a cutoff to p_1 as well as to I;
                            0.0 reproduces the current branches, which dropped it.
        min_sep             Minimum sequence separation drawn on the map
        line_width_scale    Multiplier on the arc line width
        residue_style       'marker' -> coloured dots (CG branch, newer)
                            'letter' -> coloured one-letter codes (used by the paper figures)
    """

    name: str
    cutoff: float
    pair_selection: str
    contact_scheme: str
    has_backbone: bool
    slice_protein: bool
    residue_index_offset: int
    reference_a: float
    reference_b: float
    inter_cutoff: Tuple[float, float, float, float]
    value_list: Tuple[int, int, int, int, int]
    prob_floor: float
    min_sep: int
    line_width_scale: float
    residue_style: str


# All-atom profile. Constants as they stand on 0811_final_all_atom_version.
ALL_ATOM = Profile(
    name='all_atom',
    cutoff=0.8,
    pair_selection='mdtraj_ca',
    contact_scheme='ca',
    has_backbone=True,
    slice_protein=True,
    residue_index_offset=0,
    reference_a=1.64,
    reference_b=-1.33,
    inter_cutoff=(1.5, 0.5, -1, -2),
    value_list=(2, 1, 0, -1, -2),
    prob_floor=0.0,
    min_sep=3,
    line_width_scale=2.0,
    residue_style='marker',
)

# Coarse-grained profile. Constants as they stand on 0519_CG (CALVADOS trajectories).
# The steeper exponent matches Fig. 8A of the manuscript: CG contact probability falls off
# roughly an order of magnitude faster than all-atom over the same separation range.
CG = Profile(
    name='cg',
    cutoff=1.2,
    pair_selection='all_pairs',
    contact_scheme='closest-heavy',
    has_backbone=False,
    slice_protein=False,
    residue_index_offset=1,
    reference_a=13.12,
    reference_b=-2.32,
    inter_cutoff=(2, 1, -1, -2),
    value_list=(2, 1, 0, -1, -2),
    prob_floor=0.0,
    min_sep=5,
    line_width_scale=4.0,
    residue_style='marker',
)

PROFILES = {'all_atom': ALL_ATOM, 'cg': CG}

DEFAULT_PROFILE = ALL_ATOM


# Look up a profile by name so the command line and the library agree on the spelling.
def get_profile(profile):
    """
    :param profile: profile name, or an already built Profile
    :type profile: str/Profile/None
    :return: Profile
    """
    if profile is None:
        return DEFAULT_PROFILE
    if isinstance(profile, Profile):
        return profile
    try:
        return PROFILES[profile]
    except KeyError:
        raise ValueError(
            'Unknown profile {!r}. Available profiles: {}'.format(
                profile, ', '.join(sorted(PROFILES))))
