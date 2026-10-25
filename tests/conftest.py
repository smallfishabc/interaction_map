"""Shared fixtures.

The tests fall into two groups. Most run on synthetic data and are always available. A few
compare against the real August 2023 p53 output in All_atom_p53_interaction, or against a
CALVADOS mutant in CG_p53_random; those skip when the data folders are absent, so the
suite still runs on a clean checkout.
"""
import os
import sys

import pytest

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT not in sys.path:
    sys.path.insert(0, PROJECT)

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd

P53_DIR = os.path.join(PROJECT, 'All_atom_p53_interaction', 'p53-summary')
P53_RUN = os.path.join(P53_DIR, 'BB', 'S_0')
CG_DIR = os.path.join(PROJECT, 'CG_p53_random')

requires_p53 = pytest.mark.skipif(
    not os.path.exists(os.path.join(P53_RUN, 'interaction_0317.csv')),
    reason='All_atom_p53_interaction reference data not present')

requires_cg = pytest.mark.skipif(
    not os.path.exists(os.path.join(CG_DIR, 'A39G_K24Y.pdb')),
    reason='CG_p53_random data not present')


@pytest.fixture(scope='session')
def p53_sequence():
    with open(os.path.join(P53_DIR, 'seq.fasta')) as handle:
        return handle.read().strip()


@pytest.fixture(scope='session')
def p53_reference():
    """The August 2023 interaction table, 1711 pairs."""
    return pd.read_csv(os.path.join(P53_RUN, 'interaction_0317.csv'), index_col=0)


@pytest.fixture(scope='session')
def p53_png():
    return os.path.join(P53_RUN, 'p53.png')


@pytest.fixture
def pair_table():
    """A small contact table with 1-based positions and separations from 3 up."""
    rows = []
    for i in range(1, 21):
        for j in range(i + 3, 21):
            rows.append((i, j, 0.05, j - i))
    return pd.DataFrame(rows, columns=['r_1', 'r_2', 'cont_prob', 'distance'])


def probe_table(targets, separation=10, a=1.64, b=-1.33):
    """Pairs whose ln(p_obs / p_ref) equals each requested value exactly."""
    ref = a * separation ** b
    return pd.DataFrame({
        'r_1': list(range(1, len(targets) + 1)),
        'r_2': [i + separation for i in range(1, len(targets) + 1)],
        'distance': [separation] * len(targets),
        'cont_prob': [ref * np.exp(t) for t in targets],
    })


def synthetic_trajectory(n_residues=12, n_frames=5, atoms_per_residue=1, seed=0):
    """A minimal in-memory mdtraj trajectory.

    atoms_per_residue=1 mimics a coarse-grained bead model, where atom and residue indices
    coincide. The single atom is named CA so the 'CA' contact scheme can find it.
    """
    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    for _ in range(n_residues):
        residue = top.add_residue('ALA', chain)
        for name in ('CA', 'CB')[:atoms_per_residue]:
            top.add_atom(name, md.element.carbon, residue)

    rng = np.random.default_rng(seed)
    n_atoms = n_residues * atoms_per_residue
    # A loose self-avoiding-ish walk so some pairs fall inside any sensible cutoff.
    xyz = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    for frame in range(n_frames):
        position = np.zeros(3)
        for atom in range(n_atoms):
            position = position + rng.normal(0.0, 0.38, 3)
            xyz[frame, atom] = position
    return md.Trajectory(xyz, top)
