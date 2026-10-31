# -*- coding: utf-8 -*-
"""
Pick mutation candidates from an interaction map.

A strong interaction between two residues makes their neighbours look connected as well,
so reading single pairs off the map tends to point at a cluster rather than at a residue.
The sliding window answers that: each pair is scored together with its neighbours, three
residues either side, with the neighbours counting for less.

    I_chunk = sum over i,j in the 3x3 window of  I_ij / 2 ** (|i - r1| + |j - r2|)

which is Eq. 3 of the paper. A pair at the centre counts fully, a pair one step away counts
a half, two steps a quarter.

Ranking by that, filtered by what the residues are, gives the list of positions worth
mutating. Which mutation to make depends on what is being changed: an attractive contact
between aromatics is broken by putting a charge there, a repulsive contact between like
charges is relieved by removing one.

Independent of the analysis modules, as `calibrate` is: it reads a table that has already
been normalized and does not import the pipeline that produced it.

@author: ShaharGroup-fyu
"""
from __future__ import print_function

import numpy as np
import pandas as pd

# Residue classes. The same table the published analysis used, kept here so nothing has to
# be downloaded at run time.
RESIDUE_TYPE = {
    'Y': 'aromatic', 'F': 'aromatic', 'W': 'aromatic',
    'R': 'positive', 'H': 'positive', 'K': 'positive',
    'D': 'negative', 'E': 'negative',
    'S': 'polar', 'T': 'polar', 'Q': 'polar', 'N': 'polar',
    'A': 'hydrophobic', 'V': 'hydrophobic', 'I': 'hydrophobic', 'L': 'hydrophobic',
    'M': 'hydrophobic', 'C': 'hydrophobic', 'G': 'hydrophobic', 'P': 'hydrophobic',
}

# Kyte-Doolittle hydropathy. Positive is hydrophobic.
KYTE_DOOLITTLE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
    'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3,
}

CHARGE = {aa: (1 if t == 'positive' else -1 if t == 'negative' else 0)
          for aa, t in RESIDUE_TYPE.items()}

# Glycine's hydropathy, the default boundary between a hydrophobic and a hydrophilic chunk.
HYDROPHOBICITY_CUTOFF = -1.2


def residue_table():
    """The residue classes and scales as a dataframe, for inspection."""
    rows = [{'residue': aa, 'residue_type': RESIDUE_TYPE[aa],
             'kyte_doolittle': KYTE_DOOLITTLE[aa], 'charge': CHARGE[aa]}
            for aa in sorted(RESIDUE_TYPE)]
    return pd.DataFrame(rows)


def select_pairs(interaction, contact_probability_cutoff=0.01, min_separation=5,
                 sequence_length=None, drop_termini=True):
    """
    Narrow an interaction table to the pairs worth ranking.

    :param interaction: normalized table, with r_1, r_2, cont_prob, distance,
                        relative_strength and plot_value
    :param contact_probability_cutoff: pairs seen less often than this are dropped
    :param min_separation: pairs closer than this along the sequence are dropped
    :param sequence_length: needed to drop the last residue when drop_termini is set
    :param drop_termini: drop pairs involving the first or last residue, whose
                         interactions are one-sided and read as stronger than they are
    :return: the surviving rows, most strongly interacting first
    """
    keep = ((interaction['cont_prob'] >= contact_probability_cutoff)
            & (interaction['distance'] >= min_separation))
    if drop_termini:
        keep &= interaction['r_1'] != 1
        if sequence_length is not None:
            keep &= interaction['r_2'] != sequence_length
    return interaction[keep].sort_values('relative_strength', ascending=False).copy()


def strength_lookup(selected):
    """
    Index an interaction table by residue pair, ready for chunk_strength.

    :param selected: table with r_1, r_2 and relative_strength
    :return: Series keyed by (r_1, r_2)
    """
    lookup = selected.set_index(['r_1', 'r_2'])['relative_strength']
    return lookup[~lookup.index.duplicated()]


def chunk_strength(lookup, r1, r2, window=1):
    """
    Eq. 3: the interaction strength of a residue pair together with its neighbours.

    Pairs missing from the lookup, because they fell below a cutoff, contribute nothing.

    :param lookup: Series keyed by (r_1, r_2), from strength_lookup
    :param r1, r2: the pair at the centre of the window
    :param window: how far either side to reach. 1 gives the 3x3 window of the paper.
    :return: the weighted sum
    """
    total = 0.0
    for i in range(r1 - window, r1 + window + 1):
        for j in range(r2 - window, r2 + window + 1):
            try:
                value = float(lookup.loc[(i, j)])
            except (KeyError, TypeError):
                continue
            total += value / 2.0 ** (abs(i - r1) + abs(j - r2))
    return float(total)


def add_chunk_strength(selected, sequence, window=1):
    """
    Score every pair in a selection with Eq. 3, and describe the chunks around it.

    Adds, for each side of the pair: the three-residue chunk itself, its summed
    Kyte-Doolittle hydropathy, its net charge and its aromatic count. Those are what the
    choice of mutation is made on.

    :param selected: output of select_pairs
    :param sequence: the protein sequence, 1-based positions matching r_1 and r_2
    :param window: half width of the window
    :return: the same table with the extra columns, strongest chunk first
    """
    table = selected.copy()
    if table.empty:
        for column in ('chunk_strength', 'r_1_chunk', 'r_2_chunk', 'r_1_hydro',
                       'r_2_hydro', 'r_1_charge', 'r_2_charge', 'r_1_aromatic',
                       'r_2_aromatic'):
            table[column] = []
        return table

    lookup = strength_lookup(table)

    def chunk_of(position):
        start = max(position - window, 1)
        stop = min(position + window, len(sequence))
        return sequence[start - 1:stop]

    table['chunk_strength'] = [chunk_strength(lookup, int(r1), int(r2), window)
                               for r1, r2 in zip(table['r_1'], table['r_2'])]
    for side in ('r_1', 'r_2'):
        chunks = [chunk_of(int(p)) for p in table[side]]
        table[side + '_chunk'] = chunks
        table[side + '_hydro'] = [sum(KYTE_DOOLITTLE.get(a, 0.0) for a in c) for c in chunks]
        table[side + '_charge'] = [sum(CHARGE.get(a, 0) for a in c) for c in chunks]
        table[side + '_aromatic'] = [sum(RESIDUE_TYPE.get(a) == 'aromatic' for a in c)
                                     for c in chunks]
    return table.sort_values('chunk_strength', ascending=False)


# Which residue to put in, for each thing being changed. An attractive contact is broken
# by making the two sides repel; a repulsive one is relieved by removing the repulsion.
SUBSTITUTION = {
    ('attractive', 'aromatic'): 'E',
    ('attractive', 'hydrophobic'): 'E',
    ('attractive', 'charge'): 'E',
    ('repulsive', 'charge'): 'S',
    ('repulsive', 'hydrophilic'): 'W',
}


def filter_candidates(scored, interaction_type='attractive', residue_types='aromatic',
                      hydrophobicity_cutoff=HYDROPHOBICITY_CUTOFF):
    """
    Keep the pairs where both sides are the kind of residue being targeted.

    :param interaction_type: 'attractive' or 'repulsive'
    :param residue_types: 'charge', 'aromatic', 'hydrophobic' or 'hydrophilic'
    :return: the candidates, strongest first for attractive, strongest repulsion first
             for repulsive
    """
    if interaction_type not in ('attractive', 'repulsive'):
        raise ValueError("interaction_type must be 'attractive' or 'repulsive'")
    table = scored[scored['plot_value'] > 0] if interaction_type == 'attractive' \
        else scored[scored['plot_value'] < 0]

    if residue_types == 'charge':
        # Opposite charges attract, like charges repel.
        product = table['r_1_charge'] * table['r_2_charge']
        table = table[product < 0] if interaction_type == 'attractive' else table[product > 0]
    elif residue_types == 'aromatic':
        table = table[(table['r_1_aromatic'] * table['r_2_aromatic']) > 0]
    elif residue_types == 'hydrophobic':
        above = ((table['r_1_hydro'] - hydrophobicity_cutoff)
                 * (table['r_2_hydro'] - hydrophobicity_cutoff) > 0)
        table = table[above]
        table = table[(table['r_1_hydro'] - hydrophobicity_cutoff) > 0]
    elif residue_types == 'hydrophilic':
        same = ((table['r_1_hydro'] - hydrophobicity_cutoff)
                * (table['r_2_hydro'] - hydrophobicity_cutoff) > 0)
        table = table[same]
        table = table[(table['r_1_hydro'] - hydrophobicity_cutoff) < 0]
    else:
        raise ValueError("residue_types must be one of 'charge', 'aromatic', "
                         "'hydrophobic', 'hydrophilic'")

    ascending = interaction_type == 'repulsive'
    return table.sort_values('chunk_strength', ascending=ascending)


def design_mutants(candidates, sequence, name='mutant', interaction_type='attractive',
                   residue_types='aromatic', substitution=None, limit=None, both_sides=False):
    """
    Turn ranked candidates into mutant sequences.

    One position is changed per mutant by default, taken from each side of the ranked
    pairs in turn. `both_sides=True` changes both positions of a pair in one mutant.

    :param substitution: the residue to put in. Defaults to the usual choice for the
                         interaction and residue type being targeted.
    :param limit: keep at most this many mutants
    :return: dataframe of name and sequence, in rank order
    """
    replacement = substitution or SUBSTITUTION.get((interaction_type, residue_types))
    if replacement is None:
        raise ValueError('no default substitution for {0}/{1}; pass substitution='
                         .format(interaction_type, residue_types))

    seen = set()
    rows = []
    for _, row in candidates.iterrows():
        positions = [(int(row['r_1']), int(row['r_2']))] if both_sides \
            else [(int(row['r_1']),), (int(row['r_2']),)]
        for group in positions:
            if any(p < 1 or p > len(sequence) for p in group):
                continue
            key = tuple(sorted(group))
            if key in seen:
                continue
            seen.add(key)
            mutated = list(sequence)
            labels = []
            for p in group:
                labels.append('{0}{1}{2}'.format(sequence[p - 1], p, replacement))
                mutated[p - 1] = replacement
            rows.append({'name': '{0}_{1}'.format(name, '_'.join(labels)),
                         'sequence': ''.join(mutated),
                         'positions': '+'.join(str(p) for p in group),
                         'chunk_strength': float(row['chunk_strength'])})
            if limit and len(rows) >= limit:
                return pd.DataFrame(rows)
    return pd.DataFrame(rows)


def scan(interaction, sequence, name='mutant', interaction_type='attractive',
         residue_types='aromatic', contact_probability_cutoff=0.01, min_separation=5,
         window=1, limit=None, hydrophobicity_cutoff=HYDROPHOBICITY_CUTOFF):
    """
    The whole thing: narrow, score with Eq. 3, filter by residue type, design mutants.

    :return: (mutants, candidates) so the ranking behind the list can be inspected
    """
    selected = select_pairs(interaction, contact_probability_cutoff, min_separation,
                            sequence_length=len(sequence))
    scored = add_chunk_strength(selected, sequence, window=window)
    candidates = filter_candidates(scored, interaction_type, residue_types,
                                   hydrophobicity_cutoff)
    mutants = design_mutants(candidates, sequence, name, interaction_type, residue_types,
                             limit=limit)
    return mutants, candidates
