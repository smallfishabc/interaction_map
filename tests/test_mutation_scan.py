"""Mutation scanning: the Eq. 3 sliding window and the candidate selection."""
import numpy as np
import pandas as pd
import pytest

from interaction_map.util import mutation_scan as ms

from conftest import requires_p53

SEQ = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPD'


def lookup(pairs):
    return ms.strength_lookup(
        pd.DataFrame([{'r_1': i, 'r_2': j, 'relative_strength': v} for i, j, v in pairs]))


# --------------------------------------------------------------- Eq. 3

def test_centre_of_the_window_counts_in_full():
    assert ms.chunk_strength(lookup([(20, 40, 1.0)]), 20, 40) == 1.0


@pytest.mark.parametrize('pair,weight', [
    ((21, 40), 0.5), ((19, 40), 0.5), ((20, 41), 0.5), ((20, 39), 0.5),   # one step
    ((21, 41), 0.25), ((19, 39), 0.25), ((21, 39), 0.25), ((19, 41), 0.25),  # two steps
])
def test_neighbours_are_halved_per_step(pair, weight):
    """Eq. 3: the weight is 1 / 2 ** (|i - r1| + |j - r2|)."""
    assert ms.chunk_strength(lookup([(pair[0], pair[1], 1.0)]), 20, 40) == weight


def test_full_window_sums_the_weights():
    pairs = [(i, j, 1.0) for i in range(19, 22) for j in range(39, 42)]
    expected = sum(1.0 / 2 ** (abs(i - 20) + abs(j - 40)) for i, j, _ in pairs)
    assert ms.chunk_strength(lookup(pairs), 20, 40) == pytest.approx(expected)
    assert expected == pytest.approx(4.0)


def test_pairs_outside_the_window_are_ignored():
    assert ms.chunk_strength(lookup([(23, 40, 1.0)]), 20, 40) == 0.0


def test_missing_pairs_contribute_nothing():
    """Pairs dropped by a cutoff are absent, not zero-strength."""
    assert ms.chunk_strength(lookup([(20, 40, 2.0)]), 20, 40) == 2.0


def test_a_wider_window_reaches_further():
    assert ms.chunk_strength(lookup([(22, 40, 1.0)]), 20, 40, window=2) == 0.25
    assert ms.chunk_strength(lookup([(22, 40, 1.0)]), 20, 40, window=1) == 0.0


def test_strength_is_signed():
    """A repulsive neighbourhood must not look attractive."""
    assert ms.chunk_strength(lookup([(20, 40, -2.0)]), 20, 40) == -2.0


# --------------------------------------------------------------- residue tables

def test_every_standard_residue_is_classified():
    for residue in 'ACDEFGHIKLMNPQRSTVWY':
        assert residue in ms.RESIDUE_TYPE
        assert residue in ms.KYTE_DOOLITTLE
        assert residue in ms.CHARGE


def test_charges_follow_the_residue_classes():
    assert ms.CHARGE['K'] == 1 and ms.CHARGE['R'] == 1 and ms.CHARGE['H'] == 1
    assert ms.CHARGE['D'] == -1 and ms.CHARGE['E'] == -1
    assert ms.CHARGE['A'] == 0


def test_hydropathy_has_the_expected_extremes():
    """Kyte-Doolittle: isoleucine most hydrophobic, arginine least."""
    assert max(ms.KYTE_DOOLITTLE, key=ms.KYTE_DOOLITTLE.get) == 'I'
    assert min(ms.KYTE_DOOLITTLE, key=ms.KYTE_DOOLITTLE.get) == 'R'


def test_residue_table_covers_everything():
    table = ms.residue_table()
    assert len(table) == 20
    assert set(table.columns) == {'residue', 'residue_type', 'kyte_doolittle', 'charge'}


# --------------------------------------------------------------- selection

def interaction_table():
    rows = []
    for r1 in range(1, 40):
        for r2 in range(r1 + 3, 61):
            rows.append({'r_1': r1, 'r_2': r2, 'cont_prob': 0.05,
                         'distance': r2 - r1, 'relative_strength': 0.8, 'plot_value': 1})
    return pd.DataFrame(rows)


def test_weak_contacts_are_dropped():
    table = interaction_table()
    table.loc[table.index[:10], 'cont_prob'] = 0.001
    kept = ms.select_pairs(table, contact_probability_cutoff=0.01, min_separation=5)
    assert (kept['cont_prob'] >= 0.01).all()


def test_short_range_pairs_are_dropped():
    kept = ms.select_pairs(interaction_table(), min_separation=10)
    assert kept['distance'].min() >= 10


def test_terminal_pairs_are_dropped():
    """Interactions at the chain ends are one-sided and read stronger than they are."""
    kept = ms.select_pairs(interaction_table(), min_separation=5, sequence_length=60)
    assert (kept['r_1'] != 1).all()
    assert (kept['r_2'] != 60).all()


def test_termini_can_be_kept():
    kept = ms.select_pairs(interaction_table(), min_separation=5, sequence_length=60,
                           drop_termini=False)
    assert (kept['r_1'] == 1).any()


# --------------------------------------------------------------- chunk properties

def test_chunk_is_the_three_residues_around_the_position():
    scored = ms.add_chunk_strength(ms.select_pairs(interaction_table(), min_separation=5,
                                                   sequence_length=len(SEQ)), SEQ)
    row = scored[scored['r_1'] == 24].iloc[0]
    assert row['r_1_chunk'] == SEQ[22:25]          # positions 23, 24, 25


def test_chunk_properties_are_summed_over_the_chunk():
    scored = ms.add_chunk_strength(ms.select_pairs(interaction_table(), min_separation=5,
                                                   sequence_length=len(SEQ)), SEQ)
    row = scored[scored['r_1'] == 24].iloc[0]
    chunk = row['r_1_chunk']
    assert row['r_1_hydro'] == pytest.approx(sum(ms.KYTE_DOOLITTLE[a] for a in chunk))
    assert row['r_1_charge'] == sum(ms.CHARGE[a] for a in chunk)
    assert row['r_1_aromatic'] == sum(ms.RESIDUE_TYPE[a] == 'aromatic' for a in chunk)


def test_chunks_are_clipped_at_the_sequence_ends():
    scored = ms.add_chunk_strength(
        ms.select_pairs(interaction_table(), min_separation=5, drop_termini=False), SEQ)
    row = scored[scored['r_1'] == 1].iloc[0]
    assert row['r_1_chunk'] == SEQ[0:2]


def test_an_empty_selection_still_returns_the_columns():
    empty = ms.add_chunk_strength(interaction_table().iloc[0:0], SEQ)
    for column in ('chunk_strength', 'r_1_chunk', 'r_1_hydro', 'r_2_aromatic'):
        assert column in empty.columns


# --------------------------------------------------------------- filtering

def scored_table():
    return ms.add_chunk_strength(
        ms.select_pairs(interaction_table(), min_separation=5, sequence_length=len(SEQ)), SEQ)


def test_attractive_and_repulsive_are_separated():
    scored = scored_table()
    scored.loc[scored.index[:50], 'plot_value'] = -1
    attractive = ms.filter_candidates(scored, 'attractive', 'aromatic')
    repulsive = ms.filter_candidates(scored, 'repulsive', 'charge')
    assert (attractive['plot_value'] > 0).all()
    assert (repulsive['plot_value'] < 0).all()


def test_aromatic_filter_needs_an_aromatic_on_both_sides():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    assert ((candidates['r_1_aromatic'] > 0) & (candidates['r_2_aromatic'] > 0)).all()


def test_attractive_charge_filter_wants_opposite_charges():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'charge')
    assert ((candidates['r_1_charge'] * candidates['r_2_charge']) < 0).all()


def test_repulsive_charge_filter_wants_like_charges():
    scored = scored_table()
    scored['plot_value'] = -1
    candidates = ms.filter_candidates(scored, 'repulsive', 'charge')
    assert ((candidates['r_1_charge'] * candidates['r_2_charge']) > 0).all()


def test_unknown_filters_are_rejected():
    with pytest.raises(ValueError):
        ms.filter_candidates(scored_table(), 'sideways', 'aromatic')
    with pytest.raises(ValueError):
        ms.filter_candidates(scored_table(), 'attractive', 'sticky')


# --------------------------------------------------------------- mutants

def test_mutant_carries_one_substitution_at_the_named_position():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    mutants = ms.design_mutants(candidates, SEQ, name='p53', limit=3)
    for _, row in mutants.iterrows():
        position = int(row['positions'])
        assert len(row['sequence']) == len(SEQ)
        assert row['sequence'][position - 1] == 'E'
        differing = [k for k in range(len(SEQ)) if SEQ[k] != row['sequence'][k]]
        assert differing == [position - 1]


def test_mutant_name_records_the_substitution():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    row = ms.design_mutants(candidates, SEQ, name='p53', limit=1).iloc[0]
    position = int(row['positions'])
    assert row['name'] == 'p53_{0}{1}E'.format(SEQ[position - 1], position)


def test_double_mutants_change_both_sides():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    row = ms.design_mutants(candidates, SEQ, name='p53', both_sides=True, limit=1).iloc[0]
    differing = [k for k in range(len(SEQ)) if SEQ[k] != row['sequence'][k]]
    assert len(differing) == 2


def test_no_duplicate_mutants():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    mutants = ms.design_mutants(candidates, SEQ, name='p53')
    assert len(set(mutants['sequence'])) == len(mutants)


def test_limit_is_respected():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    assert len(ms.design_mutants(candidates, SEQ, name='p53', limit=4)) == 4


def test_substitution_can_be_chosen():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    row = ms.design_mutants(candidates, SEQ, name='p53', substitution='K', limit=1).iloc[0]
    assert row['sequence'][int(row['positions']) - 1] == 'K'


def test_an_unknown_combination_needs_an_explicit_substitution():
    candidates = ms.filter_candidates(scored_table(), 'attractive', 'aromatic')
    with pytest.raises(ValueError):
        ms.design_mutants(candidates, SEQ, interaction_type='attractive',
                          residue_types='nonsense')


# --------------------------------------------------------------- against the real map

@requires_p53
def test_scan_finds_the_known_p53_interaction_region(p53_reference, p53_sequence):
    """The strong attractive region of p53 NTAD involves W23/K24 against W53/F54.

    Those positions appear among the mutants reported in the paper, so a scan that no
    longer reaches them has changed.
    """
    mutants, candidates = ms.scan(p53_reference, p53_sequence, name='p53',
                                  interaction_type='attractive', residue_types='aromatic')
    assert not candidates.empty
    top = candidates.iloc[0]
    assert {int(top['r_1']), int(top['r_2'])} <= {23, 24, 53, 54}
    positions = {int(p) for p in mutants['positions']}
    assert {23, 24, 53, 54} <= positions


@requires_p53
def test_scan_ranks_by_chunk_strength(p53_reference, p53_sequence):
    _mutants, candidates = ms.scan(p53_reference, p53_sequence, name='p53')
    strengths = candidates['chunk_strength'].values
    assert (np.diff(strengths) <= 1e-12).all()
