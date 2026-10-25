"""Normalization: the reference curve, the bins, and the edge cases around them."""
import numpy as np
import pandas as pd
import pytest

import normalization as nl
import profiles

from conftest import probe_table


def test_fitting_function_is_the_power_law():
    row = pd.Series({'distance': 10})
    assert nl.fitting_function(row, 1.64, -1.33) == 1.64 * 10 ** -1.33


def test_reference_probability_matches_fitting_function():
    """The helper and the row-wise form must agree, on scalars and on arrays."""
    for sep in (3, 7, 25, 60):
        row = pd.Series({'distance': sep})
        assert nl.reference_probability(sep, 1.64, -1.33) == nl.fitting_function(row, 1.64, -1.33)
    seps = np.array([3, 7, 25, 60])
    expected = np.array([nl.fitting_function(pd.Series({'distance': s}), 1.64, -1.33) for s in seps])
    assert np.array_equal(nl.reference_probability(seps, 1.64, -1.33), expected)


def test_adds_the_three_derived_columns(pair_table):
    out = nl.normalization(pair_table, profile='all_atom')
    for column in ('gs_standard', 'relative_strength', 'plot_value'):
        assert column in out.columns


def test_uses_the_profile_constants(pair_table):
    aa = nl.normalization(pair_table.copy(), profile='all_atom')
    cg = nl.normalization(pair_table.copy(), profile='cg')
    sep = aa['distance'].iloc[0]
    assert aa['gs_standard'].iloc[0] == 1.64 * sep ** -1.33
    assert cg['gs_standard'].iloc[0] == 13.12 * sep ** -2.32


def test_explicit_arguments_override_the_profile(pair_table):
    out = nl.normalization(pair_table, profile='all_atom', a1=2.0, b1=-1.0)
    sep = out['distance'].iloc[0]
    assert out['gs_standard'].iloc[0] == 2.0 * sep ** -1.0


@pytest.mark.parametrize('target,expected', [
    (2.0, 2), (1.6, 2),          # above the strong attractive edge
    (1.5, 2),                    # exactly on it
    (1.0, 1), (0.5, 1),          # weak attractive, inclusive at 0.5
    (0.4, 0), (0.0, 0), (-0.9, 0),
    (-1.0, -1), (-1.9, -1),      # weak repulsive, inclusive at -1
    (-2.0, -2), (-3.0, -2),      # strong repulsive, inclusive at -2
])
def test_all_atom_bin_edges(target, expected):
    out = nl.normalization(probe_table([target]), profile='all_atom')
    assert int(out['plot_value'].iloc[0]) == expected


@pytest.mark.parametrize('target,expected', [
    (2.5, 2), (2.0, 2),
    (1.5, 1), (1.0, 1),
    (0.5, 0), (-0.5, 0),
    (-1.0, -1), (-1.5, -1),
    (-2.0, -2), (-3.0, -2),
])
def test_cg_bin_edges(target, expected):
    out = nl.normalization(probe_table([target]), profile='cg', a1=1.64, b1=-1.33)
    assert int(out['plot_value'].iloc[0]) == expected


def test_relative_strength_is_the_log_ratio():
    out = nl.normalization(probe_table([0.75]), profile='all_atom')
    expected = np.log(out['cont_prob'].iloc[0] / out['gs_standard'].iloc[0])
    assert out['relative_strength'].iloc[0] == pytest.approx(expected, abs=1e-15)


def test_zero_contact_probability_is_neutral_not_infinite():
    """p_obs == 0 must not produce -inf or a repulsive call."""
    table = pd.DataFrame({'r_1': [1], 'r_2': [11], 'distance': [10], 'cont_prob': [0.0]})
    out = nl.normalization(table, profile='all_atom')
    assert out['relative_strength'].iloc[0] == 0.0
    assert int(out['plot_value'].iloc[0]) == 0
    assert np.isfinite(out['relative_strength']).all()


def test_prob_floor_defaults_to_off():
    """The shipped behaviour applies no contact-probability cutoff."""
    assert profiles.ALL_ATOM.prob_floor == 0.0
    assert profiles.CG.prob_floor == 0.0
    out = nl.normalization(probe_table([2.0]), profile='all_atom')
    assert int(out['plot_value'].iloc[0]) == 2


def test_prob_floor_suppresses_low_probability_pairs():
    out = nl.normalization(probe_table([2.0, -3.0]), profile='all_atom', prob_floor=1.0)
    assert set(out['plot_value']) == {0}


def test_prob_floor_keeps_pairs_above_it():
    table = probe_table([2.0])
    floor = table['cont_prob'].iloc[0] / 2
    out = nl.normalization(table, profile='all_atom', prob_floor=floor)
    assert int(out['plot_value'].iloc[0]) == 2


def test_cg_reference_exceeds_one_at_short_separation():
    """Documented behaviour, kept deliberately.

    Below separation 4 the CG curve returns a 'probability' above 1, so bonded neighbours,
    which are always in contact, are classified strongly repulsive. min_sep removes them
    from the figure but not from the table.
    """
    cg = profiles.CG
    for sep in (1, 2, 3):
        assert nl.reference_probability(sep, cg.reference_a, cg.reference_b) > 1.0
    assert nl.reference_probability(4, cg.reference_a, cg.reference_b) < 1.0

    bonded = pd.DataFrame({'r_1': [1], 'r_2': [2], 'distance': [1], 'cont_prob': [1.0]})
    out = nl.normalization(bonded, profile='cg')
    assert int(out['plot_value'].iloc[0]) == -2
    assert out['distance'].iloc[0] < cg.min_sep


def test_normalization_writes_into_the_frame_it_is_given(pair_table):
    """Documents the in-place behaviour, which callers rely on."""
    returned = nl.normalization(pair_table, profile='all_atom')
    assert returned is pair_table
