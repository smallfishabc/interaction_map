"""The reference curve fit.

The arithmetic is the one from the GS notebooks, so the first test here checks it against
that expression directly rather than against a restatement of it.
"""
import json
import os

import numpy as np
import pandas as pd
import pytest

from ensemble_map.util import calibrate as cal

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')
PDB = os.path.join(FIXTURES, 'A39G_K24Y.pdb')
DCD = os.path.join(FIXTURES, 'A39G_K24Y.dcd')


# --------------------------------------------------------------- the curve

def test_power_law_on_a_scalar():
    assert cal.power_law(10, 1.64, -1.33) == 1.64 * 10 ** -1.33


def test_power_law_on_an_array():
    got = cal.power_law(np.array([3, 10]), 1.64, -1.33)
    assert np.allclose(got, [1.64 * 3 ** -1.33, 1.64 * 10 ** -1.33])


# --------------------------------------------------------------- the averaging

def test_matches_the_notebook_expression():
    """Guards the one piece of arithmetic that was inherited rather than written."""
    rng = np.random.default_rng(0)
    separations = np.repeat(np.arange(3, 25), 7)
    distances = rng.uniform(0.2, 2.0, size=(len(separations), 500))
    cutoff = 0.8

    pair_df = pd.DataFrame(distances, index=separations)
    notebook = ((pair_df < cutoff).sum(axis=1) / len(pair_df.columns)) \
        .reset_index().rename(columns={'index': 'distance'}) \
        .groupby('distance').mean().squeeze()

    mine, _counts = cal.contact_probability_by_separation(separations, distances, cutoff)
    assert np.array_equal(notebook.values, mine.values)
    assert np.array_equal(notebook.index.values, mine.index.values)


def test_counts_the_pairs_behind_each_separation():
    separations = np.array([3, 3, 3, 4, 4])
    distances = np.zeros((5, 10))
    _p, counts = cal.contact_probability_by_separation(separations, distances, 1.0)
    assert counts.loc[3] == 3
    assert counts.loc[4] == 2


def test_probability_is_the_fraction_of_frames_inside_the_cutoff():
    # One pair, four frames, two of them in contact.
    distances = np.array([[0.1, 0.1, 9.0, 9.0]])
    p, _ = cal.contact_probability_by_separation(np.array([5]), distances, 1.0)
    assert p.loc[5] == 0.5


def test_a_pair_never_in_contact_contributes_zero():
    distances = np.array([[9.0, 9.0, 9.0]])
    p, _ = cal.contact_probability_by_separation(np.array([7]), distances, 1.0)
    assert p.loc[7] == 0.0


# --------------------------------------------------------------- the fit

@pytest.mark.parametrize('a,b', [(1.64, -1.33), (13.12, -2.32), (2.175, -1.555)])
def test_recovers_a_curve_it_was_given(a, b):
    separations = np.arange(3, 61)
    probability = pd.Series(cal.power_law(separations, a, b), index=separations)
    fit = cal.fit_power_law(probability)
    assert fit['a'] == pytest.approx(a, rel=1e-6)
    assert fit['b'] == pytest.approx(b, rel=1e-6)
    assert fit['r_squared'] == pytest.approx(1.0, abs=1e-9)


def test_residuals_are_observed_minus_predicted():
    separations = np.arange(3, 20)
    probability = pd.Series(cal.power_law(separations, 2.0, -1.5), index=separations)
    fit = cal.fit_power_law(probability)
    assert np.allclose(np.array(fit['observed']) - np.array(fit['predicted']),
                       fit['residuals'])


def test_separations_with_no_contacts_are_left_out_of_the_fit():
    separations = np.arange(3, 30)
    values = cal.power_law(separations, 2.0, -1.5)
    values[-5:] = 0.0                      # nothing observed at long separation
    fit = cal.fit_power_law(pd.Series(values, index=separations))
    assert fit['n_separations_dropped'] == 5
    assert fit['n_separations_fitted'] == len(separations) - 5
    assert fit['a'] == pytest.approx(2.0, rel=1e-6)


def test_too_few_usable_points_is_an_error_not_a_guess():
    probability = pd.Series([0.5, 0.0, 0.0], index=[3, 4, 5])
    with pytest.raises(ValueError):
        cal.fit_power_law(probability)


# --------------------------------------------------------------- the ceiling

def test_ceiling_finds_where_the_curve_stops_being_a_probability():
    """A negative exponent climbs without limit as separation falls."""
    assert cal.probability_ceiling({'a': 13.12, 'b': -2.32}) == 3
    assert cal.probability_ceiling({'a': 1.64, 'b': -1.33}) == 1


def test_ceiling_is_consistent_with_the_curve():
    for a, b in [(13.12, -2.32), (1.64, -1.33), (6.008, -2.542)]:
        ceiling = cal.probability_ceiling({'a': a, 'b': b})
        if ceiling is None:
            continue
        assert cal.power_law(ceiling, a, b) > 1.0
        assert cal.power_law(ceiling + 1, a, b) <= 1.0


def test_no_ceiling_when_the_curve_never_exceeds_one():
    assert cal.probability_ceiling({'a': 0.5, 'b': -1.3}) is None


# --------------------------------------------------------------- end to end

@pytest.fixture(scope='module')
def fixture_fit():
    return cal.calibrate(PDB, DCD, cutoff=1.2, scheme='closest-heavy', all_pairs=True,
                         source='CG fixture, plumbing check')


def test_fit_on_the_fixture_is_pinned(fixture_fit):
    """Golden values. The fixture is a p53 mutant, not a reference linker, so these are
    not a reference curve -- they are here to catch the fit changing."""
    assert fixture_fit['a'] == pytest.approx(13.334, abs=1e-3)
    assert fixture_fit['b'] == pytest.approx(-2.3389, abs=1e-4)
    assert fixture_fit['r_squared'] == pytest.approx(0.99352, abs=1e-5)


def test_fit_records_what_it_describes(fixture_fit):
    assert fixture_fit['cutoff'] == 1.2
    assert fixture_fit['scheme'] == 'closest-heavy'
    assert fixture_fit['all_pairs'] is True
    assert fixture_fit['separation_min'] == 3
    assert fixture_fit['separation_max'] == 60
    assert fixture_fit['n_frames'] == 1000
    assert fixture_fit['n_residues'] == 61
    assert fixture_fit['topology'] == 'A39G_K24Y.pdb'
    assert fixture_fit['source'].startswith('CG fixture')
    assert fixture_fit['fitted'].endswith('Z')


def test_minimum_separation_is_respected():
    fit = cal.calibrate(PDB, DCD, cutoff=1.2, scheme='closest-heavy', all_pairs=True,
                        min_separation=8)
    assert fit['separation_min'] == 8


def test_maximum_separation_is_respected():
    fit = cal.calibrate(PDB, DCD, cutoff=1.2, scheme='closest-heavy', all_pairs=True,
                        max_separation=30)
    assert fit['separation_max'] == 30


def test_a_single_trajectory_is_accepted_as_a_bare_name(fixture_fit):
    assert fixture_fit['n_trajectories'] == 1
    assert fixture_fit['trajectories'] == ['A39G_K24Y.dcd']


# --------------------------------------------------------------- reporting

def test_summary_reports_the_curve_and_its_range(fixture_fit):
    text = cal.summarise(fixture_fit)
    assert 'p_ref(s)' in text
    assert 'R squared' in text
    assert 'separation 3 to 60' in text
    assert '1.2 nm' in text


def test_summary_reports_where_the_curve_exceeds_one(fixture_fit):
    assert 'p_ref above 1' in cal.summarise(fixture_fit)


def test_summary_says_nothing_when_the_curve_stays_below_one():
    separations = np.arange(3, 40)
    fit = cal.fit_power_law(pd.Series(cal.power_law(separations, 0.5, -1.3), index=separations))
    fit.update({'cutoff': 0.8, 'scheme': 'ca', 'n_frames': 10, 'n_residues': 20,
                'n_trajectories': 1, 'separation_min': 3, 'separation_max': 39})
    assert 'p_ref above 1' not in cal.summarise(fit)


def test_saved_fit_round_trips(fixture_fit, tmp_path):
    path = str(tmp_path / 'fit.json')
    cal.save(fixture_fit, path)
    assert cal.load(path) == fixture_fit


def test_saved_fit_is_readable_json(fixture_fit, tmp_path):
    path = str(tmp_path / 'fit.json')
    cal.save(fixture_fit, path)
    with open(path) as handle:
        assert json.load(handle)['a'] == fixture_fit['a']


# --------------------------------------------------------------- command line

def test_command_line_runs_and_writes(tmp_path, capsys):
    out = str(tmp_path / 'cli.json')
    fit = cal.main([PDB, DCD, '--cutoff', '1.2', '--all-pairs', '--out', out])
    assert fit['a'] == pytest.approx(13.334, abs=1e-3)
    assert os.path.exists(out)
    assert 'p_ref(s)' in capsys.readouterr().out


def test_command_line_defaults_the_scheme_to_the_model():
    """--all-pairs implies a bead model, where 'ca' would find nothing."""
    fit = cal.main([PDB, DCD, '--cutoff', '1.2', '--all-pairs'])
    assert fit['scheme'] == 'closest-heavy'


def test_command_line_accepts_a_separation_range(tmp_path):
    fit = cal.main([PDB, DCD, '--cutoff', '1.2', '--all-pairs', '--sep', '5', '25'])
    assert fit['separation_min'] == 5
    assert fit['separation_max'] == 25


# --------------------------------------------------------------- independence

def test_calibration_does_not_import_the_analysis():
    """It exists to produce the constants a profile is given. If it read a profile to get
    them, it would be measuring its own assumption."""
    import inspect
    source = inspect.getsource(cal)
    for name in ('profiles', 'normalization', 'default_function', 'contact_map_generation'):
        assert 'import {0}'.format(name) not in source, name
