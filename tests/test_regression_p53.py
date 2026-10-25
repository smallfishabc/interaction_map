"""End-to-end regression against the August 2023 p53 output.

The reference was produced on a different numpy/libm build. It cannot even regenerate its
own relative_strength column from its own gs_standard column to better than 5.3e-13, and
every formulation of the power law reproduces it equally well, so the continuous columns
are compared with a tolerance far below anything that could move a result. What decides
the figure -- the interaction level and the pair identity -- must match exactly.
"""
import numpy as np
import pytest

from interaction_map import normalization as nl
from interaction_map import profiles

from conftest import requires_p53

RTOL = 1e-11

pytestmark = requires_p53


@pytest.fixture(scope='module')
def replayed(p53_reference):
    """Normalization rerun on the stored contact probabilities."""
    source = p53_reference[['r_1', 'r_2', 'cont_prob', 'distance']].copy()
    return nl.normalization(source, profile='all_atom')


def test_reference_has_the_expected_shape(p53_reference):
    assert len(p53_reference) == 1711
    assert list(p53_reference.columns) == ['r_1', 'r_2', 'cont_prob', 'distance',
                                           'gs_standard', 'relative_strength', 'plot_value']


def test_interaction_levels_are_identical(replayed, p53_reference):
    assert (replayed['plot_value'].values == p53_reference['plot_value'].values).all()


def test_level_counts_are_unchanged(replayed):
    counts = {level: int((replayed['plot_value'] == level).sum())
              for level in (2, 1, 0, -1, -2)}
    assert counts == {2: 60, 1: 129, 0: 625, -1: 166, -2: 731}


@pytest.mark.parametrize('column', ['gs_standard', 'relative_strength'])
def test_continuous_columns_match_within_tolerance(replayed, p53_reference, column):
    got = replayed[column].values
    want = p53_reference[column].values
    rel = np.abs(got - want) / np.maximum(np.abs(want), 1e-300)
    assert rel.max() <= RTOL, '{0}: max relative diff {1:.3e}'.format(column, rel.max())


def test_calibrated_constants_have_not_moved():
    p = profiles.ALL_ATOM
    assert (p.reference_a, p.reference_b) == (1.64, -1.33)
    assert tuple(p.inter_cutoff) == (1.5, 0.5, -1, -2)
    assert p.cutoff == 0.8


def test_the_reference_file_is_internally_inconsistent(p53_reference):
    """Documents why the tolerance above exists, rather than leaving it unexplained.

    Recomputing relative_strength from the file's own two columns does not reproduce the
    file's own third column. The discrepancy predates this branch.
    """
    with np.errstate(divide='ignore'):
        recomputed = np.where(p53_reference['cont_prob'] == 0, 0,
                              np.log(p53_reference['cont_prob'] / p53_reference['gs_standard']))
    worst = np.abs(recomputed - p53_reference['relative_strength'].values).max()
    assert worst > 0
    assert worst < 1e-11
