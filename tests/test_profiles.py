"""The two profiles carry calibrated constants. These tests pin them.

If one of these fails, a published number has moved.
"""
import pytest

from interaction_map import profiles


def test_all_atom_constants():
    p = profiles.ALL_ATOM
    assert p.name == 'all_atom'
    assert p.cutoff == 0.8
    assert (p.reference_a, p.reference_b) == (1.64, -1.33)
    assert tuple(p.inter_cutoff) == (1.5, 0.5, -1, -2)
    assert tuple(p.value_list) == (2, 1, 0, -1, -2)
    assert p.pair_selection == 'mdtraj_ca'
    assert p.slice_protein is True
    assert p.residue_index_offset == 0
    assert p.min_sep == 3
    assert p.line_width_scale == 2.0
    assert p.prob_floor == 0.0


def test_cg_constants():
    p = profiles.CG
    assert p.name == 'cg'
    assert p.cutoff == 1.2
    assert (p.reference_a, p.reference_b) == (13.12, -2.32)
    assert tuple(p.inter_cutoff) == (2, 1, -1, -2)
    assert tuple(p.value_list) == (2, 1, 0, -1, -2)
    assert p.pair_selection == 'all_pairs'
    assert p.slice_protein is False
    assert p.residue_index_offset == 1
    assert p.min_sep == 5
    assert p.line_width_scale == 4.0
    assert p.prob_floor == 0.0


def test_the_two_profiles_differ_where_they_should():
    aa, cg = profiles.ALL_ATOM, profiles.CG
    for field in ('cutoff', 'pair_selection', 'slice_protein', 'residue_index_offset',
                  'reference_a', 'reference_b', 'inter_cutoff', 'min_sep',
                  'line_width_scale'):
        assert getattr(aa, field) != getattr(cg, field), field


def test_get_profile_accepts_name_instance_and_none():
    assert profiles.get_profile('all_atom') is profiles.ALL_ATOM
    assert profiles.get_profile('cg') is profiles.CG
    assert profiles.get_profile(profiles.CG) is profiles.CG
    assert profiles.get_profile(None) is profiles.ALL_ATOM


def test_get_profile_rejects_unknown_name():
    with pytest.raises(ValueError) as excinfo:
        profiles.get_profile('coarse')
    # The message should list what is available, not just complain.
    assert 'all_atom' in str(excinfo.value)
    assert 'cg' in str(excinfo.value)


def test_profiles_are_immutable():
    """Frozen, so one analysis cannot silently retune another."""
    with pytest.raises(Exception):
        profiles.ALL_ATOM.cutoff = 1.0


def test_registry_matches_the_instances():
    assert profiles.PROFILES == {'all_atom': profiles.ALL_ATOM, 'cg': profiles.CG}
    for name, profile in profiles.PROFILES.items():
        assert profile.name == name
