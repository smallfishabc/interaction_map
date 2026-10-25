"""The conditional analysis: frame selection, the split, and the ensemble comparison.

The selector is the most consequential knob here, so most of these tests are about which
residue pairs it looks at and how the numbering is translated.
"""
import numpy as np
import pytest

from interaction_map import profiles
from interaction_map import search_chunk as sc

from conftest import synthetic_trajectory


# --------------------------------------------------------------- the selector

def test_default_window_is_the_published_three_by_three():
    selector = sc.ChunkContact(13, 24)
    pairs = selector.pairs()
    assert len(pairs) == 9
    assert sorted(pairs) == sorted([[x, y] for x in (12, 13, 14) for y in (23, 24, 25)])


def test_chunk_false_narrows_to_the_single_pair():
    assert sc.ChunkContact(13, 24, chunk=False).pairs() == [[13, 24]]


def test_window_two_gives_a_five_by_five():
    assert len(sc.ChunkContact(13, 24, window=2).pairs()) == 25


def test_window_zero_is_the_single_pair():
    assert sc.ChunkContact(13, 24, window=0).pairs() == [[13, 24]]


def test_all_atom_positions_pass_through_unshifted():
    """The ACE cap already makes mdtraj's index the 1-based position."""
    selector = sc.ChunkContact(13, 24, profile='all_atom')
    assert selector.profile.residue_index_offset == 0
    assert [13, 24] in selector.pairs()


def test_cg_positions_are_shifted_back_to_raw_indices():
    """CG stores position = raw index + 1, so the selector must undo that.

    The window still spans one either side, so neighbouring pairs remain inside it. What
    matters is that it is centred on the shifted pair.
    """
    pairs = sc.ChunkContact(13, 24, profile='cg').pairs()
    assert sorted(pairs) == sorted([[x, y] for x in (11, 12, 13) for y in (22, 23, 24)])
    # Same request under all_atom is centred one higher.
    all_atom = sc.ChunkContact(13, 24, profile='all_atom').pairs()
    assert sorted(all_atom) == sorted([[x, y] for x in (12, 13, 14) for y in (23, 24, 25)])


def test_the_two_profiles_centre_the_window_differently():
    cg_centre = sc.ChunkContact(13, 24, window=0, profile='cg').pairs()
    aa_centre = sc.ChunkContact(13, 24, window=0, profile='all_atom').pairs()
    assert cg_centre == [[12, 23]]
    assert aa_centre == [[13, 24]]


def test_cutoff_comes_from_the_profile():
    assert sc.ChunkContact(5, 15, profile='all_atom').cutoff == 0.8
    assert sc.ChunkContact(5, 15, profile='cg').cutoff == 1.2
    assert sc.ChunkContact(5, 15, profile='cg', cutoff=0.9).cutoff == 0.9


def test_describe_records_what_was_asked_for():
    described = sc.ChunkContact(13, 24, window=1, profile='all_atom').describe()
    assert described == {'r1': 13, 'r2': 24, 'window': 1, 'cutoff': 0.8,
                         'chunk': True, 'profile': 'all_atom', 'n_pairs': 9}


def test_describe_reports_window_zero_when_not_chunking():
    assert sc.ChunkContact(13, 24, chunk=False).describe()['window'] == 0


# --------------------------------------------------------------- frame selection

def test_frames_where_returns_one_boolean_per_frame():
    traj = synthetic_trajectory(n_residues=20, n_frames=11, atoms_per_residue=1)
    mask = sc.frames_where(traj, sc.ChunkContact(5, 15, profile='cg'))
    assert mask.dtype == bool
    assert mask.shape == (11,)


def test_a_generous_cutoff_selects_every_frame():
    traj = synthetic_trajectory(n_residues=20, n_frames=8, atoms_per_residue=1)
    mask = sc.frames_where(traj, sc.ChunkContact(5, 15, cutoff=1e6, profile='cg'))
    assert mask.all()


def test_a_tiny_cutoff_selects_no_frame():
    traj = synthetic_trajectory(n_residues=20, n_frames=8, atoms_per_residue=1)
    mask = sc.frames_where(traj, sc.ChunkContact(5, 15, cutoff=1e-9, profile='cg'))
    assert not mask.any()


def test_a_wider_window_can_only_select_more_frames():
    traj = synthetic_trajectory(n_residues=24, n_frames=30, atoms_per_residue=1, seed=3)
    narrow = sc.frames_where(traj, sc.ChunkContact(6, 18, window=0, cutoff=1.2, profile='cg'))
    wide = sc.frames_where(traj, sc.ChunkContact(6, 18, window=2, cutoff=1.2, profile='cg'))
    assert (wide | narrow == wide).all()


def test_inter_frames_matches_the_selector_it_builds():
    traj = synthetic_trajectory(n_residues=20, n_frames=9, atoms_per_residue=1)
    direct = sc.frames_where(traj, sc.ChunkContact(5, 15, cutoff=1.2, profile='cg'))
    through = sc.inter_frames(5, 15, traj, cutoff=1.2, profile='cg')
    assert (direct == through).all()


# --------------------------------------------------------------- the split

def test_split_maps_cover_every_frame_between_them():
    traj = synthetic_trajectory(n_residues=20, n_frames=20, atoms_per_residue=1, seed=5)
    mask = sc.frames_where(traj, sc.ChunkContact(5, 15, cutoff=1.2, profile='cg'))
    if mask.all() or not mask.any():
        pytest.skip('the synthetic trajectory did not produce a two sided split')
    with_it, without_it = sc.new_interaction_map('t', traj, mask, profile='cg')
    assert int(with_it.contact['n_frames'].iloc[0]) == int(mask.sum())
    assert int(without_it.contact['n_frames'].iloc[0]) == int((~mask).sum())
    assert (int(mask.sum()) + int((~mask).sum())) == traj.n_frames


def test_split_maps_are_named_apart():
    traj = synthetic_trajectory(n_residues=20, n_frames=12, atoms_per_residue=1, seed=5)
    mask = np.zeros(12, dtype=bool)
    mask[:6] = True
    with_it, without_it = sc.new_interaction_map('demo', traj, mask, profile='cg')
    assert with_it.name == 'demowith'
    assert without_it.name == 'demowithout'


# --------------------------------------------------------------- observables

def test_rg_matches_mdtraj():
    import mdtraj as md
    import statistics as st
    traj = synthetic_trajectory(n_residues=12, n_frames=7, atoms_per_residue=1)
    assert sc.calc_Rg(traj) == st.mean(md.compute_rg(traj))


def test_end_to_end_distance_uses_the_first_and_last_alpha_carbon():
    """statistics.mean is kept rather than numpy's, so compare against the same thing."""
    import mdtraj as md
    import statistics as st
    traj = synthetic_trajectory(n_residues=12, n_frames=5, atoms_per_residue=1)
    alpha = traj.topology.select_atom_indices(selection='alpha')
    distances = md.compute_distances(traj, [[alpha[0], alpha[-1]]])
    expected = st.mean([float(d[0]) for d in distances])
    assert sc.calc_Re(traj, traj) == expected


# The secondary structure and hydrogen bond observables need a real all-atom topology with
# bonds, which the synthetic trajectory does not have. Stub them out to test the part that
# changed, which is the assembly of the result and the frame counts.
@pytest.fixture
def light_observables(monkeypatch):
    monkeypatch.setattr(sc, 'calc_Heli', lambda t, length: 0.1)
    monkeypatch.setattr(sc, 'calc_Beta', lambda t, length: 0.2)
    monkeypatch.setattr(sc, 'calc_HB', lambda r, j, length: 0.3)


def test_ensemble_reports_the_frame_count(light_observables):
    """The original output did not say how big each side of the split was."""
    traj = synthetic_trajectory(n_residues=12, n_frames=7, atoms_per_residue=1)
    values = sc.calculate_ensemble(traj, traj, 12)
    assert values['Frames'] == 7
    assert set(values) == {'Rg', 'Re', 'Heli', 'Beta', 'HB', 'Frames'}


def test_compare_ensemble_splits_the_frames_between_the_two_sides(light_observables):
    traj = synthetic_trajectory(n_residues=12, n_frames=10, atoms_per_residue=1)
    mask = np.zeros(10, dtype=bool)
    mask[:4] = True
    result = sc.compare_ensemble(traj, traj, 12, mask)
    assert result['interaction']['Frames'] == 4
    assert result['non_interaction']['Frames'] == 6
    assert set(result) == {'interaction', 'non_interaction'}


def test_compare_ensemble_measures_the_two_sides_separately(light_observables):
    """Rg is a real calculation here, so the two sides must not come out identical."""
    traj = synthetic_trajectory(n_residues=12, n_frames=10, atoms_per_residue=1, seed=7)
    mask = np.zeros(10, dtype=bool)
    mask[:4] = True
    result = sc.compare_ensemble(traj, traj, 12, mask)
    assert result['interaction']['Rg'] != result['non_interaction']['Rg']
