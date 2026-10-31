"""Contact map generation: pair selection, residue numbering, and the numbering guard."""
import os

import numpy as np
import pandas as pd
import pytest

from ensemble_map import contact_map_generation as cmg
from ensemble_map import profiles

from conftest import requires_p53, synthetic_trajectory

SEQ20 = 'MEEPQSDPSVEPPLSQETFS'


# --------------------------------------------------------------- numbering guard

def test_check_residue_numbering_accepts_a_matching_table(capsys):
    table = pd.DataFrame({'r_1': [1, 2], 'r_2': [19, 20]})
    assert cmg.check_residue_numbering(table, SEQ20) is True
    assert capsys.readouterr().out == ''


def test_check_residue_numbering_rejects_a_shifted_table(capsys):
    """A zero-based table against a capped sequence is the failure this guards."""
    table = pd.DataFrame({'r_1': [0, 1], 'r_2': [18, 19]})
    assert cmg.check_residue_numbering(table, SEQ20) is False
    assert 'WARNING' in capsys.readouterr().out


def test_check_residue_numbering_rejects_a_truncated_table():
    table = pd.DataFrame({'r_1': [1], 'r_2': [15]})
    assert cmg.check_residue_numbering(table, SEQ20) is False


# --------------------------------------------------------------- pair selection

def test_all_atom_selection_drops_short_separations():
    """mdtraj's contacts='all' excludes |i - j| < 3, which the all-atom path relies on."""
    traj = synthetic_trajectory(n_residues=12, n_frames=4, atoms_per_residue=1)
    contact = cmg.ContactProbData('t', 0.8, traj, profile='all_atom')
    assert contact.contact['distance'].min() >= 3


def test_cg_selection_keeps_every_pair():
    traj = synthetic_trajectory(n_residues=12, n_frames=4, atoms_per_residue=1)
    contact = cmg.ContactProbData('t', 1.2, traj, profile='cg')
    n = traj.n_residues
    assert len(contact.contact) == n * (n - 1) // 2
    assert contact.contact['distance'].min() == 1


def test_unknown_pair_selection_is_rejected():
    traj = synthetic_trajectory(n_residues=8, n_frames=2, atoms_per_residue=1)
    broken = profiles.Profile(**dict(profiles.ALL_ATOM.__dict__, pair_selection='magic'))
    with pytest.raises(ValueError):
        cmg.ContactProbData('t', 0.8, traj, profile=broken)


# --------------------------------------------------------------- numbering offset

def test_cg_offset_makes_the_first_residue_one():
    """CALVADOS topologies have no caps, so the raw index needs shifting."""
    traj = synthetic_trajectory(n_residues=10, n_frames=3, atoms_per_residue=1)
    contact = cmg.ContactProbData('t', 1.2, traj, profile='cg')
    assert int(contact.contact['r_1'].min()) == 1
    assert int(contact.contact['r_2'].max()) == traj.n_residues


def test_all_atom_offset_is_not_applied():
    """Uncapped input therefore starts at 0, which is what the guard is for."""
    traj = synthetic_trajectory(n_residues=12, n_frames=3, atoms_per_residue=1)
    contact = cmg.ContactProbData('t', 0.8, traj, profile='all_atom')
    assert int(contact.contact['r_1'].min()) == 0
    assert cmg.check_residue_numbering(contact.contact, 'A' * traj.n_residues) is False


def test_distance_is_the_sequence_separation():
    traj = synthetic_trajectory(n_residues=12, n_frames=3, atoms_per_residue=1)
    contact = cmg.ContactProbData('t', 0.8, traj, profile='all_atom')
    table = contact.contact
    assert (table['distance'] == table['r_2'] - table['r_1']).all()


def test_indices_are_integers():
    traj = synthetic_trajectory(n_residues=10, n_frames=3, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg').contact
    assert table['r_1'].dtype.kind == 'i'
    assert table['r_2'].dtype.kind == 'i'


def test_contact_probability_is_a_probability():
    traj = synthetic_trajectory(n_residues=12, n_frames=6, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg').contact
    assert table['cont_prob'].min() >= 0.0
    assert table['cont_prob'].max() <= 1.0


def test_contact_probability_is_quantised_by_frame_count():
    """Each value is a count of frames over the total, so it lands on a 1/n_frames grid."""
    n_frames = 8
    traj = synthetic_trajectory(n_residues=10, n_frames=n_frames, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg').contact
    counts = table['cont_prob'].values * n_frames
    assert np.allclose(counts, np.round(counts))


def test_cutoff_is_monotonic():
    """A larger cutoff can only add contacts, never remove them."""
    traj = synthetic_trajectory(n_residues=12, n_frames=6, atoms_per_residue=1)
    tight = cmg.ContactProbData('t', 0.6, traj, profile='cg').contact
    loose = cmg.ContactProbData('t', 1.5, traj, profile='cg').contact
    assert (loose['cont_prob'].values >= tight['cont_prob'].values).all()


# --------------------------------------------------------------- residue identity

def test_residue_letter_uses_the_standard_code():
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue('TRP', chain)
    assert cmg.residue_letter(residue) == 'W'


def test_residue_letter_falls_back_to_a_one_character_name():
    """Coarse-grained beads are named with the one letter code already."""
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    assert cmg.residue_letter(top.add_residue('E', chain)) == 'E'


def test_residue_letter_is_x_for_caps_and_unknowns():
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    assert cmg.residue_letter(top.add_residue('ACE', chain)) == 'X'


def test_contact_table_carries_residue_letters_and_chains():
    traj = synthetic_trajectory(n_residues=12, n_frames=3, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg').contact
    for column in ('r_1_aa', 'r_2_aa', 'r_1_chain', 'r_2_chain'):
        assert column in table.columns
    assert set(table['r_1_aa']) == {'A'}
    assert set(table['r_1_chain']) == {0}


@requires_p53
def test_residue_letters_agree_with_the_sequence(p53_reference, p53_sequence):
    """An independent check on the 1-based numbering: the letter at r_1 must be the
    residue the sequence has at that position."""
    import mdtraj as md
    pdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'All_atom_p53_interaction', 'p53-summary', 'BB', 'S_0', '__START_0.pdb')
    traj = md.load(pdb)
    sliced = traj.atom_slice(traj.top.select('protein'))
    table = cmg.ContactProbData('t', 0.8, sliced, profile='all_atom').contact
    for _, row in table.head(50).iterrows():
        assert row['r_1_aa'] == p53_sequence[int(row['r_1']) - 1]
        assert row['r_2_aa'] == p53_sequence[int(row['r_2']) - 1]


# --------------------------------------------------------------- repeats

def test_trajectory_file_list_expands_a_repeat_count():
    assert cmg.trajectory_file_list(3) == ['__traj_0.xtc', '__traj_1.xtc', '__traj_2.xtc']


def test_trajectory_file_list_wraps_a_single_name():
    assert cmg.trajectory_file_list('one.dcd') == ['one.dcd']


def test_trajectory_file_list_passes_a_list_through():
    assert cmg.trajectory_file_list(['a.xtc', 'b.xtc']) == ['a.xtc', 'b.xtc']


def test_single_repeat_reports_no_spread():
    traj = synthetic_trajectory(n_residues=10, n_frames=6, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg',
                                replica_frames=[6]).contact
    assert set(table['n_rep']) == {1}
    assert table['p_sem'].isna().all()
    assert set(table['n_frames']) == {6}


def test_multiple_repeats_report_a_spread():
    traj = synthetic_trajectory(n_residues=10, n_frames=12, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg',
                                replica_frames=[4, 4, 4]).contact
    assert set(table['n_rep']) == {3}
    assert set(table['n_frames']) == {12}
    assert table['p_sem'].notna().all()
    assert (table['p_sem'] >= 0).all()


def test_spread_is_skipped_when_frame_counts_do_not_add_up():
    """Rather than silently slicing the wrong frames."""
    traj = synthetic_trajectory(n_residues=10, n_frames=12, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg',
                                replica_frames=[4, 4]).contact
    assert table['p_sem'].isna().all()


def test_contact_probability_is_unaffected_by_the_spread_calculation():
    """cont_prob must stay the average over every frame, however the repeats are named."""
    traj = synthetic_trajectory(n_residues=10, n_frames=12, atoms_per_residue=1)
    without = cmg.ContactProbData('t', 1.2, traj, profile='cg').contact
    with_reps = cmg.ContactProbData('t', 1.2, traj, profile='cg',
                                    replica_frames=[4, 4, 4]).contact
    assert (without['cont_prob'].values == with_reps['cont_prob'].values).all()


def test_sem_matches_the_per_replica_definition():
    traj = synthetic_trajectory(n_residues=10, n_frames=9, atoms_per_residue=1)
    table = cmg.ContactProbData('t', 1.2, traj, profile='cg',
                                replica_frames=[3, 3, 3]).contact
    # Recompute independently from the trajectory thirds.
    import mdtraj as md
    indices = traj.top.select_pairs('all', 'all')
    cont, _pairs = md.compute_contacts(traj, indices)
    hit = cont < 1.2
    means = np.vstack([hit[0:3].mean(0), hit[3:6].mean(0), hit[6:9].mean(0)])
    expected = means.std(axis=0, ddof=1) / np.sqrt(3)
    assert np.allclose(table['p_sem'].values, expected)


# --------------------------------------------------------------- the class and loader

def test_profile_is_stored_on_the_object():
    traj = synthetic_trajectory(n_residues=8, n_frames=2, atoms_per_residue=1)
    contact = cmg.ContactProbData('t', 1.2, traj, profile='cg')
    assert contact.profile is profiles.CG


def test_constructing_writes_the_csv(tmp_path):
    """A documented side effect: the object saves itself on construction."""
    os.chdir(str(tmp_path))
    traj = synthetic_trajectory(n_residues=8, n_frames=2, atoms_per_residue=1)
    cmg.ContactProbData('demo', 1.2, traj, profile='cg')
    assert (tmp_path / 'demo_1.2_contact_df_1201.csv').exists()


def test_saved_csv_round_trips_by_value(tmp_path):
    os.chdir(str(tmp_path))
    traj = synthetic_trajectory(n_residues=8, n_frames=2, atoms_per_residue=1)
    written = cmg.ContactProbData('demo', 1.2, traj, profile='cg').contact
    reloaded = cmg.ContactProbData('demo', 1.2, None, read_from_file=1, profile='cg').contact
    assert list(written.columns) == list(reloaded.columns)
    pd.testing.assert_frame_equal(written, reloaded, check_dtype=False)


def test_residue_indices_are_int64_on_every_platform(tmp_path):
    """astype('int') was the C long: int32 on Windows, int64 elsewhere, and always int64
    after a CSV round trip. Pinned to int64 so a computed table and a reloaded one agree.
    """
    os.chdir(str(tmp_path))
    traj = synthetic_trajectory(n_residues=8, n_frames=2, atoms_per_residue=1)
    written = cmg.ContactProbData('demo', 1.2, traj, profile='cg').contact
    reloaded = cmg.ContactProbData('demo', 1.2, None, read_from_file=1, profile='cg').contact
    assert written['r_1'].dtype == np.int64
    assert written['r_2'].dtype == np.int64
    assert reloaded['r_1'].dtype == written['r_1'].dtype
    assert (written['r_1'].values == reloaded['r_1'].values).all()


def test_generate_contact_takes_the_cutoff_from_the_profile(tmp_path, monkeypatch):
    """0.8 nm for all-atom, 1.2 nm for CG, unless given explicitly."""
    os.chdir(str(tmp_path))
    traj = synthetic_trajectory(n_residues=10, n_frames=3, atoms_per_residue=1)
    monkeypatch.setattr(cmg, 'load_traj_protein', lambda *a, **k: (traj, traj, traj.n_frames))
    assert cmg.generate_contact('x', profile='all_atom').cutoff == 0.8
    assert cmg.generate_contact('x', profile='cg').cutoff == 1.2
    assert cmg.generate_contact('x', cutoff=0.95, profile='cg').cutoff == 0.95


def test_load_traj_protein_skips_slicing_for_cg(tmp_path, monkeypatch):
    traj = synthetic_trajectory(n_residues=10, n_frames=3, atoms_per_residue=1)
    monkeypatch.setattr(cmg, 'load_xtc', lambda xtc, pdb: traj)
    full, sliced, frames = cmg.load_traj_protein('t.pdb', 't.xtc', profile='cg')
    assert sliced is full
    assert frames == traj.n_frames


def test_load_traj_protein_slices_for_all_atom(tmp_path, monkeypatch):
    traj = synthetic_trajectory(n_residues=10, n_frames=3, atoms_per_residue=2)
    monkeypatch.setattr(cmg, 'load_xtc', lambda xtc, pdb: traj)
    full, sliced, _frames = cmg.load_traj_protein('t.pdb', 't.xtc', profile='all_atom')
    assert sliced is not full


def test_repeat_count_expands_to_a_file_list(monkeypatch):
    """An integer means '__traj_0.xtc' .. '__traj_{n-1}.xtc', in order."""
    seen = {}
    traj = synthetic_trajectory(n_residues=6, n_frames=2, atoms_per_residue=1)

    def fake_load(xtc, pdb):
        seen['xtc'] = xtc
        return traj

    monkeypatch.setattr(cmg, 'load_xtc', fake_load)
    cmg.load_traj_protein('top.pdb', 5, profile='cg')
    assert seen['xtc'] == ['__traj_0.xtc', '__traj_1.xtc', '__traj_2.xtc',
                           '__traj_3.xtc', '__traj_4.xtc']


# --------------------------------------------------------------- against real data

@requires_p53
def test_p53_topology_is_capped(p53_sequence):
    """The fact the whole all-atom numbering rests on."""
    import mdtraj as md
    pdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'All_atom_p53_interaction', 'p53-summary', 'BB', 'S_0', '__START_0.pdb')
    traj = md.load(pdb)
    sliced = traj.atom_slice(traj.top.select('protein'))
    names = [r.name for r in sliced.top.residues]
    assert names[0] == 'ACE'
    assert names[-1] == 'NME'
    # cap + sequence + cap
    assert sliced.n_residues == len(p53_sequence) + 2


@requires_p53
def test_p53_reference_numbering_spans_the_sequence(p53_reference, p53_sequence):
    assert int(p53_reference['r_1'].min()) == 1
    assert int(p53_reference['r_2'].max()) == len(p53_sequence)
    assert cmg.check_residue_numbering(p53_reference, p53_sequence) is True
