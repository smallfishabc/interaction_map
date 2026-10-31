"""Schema checks, and the sequence/topology guard."""
import numpy as np
import pandas as pd
import pytest

from ensemble_map import contact_map_generation as cmg
from ensemble_map import normalization as nl
from ensemble_map import schema

from conftest import requires_p53, synthetic_trajectory


def good_contact_table():
    rows = [(1, 5, 0.4, 4), (1, 9, 0.1, 8), (2, 20, 0.02, 18)]
    return pd.DataFrame(rows, columns=['r_1', 'r_2', 'cont_prob', 'distance'])


# --------------------------------------------------------------- schema

def test_a_sound_contact_table_has_no_problems():
    assert schema.validate(good_contact_table(), schema.CONTACT_TABLE) == []


def test_missing_column_is_reported():
    table = good_contact_table().drop(columns=['distance'])
    problems = schema.validate(table, schema.CONTACT_TABLE)
    assert len(problems) == 1
    assert 'distance' in problems[0]


def test_empty_table_is_reported():
    empty = good_contact_table().iloc[0:0]
    problems = schema.validate(empty, schema.CONTACT_TABLE)
    assert any('no rows' in p for p in problems)


def test_probability_out_of_range_is_reported():
    table = good_contact_table()
    table.loc[0, 'cont_prob'] = 1.5
    assert any('between 0 and 1' in p for p in schema.validate(table, schema.CONTACT_TABLE))


def test_inconsistent_separation_is_reported():
    """The failure mode when two tables get concatenated with different numbering."""
    table = good_contact_table()
    table.loc[0, 'distance'] = 99
    assert any('r_2 - r_1' in p for p in schema.validate(table, schema.CONTACT_TABLE))


def test_reversed_pair_is_reported():
    table = good_contact_table()
    table.loc[0, ['r_1', 'r_2']] = [9, 1]
    table.loc[0, 'distance'] = -8
    assert any('greater than' in p for p in schema.validate(table, schema.CONTACT_TABLE))


def test_missing_values_are_reported():
    table = good_contact_table()
    table.loc[1, 'cont_prob'] = np.nan
    assert any('missing values' in p for p in schema.validate(table, schema.CONTACT_TABLE))


def test_interaction_table_requires_the_derived_columns():
    problems = schema.validate(good_contact_table(), schema.INTERACTION_TABLE)
    assert any('gs_standard' in p for p in problems)


def test_normalized_table_satisfies_the_interaction_schema():
    out = nl.normalization(good_contact_table(), profile='all_atom')
    assert schema.validate(out, schema.INTERACTION_TABLE) == []


def test_unknown_interaction_level_is_reported():
    out = nl.normalization(good_contact_table(), profile='all_atom')
    out.loc[0, 'plot_value'] = 7
    assert any('levels' in p for p in schema.validate(out, schema.INTERACTION_TABLE))


def test_require_raises_on_a_bad_table():
    with pytest.raises(ValueError):
        schema.require(good_contact_table().drop(columns=['r_1']), schema.CONTACT_TABLE)


def test_require_returns_the_frame_when_sound():
    table = good_contact_table()
    assert schema.require(table, schema.CONTACT_TABLE) is table


def test_warn_if_invalid_prints_and_reports(capsys):
    assert schema.warn_if_invalid(good_contact_table(), schema.CONTACT_TABLE) is True
    assert capsys.readouterr().out == ''
    assert schema.warn_if_invalid(good_contact_table().drop(columns=['r_1']),
                                  schema.CONTACT_TABLE) is False
    assert 'WARNING' in capsys.readouterr().out


# --------------------------------------------------------------- sequence guard

def test_topology_sequence_is_empty_for_bead_models():
    """CG beads are named with one letter codes mdtraj does not treat as residues."""
    traj = synthetic_trajectory(n_residues=6, n_frames=2, atoms_per_residue=1)
    # The synthetic topology uses ALA, which does have a code.
    assert cmg.topology_sequence(traj) == 'A' * 6


def test_sequence_check_passes_when_they_agree():
    traj = synthetic_trajectory(n_residues=6, n_frames=2, atoms_per_residue=1)
    assert cmg.check_sequence_against_topology(traj, 'AAAAAA') is True


def test_sequence_check_reports_a_length_mismatch(capsys):
    traj = synthetic_trajectory(n_residues=6, n_frames=2, atoms_per_residue=1)
    assert cmg.check_sequence_against_topology(traj, 'AAAA') is False
    assert 'residues' in capsys.readouterr().out


def test_sequence_check_reports_the_first_differing_residue(capsys):
    traj = synthetic_trajectory(n_residues=6, n_frames=2, atoms_per_residue=1)
    assert cmg.check_sequence_against_topology(traj, 'AAWAAA') is False
    out = capsys.readouterr().out
    assert 'residue 3' in out


def test_sequence_check_is_skipped_when_the_topology_cannot_be_read():
    """Returns None rather than a false alarm."""
    import mdtraj as md
    top = md.Topology()
    chain = top.add_chain()
    for _ in range(4):
        residue = top.add_residue('XXX', chain)
        top.add_atom('CA', md.element.carbon, residue)
    traj = md.Trajectory(np.zeros((2, 4, 3), dtype=np.float32), top)
    assert cmg.check_sequence_against_topology(traj, 'ABCD') is None


@requires_p53
def test_p53_topology_sequence_matches_seq_fasta(p53_sequence):
    """The capped all-atom topology yields exactly the sequence in seq.fasta."""
    import os
    import mdtraj as md
    pdb = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'All_atom_p53_interaction', 'p53-summary', 'BB', 'S_0', '__START_0.pdb')
    traj = md.load(pdb)
    sliced = traj.atom_slice(traj.top.select('protein'))
    assert cmg.topology_sequence(sliced) == p53_sequence
    assert cmg.check_sequence_against_topology(sliced, p53_sequence) is True


# --------------------------------------------------------------- sequence reading

def test_sequence_cleaning_matches_the_original_implementation():
    """The original walked backwards with a regex; rstrip must give the same answer."""
    import re
    from ensemble_map import readpath

    def original(seq):
        while re.search(r'\s', seq[-1]) is not None:
            seq = seq[:-1]
        return seq

    for case in ['MEEPQ', 'MEEPQ\n', 'MEEPQ\r\n', 'MEEPQ  \n\n', '  MEEPQ\n', 'MEEPQ\tX\n']:
        assert readpath._clean(case) == original(case), case


def test_an_empty_sequence_file_still_fails():
    """The original ran off the front of the string and raised IndexError. Failing is
    still right; only the message changed."""
    from ensemble_map import readpath
    for case in ['', '   ', '\n']:
        with pytest.raises(ValueError):
            readpath._clean(case)
