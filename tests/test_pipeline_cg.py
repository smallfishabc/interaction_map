"""End to end run on the committed coarse-grained fixture.

A39G_K24Y is one CALVADOS p53 double mutant: 61 beads, 1000 frames, 805 KB. Small enough
to live in the repository, so unlike the tests that need All_atom_p53_interaction these
run anywhere. The same files are what TUTORIAL.md walks through.

The pinned numbers below are golden values. If one of them moves, the analysis changed.
"""
import os

import numpy as np
import pytest

import contact_map_generation as cmg
import interaction_plot as ip
import normalization as nl
import profiles
import readpath
import schema
import search_chunk as sc

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')
NAME = 'A39G_K24Y'
PDB = os.path.join(FIXTURES, NAME + '.pdb')
DCD = os.path.join(FIXTURES, NAME + '.dcd')
FASTA = os.path.join(FIXTURES, NAME + '.fasta')

N_RESIDUES = 61
N_FRAMES = 1000
N_PAIRS = N_RESIDUES * (N_RESIDUES - 1) // 2          # 1830, CG keeps every pair


@pytest.fixture(scope='module')
def sequence():
    return readpath.read_sequence_file(FASTA)


@pytest.fixture(scope='module')
def trajectory():
    import mdtraj as md
    return md.load(DCD, top=PDB)


@pytest.fixture(scope='module')
def contacts(tmp_path_factory):
    # Building a contact map writes a csv next to the working directory, so run in tmp.
    os.chdir(str(tmp_path_factory.mktemp('cg')))
    return cmg.generate_contact(NAME, PDB, DCD, profile='cg').contact


@pytest.fixture(scope='module')
def interaction(contacts):
    return nl.normalization(contacts.copy(), profile='cg')


# --------------------------------------------------------------- the fixture itself

def test_fixture_files_are_present():
    for path in (PDB, DCD, FASTA):
        assert os.path.exists(path), path


def test_fixture_shape(trajectory, sequence):
    assert trajectory.n_frames == N_FRAMES
    assert trajectory.n_residues == N_RESIDUES
    assert trajectory.n_atoms == N_RESIDUES           # one bead per residue
    assert len(sequence) == N_RESIDUES


def test_sequence_is_found_without_being_named(sequence):
    """These folders carry <name>.fasta rather than seq.txt."""
    found, source = readpath.find_sequence(FIXTURES, NAME)
    assert found == sequence
    assert source.endswith(NAME + '.fasta')


# --------------------------------------------------------------- contact map

def test_contact_map_covers_every_pair(contacts):
    assert len(contacts) == N_PAIRS
    assert int(contacts['r_1'].min()) == 1
    assert int(contacts['r_2'].max()) == N_RESIDUES


def test_contact_map_satisfies_its_schema(contacts):
    assert schema.validate(contacts, schema.CONTACT_TABLE) == []


def test_numbering_lines_up_with_the_sequence(contacts, sequence):
    assert cmg.check_residue_numbering(contacts, sequence) is True


def test_frame_count_is_recorded(contacts):
    assert set(contacts['n_frames']) == {N_FRAMES}


def test_single_trajectory_reports_no_spread(contacts):
    """One file means one repeat, so there is no between-repeat spread to report."""
    assert set(contacts['n_rep']) == {1}
    assert contacts['p_sem'].isna().all()


def test_contact_probability_total_is_pinned(contacts):
    """Golden value. Catches any change to how contacts are counted."""
    assert float(contacts['cont_prob'].sum()) == pytest.approx(272.213, abs=1e-9)


# --------------------------------------------------------------- normalization

def test_interaction_levels_are_pinned(interaction):
    counts = {int(k): int(v) for k, v in interaction['plot_value'].value_counts().items()}
    assert counts == {0: 1715, -2: 60, -1: 33, 1: 22}


def test_the_sixty_strong_repulsive_pairs_are_the_bonded_neighbours(interaction):
    """Documented coarse-grained behaviour, kept deliberately.

    The CG reference curve returns a value above 1 below separation 4, so bonded
    neighbours, which are always in contact, come out strongly repulsive. min_sep keeps
    them out of the figure.
    """
    strong = interaction[interaction['plot_value'] == -2]
    assert len(strong) == 60
    assert set(strong['distance']) == {1}
    assert (strong['cont_prob'] == 1.0).all()


def test_nothing_strongly_attractive_in_this_mutant(interaction):
    assert int((interaction['plot_value'] == 2).sum()) == 0


# --------------------------------------------------------------- figure

def test_map_renders_and_writes_both_formats(interaction, sequence, tmp_path):
    import matplotlib.pyplot as plt
    os.chdir(str(tmp_path))
    fig, ax = ip.interaction_map(sequence, len(sequence), interaction, NAME,
                                 profile='cg', show=False, save=True)
    assert (tmp_path / (NAME + '.png')).exists()
    assert (tmp_path / (NAME + '.svg')).exists()
    plt.close(fig)


def test_only_pairs_past_min_sep_are_drawn(interaction, sequence):
    """The 60 bonded pairs must not appear, which is what min_sep = 5 is for."""
    layout = ip.create_network(sequence, len(sequence))
    drawn = []

    class Recorder(object):
        def annotate(self, _t, xy=None, xytext=None, arrowprops=None, **kw):
            drawn.append((xy, xytext))

    for level in (2, 1, -1, -2):
        ip.interaction_plotting(interaction, layout, Recorder(), level, profile='cg')
    separations = [round(b[0] - a[0] - 0.4) for a, b in drawn]
    assert separations, 'nothing was drawn at all'
    assert min(separations) >= profiles.CG.min_sep


# --------------------------------------------------------------- conditional analysis

def test_selector_works_on_a_coarse_grained_topology(trajectory):
    """CALVADOS beads are not named CA, so the contact scheme has to come from the profile.

    Asking for 'ca' here makes mdtraj discard every pair and return an empty array.
    """
    assert profiles.CG.contact_scheme == 'closest-heavy'
    assert profiles.ALL_ATOM.contact_scheme == 'ca'
    mask = sc.ChunkContact(20, 40, profile='cg').mask(trajectory)
    assert mask.shape == (N_FRAMES,)
    assert 0 < int(mask.sum()) < N_FRAMES


@pytest.mark.parametrize('pair,expected', [((20, 40), 46), ((30, 50), 34)])
def test_single_interaction_frame_counts_are_pinned(trajectory, pair, expected):
    mask = sc.ChunkContact(pair[0], pair[1], profile='cg').mask(trajectory)
    assert int(mask.sum()) == expected


def test_any_is_the_union_and_all_is_the_intersection(trajectory):
    a = sc.ChunkContact(20, 40, profile='cg').mask(trajectory)
    b = sc.ChunkContact(30, 50, profile='cg').mask(trajectory)
    any_mask = sc.interaction_set([(20, 40), (30, 50)], mode='any', profile='cg').mask(trajectory)
    all_mask = sc.interaction_set([(20, 40), (30, 50)], mode='all', profile='cg').mask(trajectory)
    assert (any_mask == (a | b)).all()
    assert (all_mask == (a & b)).all()
    # and the counts add up the way a union and an intersection must
    assert int(any_mask.sum()) == int(a.sum()) + int(b.sum()) - int(all_mask.sum())


def test_combined_frame_counts_are_pinned(trajectory):
    any_mask = sc.interaction_set([(20, 40), (30, 50)], mode='any', profile='cg').mask(trajectory)
    all_mask = sc.interaction_set([(20, 40), (30, 50)], mode='all', profile='cg').mask(trajectory)
    assert int(any_mask.sum()) == 78
    assert int(all_mask.sum()) == 2


def test_all_is_never_larger_than_any(trajectory):
    pairs = [(20, 40), (30, 50), (10, 25)]
    any_mask = sc.interaction_set(pairs, mode='any', profile='cg').mask(trajectory)
    all_mask = sc.interaction_set(pairs, mode='all', profile='cg').mask(trajectory)
    assert int(all_mask.sum()) <= int(any_mask.sum())


def test_a_set_of_one_matches_the_bare_selector(trajectory):
    single = sc.ChunkContact(20, 40, profile='cg').mask(trajectory)
    as_set = sc.interaction_set([(20, 40)], mode='all', profile='cg').mask(trajectory)
    assert (single == as_set).all()


def test_set_label_records_mode_and_pairs():
    selector = sc.interaction_set([(20, 40), (30, 50)], mode='all', profile='cg')
    assert selector.label() == 'all_20_40__30_50'


def test_set_rejects_an_unknown_mode():
    with pytest.raises(ValueError):
        sc.interaction_set([(20, 40)], mode='either', profile='cg')


def test_set_rejects_being_empty():
    with pytest.raises(ValueError):
        sc.InteractionSet([], mode='any')


# --------------------------------------------------------------- ensemble properties

def test_end_to_end_distance_works_without_alpha_carbons(trajectory):
    """select_atom_indices('alpha') returns nothing for a bead model.

    The chain ends have to be the first and last bead instead, or this raises IndexError.
    """
    assert profiles.CG.has_backbone is False
    value = sc.calc_Re(trajectory, trajectory, profile='cg')
    assert np.isfinite(value) and value > 0


def test_backbone_observables_are_not_invented_for_a_bead_model(trajectory, sequence):
    """DSSP needs a backbone and the hydrogen bond search needs hydrogens.

    Neither exists here, so these are reported as not measurable rather than as zero.
    """
    values = sc.calculate_ensemble(trajectory, trajectory, len(sequence), profile='cg')
    assert np.isfinite(values['Rg'])
    assert np.isfinite(values['Re'])
    assert values['Frames'] == N_FRAMES
    for name in ('Heli', 'Beta', 'HB'):
        assert np.isnan(values[name]), name


def test_the_split_compacts_the_chain(trajectory, sequence):
    """Forming the contacts should not leave the chain the same size.

    Not a law, but on this fixture both conditions compact it, and a split that changed
    nothing at all would mean the mask was not being applied.
    """
    mask = sc.interaction_set([(20, 40), (30, 50)], mode='any', profile='cg').mask(trajectory)
    result = sc.compare_ensemble(trajectory, trajectory, len(sequence), mask, profile='cg')
    assert result['interaction']['Rg'] < result['non_interaction']['Rg']
    assert result['interaction']['Frames'] == 78
    assert result['non_interaction']['Frames'] == N_FRAMES - 78


def test_split_maps_honour_the_output_directory(trajectory, tmp_path):
    """They used to be written into the working directory, which is the data folder."""
    mask = sc.interaction_set([(20, 40)], profile='cg').mask(trajectory)
    sc.new_interaction_map('demo', trajectory, mask, profile='cg', output_dir=str(tmp_path))
    written = sorted(p.name for p in tmp_path.glob('*.csv'))
    assert written == ['demowith_1.2_contact_df_1201.csv',
                       'demowithout_1.2_contact_df_1201.csv']


def test_fixture_folder_is_not_written_to():
    """The tutorial runs against tests/fixtures; nothing may be created there."""
    allowed = {NAME + '.pdb', NAME + '.dcd', NAME + '.fasta'}
    assert set(os.listdir(FIXTURES)) == allowed
