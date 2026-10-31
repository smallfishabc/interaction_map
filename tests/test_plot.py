"""Plotting invariants.

These pin the geometry and styling that the figure depends on, independently of how the
sequence track is drawn. They are what makes it safe to replace networkx: if the arcs are
still anchored to the same coordinates with the same widths and colours, and the rendered
image is unchanged, the substitution is sound.
"""
import os

import numpy as np
import pandas as pd
import pytest

from ensemble_map import interaction_plot as ip
from ensemble_map import normalization as nl
from ensemble_map import profiles

from conftest import requires_p53

SEQ = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPD'


class RecordingAxis(object):
    """Captures annotate() calls instead of drawing them."""

    def __init__(self):
        self.calls = []

    def annotate(self, text, xy=None, xytext=None, arrowprops=None, **kwargs):
        self.calls.append({'xy': xy, 'xytext': xytext,
                           'lw': arrowprops['lw'], 'color': arrowprops['color'],
                           'connectionstyle': arrowprops['connectionstyle']})


def normalized(pairs, profile='all_atom'):
    """pairs: list of (r_1, r_2, cont_prob)."""
    table = pd.DataFrame(pairs, columns=['r_1', 'r_2', 'cont_prob'])
    table['distance'] = table['r_2'] - table['r_1']
    return nl.normalization(table, profile=profile)


# --------------------------------------------------------------- sequence track

def test_network_keys_are_one_based_sequence_positions():
    layout = ip.create_network(SEQ, len(SEQ))
    assert sorted(layout) == list(range(1, len(SEQ) + 1))


def test_network_preserves_sequence_order():
    """Keys are iterated in sequence order, which create_color_coding relies on."""
    layout = ip.create_network(SEQ, len(SEQ))
    assert list(layout) == list(range(1, len(SEQ) + 1))


def test_layout_places_residue_n_at_x_equals_n():
    layout = ip.create_network(SEQ, len(SEQ), size=10)
    pos, same = ip.create_position(layout)
    for node in layout:
        assert layout[node] == (node, 10)
    assert pos == same == layout


def test_layout_size_argument_sets_the_baseline_height():
    graph = ip.create_network(SEQ, len(SEQ), size=4)
    _pos, layout = ip.create_position(graph)
    assert layout[1] == (1, 4)


def test_seq_color_classifies_by_zero_based_index():
    neg, pos, arom = ip.seq_color('DERKHFYWA')
    assert neg == [0, 1]
    assert pos == [2, 3, 4]
    assert arom == [5, 6, 7]


# --------------------------------------------------------------- arc selection

def test_arc_endpoints_are_offset_from_the_named_residues():
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    table = normalized([(5, 25, 0.9)])
    level = int(table['plot_value'].iloc[0])
    axis = RecordingAxis()
    ip.interaction_plotting(table, layout, axis, level, profile='all_atom')
    assert len(axis.calls) == 1
    call = axis.calls[0]
    assert call['xy'] == (layout[5][0] - 0.2, layout[5][1])
    assert call['xytext'] == (layout[25][0] + 0.2, layout[25][1])


def test_first_residue_is_addressable():
    """r_1 == 1 must resolve; it used to be the risky end of the numbering."""
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    table = normalized([(1, 21, 0.9)])
    axis = RecordingAxis()
    ip.interaction_plotting(table, layout, axis, int(table['plot_value'].iloc[0]),
                            profile='all_atom')
    assert axis.calls[0]['xy'][0] == pytest.approx(0.8)


def test_last_residue_is_addressable():
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    last = len(SEQ)
    table = normalized([(last - 20, last, 0.9)])
    axis = RecordingAxis()
    ip.interaction_plotting(table, layout, axis, int(table['plot_value'].iloc[0]),
                            profile='all_atom')
    assert axis.calls[0]['xytext'][0] == pytest.approx(last + 0.2)


def test_only_the_requested_level_is_drawn():
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    table = normalized([(1, 21, 0.9), (2, 22, 0.9), (3, 23, 1e-6)])
    levels = table['plot_value'].tolist()
    for level in (2, 1, -1, -2):
        axis = RecordingAxis()
        ip.interaction_plotting(table, layout, axis, level, profile='all_atom')
        assert len(axis.calls) == levels.count(level)


@pytest.mark.parametrize('profile_name,min_sep', [('all_atom', 3), ('cg', 5)])
def test_min_sep_filters_short_separations(profile_name, min_sep):
    profile = profiles.get_profile(profile_name)
    assert profile.min_sep == min_sep
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    # One pair exactly at the threshold, one just below it.
    table = normalized([(1, 1 + min_sep, 0.9), (10, 10 + min_sep - 1, 0.9)],
                       profile=profile_name)
    drawn = []
    for level in (2, 1, -1, -2):
        axis = RecordingAxis()
        ip.interaction_plotting(table, layout, axis, level, profile=profile)
        drawn.extend(axis.calls)
    separations = [round(c['xytext'][0] - c['xy'][0] - 0.4) for c in drawn]
    assert all(s >= min_sep for s in separations)


# --------------------------------------------------------------- styling

@pytest.mark.parametrize('level,color,above', [
    (2, 'green', True),
    (1, 'lightgreen', True),
    (-1, 'orange', False),
    (-2, 'red', False),
])
def test_colour_and_curvature_per_level(level, color, above):
    assert ip.INTERACTION_STYLE[level][0] == color
    curvature = ip.INTERACTION_STYLE[level][1]
    assert curvature == ('arc3,rad=-0.5' if above else 'arc3,rad=0.5')


def test_unknown_level_is_rejected():
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    with pytest.raises(ValueError):
        ip.interaction_plotting(normalized([(1, 21, 0.9)]), layout, RecordingAxis(), 0)


@pytest.mark.parametrize('profile_name,scale', [('all_atom', 2.0), ('cg', 4.0)])
def test_line_width_formula(profile_name, scale):
    profile = profiles.get_profile(profile_name)
    assert profile.line_width_scale == scale
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    table = normalized([(1, 21, 0.9)], profile=profile_name)
    level = int(table['plot_value'].iloc[0])
    axis = RecordingAxis()
    ip.interaction_plotting(table, layout, axis, level, profile=profile)
    p = table['cont_prob'].iloc[0]
    rel = table['relative_strength'].iloc[0]
    expected = scale * p * (1 + rel) if level > 0 else -scale * p * (rel - 1)
    assert axis.calls[0]['lw'] == pytest.approx(expected, rel=1e-15)


def test_line_width_is_positive_for_both_signs():
    """The two branches of the width formula must both give a drawable width."""
    graph = ip.create_network(SEQ, len(SEQ))
    _pos, layout = ip.create_position(graph)
    table = normalized([(1, 21, 0.9), (2, 40, 1e-5)])
    for level in (2, 1, -1, -2):
        axis = RecordingAxis()
        ip.interaction_plotting(table, layout, axis, level, profile='all_atom')
        for call in axis.calls:
            assert call['lw'] > 0


# --------------------------------------------------------------- residue styles

def test_marker_style_colours_by_residue_class():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ip.color_text({1: (1, 10)}, 0, 'ignored', 'W', ax, residue_style='marker')
    assert ax.lines[0].get_color() == ip.RESIDUE_COLOR['aromatic']
    plt.close(fig)


def test_letter_style_draws_the_one_letter_code():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ip.color_text({1: (1, 10)}, 0, 'r', 'D', ax, residue_style='letter')
    assert [t.get_text() for t in ax.texts] == ['D']
    plt.close(fig)


def test_every_standard_residue_has_a_colour():
    for residue in 'ACDEFGHIKLMNPQRSTVWY':
        assert residue in ip.RESIDUE_TYPE
        assert ip.RESIDUE_TYPE[residue] in ip.RESIDUE_COLOR


def test_unknown_residue_is_reported_clearly():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    with pytest.raises(ValueError) as excinfo:
        ip.color_text({1: (1, 10)}, 0, 'black', 'X', ax, residue_style='marker')
    assert 'X' in str(excinfo.value)
    plt.close(fig)


def test_unknown_style_is_rejected():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        ip.color_text({1: (1, 10)}, 0, 'black', 'A', ax, residue_style='squiggle')
    plt.close(fig)


# --------------------------------------------------------------- whole figure

def test_interaction_map_returns_figure_and_axis(tmp_path):
    import matplotlib.pyplot as plt
    table = normalized([(1, 21, 0.9), (5, 40, 0.5), (3, 30, 1e-6)])
    os.chdir(str(tmp_path))
    fig, ax = ip.interaction_map(SEQ, len(SEQ), table, 'm', show=False, save=False)
    assert fig is not None and ax is not None
    plt.close(fig)


def test_interaction_map_writes_png_and_svg(tmp_path):
    import matplotlib.pyplot as plt
    table = normalized([(1, 21, 0.9)])
    os.chdir(str(tmp_path))
    fig, _ax = ip.interaction_map(SEQ, len(SEQ), table, 'figure', show=False, save=True)
    assert (tmp_path / 'figure.png').exists()
    assert (tmp_path / 'figure.svg').exists()
    plt.close(fig)


def test_interaction_map_can_skip_saving(tmp_path):
    import matplotlib.pyplot as plt
    table = normalized([(1, 21, 0.9)])
    os.chdir(str(tmp_path))
    fig, _ax = ip.interaction_map(SEQ, len(SEQ), table, 'nothing', show=False, save=False)
    assert not (tmp_path / 'nothing.png').exists()
    plt.close(fig)


@requires_p53
def test_rendered_figure_matches_the_published_png(tmp_path, p53_sequence, p53_reference, p53_png):
    """The strongest guard available: any change to the drawing path must keep this.

    Uses the letter renderer, which is what produced the published figures.
    """
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    letter = profiles.Profile(**dict(profiles.ALL_ATOM.__dict__, residue_style='letter'))
    os.chdir(str(tmp_path))
    fig, _ax = ip.interaction_map(p53_sequence, len(p53_sequence), p53_reference.copy(),
                                  'p53', profile=letter, show=False, save=True)
    plt.close(fig)

    reference = mpimg.imread(p53_png)
    produced = mpimg.imread(str(tmp_path / 'p53.png'))
    assert produced.shape == reference.shape
    differing = int((np.abs(produced.astype(float) - reference.astype(float)).max(axis=2) > 0).sum())
    assert differing == 0, '{0} pixels differ from the published figure'.format(differing)
