"""Tests for plotting module."""

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pytest

from idp_interaction_map.plotting import (
    create_network,
    identify_residue_types,
    get_residue_color,
    plot_interactions,
)

# Use non-interactive backend for testing
matplotlib.use("Agg")


def test_create_network(sample_sequence):
    """Test network graph creation."""
    length = len(sample_sequence)
    graph = create_network(sample_sequence, length)

    # Check graph properties
    assert isinstance(graph, nx.MultiDiGraph)
    assert len(graph.nodes()) == length

    # Check node attributes
    for i, residue in enumerate(sample_sequence, start=1):
        assert graph.nodes[i]["residue"] == residue
        assert "pos" in graph.nodes[i]


def test_identify_residue_types():
    """Test residue type categorization."""
    seq = "DERKFYW"
    negative, positive, aromatic = identify_residue_types(seq)

    assert negative == [0, 1]  # D, E
    assert positive == [2, 3]  # R, K
    assert aromatic == [4, 5, 6]  # F, Y, W


def test_get_residue_color():
    """Test color assignment for residues."""
    assert get_residue_color("D") == "red"  # negative
    assert get_residue_color("E") == "red"  # negative
    assert get_residue_color("R") == "blue"  # positive
    assert get_residue_color("K") == "blue"  # positive
    assert get_residue_color("F") == "yellow"  # aromatic
    assert get_residue_color("A") == "orange"  # hydrophobic
    assert get_residue_color("S") == "black"  # polar


def test_get_residue_color_invalid():
    """Test error handling for invalid residue."""
    with pytest.raises(ValueError, match="Unrecognized residue"):
        get_residue_color("X")


def test_plot_interactions():
    """Test interaction plotting."""
    # Create simple interaction data
    interaction_df = pd.DataFrame(
        {
            "r_1": [1, 2],
            "r_2": [5, 6],
            "cont_prob": [0.8, 0.7],
            "relative_strength": [1.5, 1.2],
            "distance": [4, 4],
            "plot_value": [2, 1],
        }
    )

    # Create simple layout
    layout = {i: (i, 10) for i in range(1, 7)}

    # Create figure
    fig, ax = plt.subplots()

    # Test favorable interactions
    strength = plot_interactions(interaction_df, layout, ax, interaction_type=2)
    assert strength == 1

    plt.close(fig)


def test_plot_interactions_unfavorable():
    """Test plotting unfavorable interactions."""
    interaction_df = pd.DataFrame(
        {
            "r_1": [1],
            "r_2": [5],
            "cont_prob": [0.3],
            "relative_strength": [-1.5],
            "distance": [4],
            "plot_value": [-2],
        }
    )

    layout = {i: (i, 10) for i in range(1, 6)}
    fig, ax = plt.subplots()

    strength = plot_interactions(interaction_df, layout, ax, interaction_type=-2)
    assert strength == -1

    plt.close(fig)


def test_plot_interactions_filters_short_range():
    """Test that short-range interactions are filtered out."""
    interaction_df = pd.DataFrame(
        {
            "r_1": [1, 2],
            "r_2": [3, 10],  # distance 2 and 8
            "cont_prob": [0.8, 0.7],
            "relative_strength": [1.5, 1.2],
            "distance": [2, 8],
            "plot_value": [2, 2],
        }
    )

    layout = {i: (i, 10) for i in range(1, 11)}
    fig, ax = plt.subplots()

    # Should only plot the second interaction (distance > 4)
    plot_interactions(interaction_df, layout, ax, interaction_type=2)

    plt.close(fig)
