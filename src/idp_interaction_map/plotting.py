"""Visualization module for protein interaction networks."""

import logging
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)

# Color scheme for amino acid types
AA_TYPE_MAP = {
    "Y": "aromatic",
    "F": "aromatic",
    "W": "aromatic",
    "R": "positive",
    "H": "positive",
    "K": "positive",
    "D": "negative",
    "E": "negative",
    "S": "polar",
    "T": "polar",
    "Q": "polar",
    "N": "polar",
    "A": "hydrophobic",
    "V": "hydrophobic",
    "I": "hydrophobic",
    "L": "hydrophobic",
    "M": "hydrophobic",
    "C": "hydrophobic",
    "G": "hydrophobic",
    "P": "hydrophobic",
}

COLOR_MAP = {
    "aromatic": "yellow",
    "positive": "blue",
    "negative": "red",
    "polar": "black",
    "hydrophobic": "orange",
}


def create_network(seq: str, length: int, size: int = 10) -> nx.MultiDiGraph:
    """
    Create NetworkX graph representing protein sequence.

    Args:
        seq: Protein sequence string
        length: Sequence length
        size: Y-coordinate for node placement

    Returns:
        NetworkX MultiDiGraph with residue nodes
    """
    graph = nx.MultiDiGraph()

    for i in range(1, length + 1):
        graph.add_node(i, residue=seq[i - 1], pos=(i, size))

    return graph


def identify_residue_types(seq: str) -> Tuple[List[int], List[int], List[int]]:
    """
    Categorize residues by charge and aromaticity.

    Args:
        seq: Protein sequence string

    Returns:
        Tuple of (negative_charged, positive_charged, aromatic) index lists
    """
    negative = []
    positive = []
    aromatic = []

    for index, residue in enumerate(seq):
        if residue in ["D", "E"]:
            negative.append(index)
        elif residue in ["R", "K", "H"]:
            positive.append(index)
        elif residue in ["F", "Y", "W"]:
            aromatic.append(index)

    return negative, positive, aromatic


def get_residue_color(residue: str) -> str:
    """
    Get color for amino acid residue.

    Args:
        residue: Single letter amino acid code

    Returns:
        Color string

    Raises:
        ValueError: If residue type is not recognized
    """
    aa_type = AA_TYPE_MAP.get(residue)
    if aa_type is None:
        raise ValueError(f"Unrecognized residue: {residue}")
    return COLOR_MAP[aa_type]


def plot_residue(
    pos: Dict, index: int, seq: str, ax: plt.Axes
) -> None:
    """
    Plot a single residue with appropriate color.

    Args:
        pos: Position dictionary from NetworkX
        index: Residue index (0-indexed)
        seq: Protein sequence
        ax: Matplotlib axes object
    """
    x, y = pos[index + 1]
    residue = seq[index]
    color = get_residue_color(residue)

    ax.plot(
        x - 0.2,
        y,
        marker="o",
        color=color,
        ms=7,
        markeredgecolor="black",
    )


def create_sequence_visualization(
    seq: str,
    graph: nx.MultiDiGraph,
    pos: Dict,
    negative: List[int],
    positive: List[int],
    aromatic: List[int],
    figsize: Tuple[int, int] = (10, 10),
    nodesize: float = 0.1,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create color-coded sequence visualization.

    Args:
        seq: Protein sequence
        graph: NetworkX graph
        pos: Position dictionary
        negative: Indices of negatively charged residues
        positive: Indices of positively charged residues
        aromatic: Indices of aromatic residues
        figsize: Figure size tuple
        nodesize: Node size for NetworkX drawing

    Returns:
        Tuple of (figure, axes)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Create residue labels
    seq_dict = {i + 1: residue for i, residue in enumerate(seq)}

    # Draw base network
    nx.draw(graph, pos, labels=seq_dict, with_labels=False, node_size=nodesize, ax=ax)

    # Plot each residue with appropriate color
    for index in range(len(seq)):
        plot_residue(pos, index, seq, ax)

    return fig, ax


def plot_interactions(
    interaction_df: pd.DataFrame,
    layout: Dict,
    ax: plt.Axes,
    interaction_type: int,
) -> float:
    """
    Plot interaction lines for a specific interaction type.

    Args:
        interaction_df: DataFrame with interaction data
        layout: Node position layout
        ax: Matplotlib axes
        interaction_type: Type of interaction to plot (2, 1, -1, -2)

    Returns:
        Overall interaction strength (1 for favorable, -1 for unfavorable)
    """
    # Configure visualization based on interaction type
    if interaction_type == 2:
        color = "green"
        connection_style = "arc3,rad=-0.5"
        strength_sign = 1
    elif interaction_type == 1:
        color = "lightgreen"
        connection_style = "arc3,rad=-0.5"
        strength_sign = 1
    elif interaction_type == -1:
        color = "orange"
        connection_style = "arc3,rad=0.5"
        strength_sign = -1
    elif interaction_type == -2:
        color = "red"
        connection_style = "arc3,rad=0.5"
        strength_sign = -1
    else:
        raise ValueError(f"Invalid interaction type: {interaction_type}")

    # Filter interactions by type and minimum distance
    selected = interaction_df[
        (interaction_df["plot_value"] == interaction_type) & (interaction_df["distance"] > 4)
    ]

    logger.info(f"Plotting {len(selected)} interactions of type {interaction_type}")

    # Plot each interaction
    for _, data in selected.iterrows():
        strength = data["cont_prob"]
        relative_strength = data["relative_strength"]
        r1, r2 = int(data["r_1"]), int(data["r_2"])

        # Calculate line width based on interaction strength
        if interaction_type > 0:
            linewidth = 4 * strength * (1 + relative_strength)
        else:
            linewidth = -4 * strength * (relative_strength - 1)

        # Get positions with offset for visualization
        x1, y1 = layout[r1][0] - 0.2, layout[r1][1]
        x2, y2 = layout[r2][0] + 0.2, layout[r2][1]

        # Draw interaction arc
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x2, y2),
            arrowprops=dict(
                arrowstyle="-",
                color=color,
                shrinkA=10,
                shrinkB=10,
                lw=linewidth,
                patchA=None,
                patchB=None,
                connectionstyle=connection_style,
            ),
        )

    return strength_sign


def create_interaction_map(
    seq: str,
    length: int,
    interaction_df: pd.DataFrame,
    output_name: str,
) -> None:
    """
    Generate and save complete interaction map visualization.

    Args:
        seq: Protein sequence
        length: Sequence length
        interaction_df: DataFrame with normalized interaction data
        output_name: Base name for output files
    """
    logger.info(f"Creating interaction map for {output_name}")

    # Create network graph
    graph = create_network(seq, length)
    pos = nx.get_node_attributes(graph, "pos")
    layout = dict((n, graph.nodes[n]["pos"]) for n in graph.nodes())

    # Identify residue types
    negative, positive, aromatic = identify_residue_types(seq)

    # Create base visualization
    fig, ax = create_sequence_visualization(
        seq, graph, pos, negative, positive, aromatic
    )

    # Plot all interaction types
    plot_interactions(interaction_df, layout, ax, 2)  # Strong favorable
    plot_interactions(interaction_df, layout, ax, 1)  # Weak favorable
    plot_interactions(interaction_df, layout, ax, -1)  # Weak unfavorable
    plot_interactions(interaction_df, layout, ax, -2)  # Strong unfavorable

    # Save output files
    plt.savefig(f"{output_name}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_name}.svg", bbox_inches="tight")

    logger.info(f"Saved interaction map to {output_name}.png and {output_name}.svg")

    plt.show()
