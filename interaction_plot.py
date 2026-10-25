# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 15:00:05 2021

@author: ShaharGroup-fyu
"""

import matplotlib.pyplot as plt
import numpy as np

import profiles

# Colour networkx used for its nodes, kept so the sequence track renders as before.
NODE_COLOR = '#1f78b4'

# Residue classes used by the marker style, following the coarse-grained branch.
RESIDUE_TYPE = {'Y': 'aromatic', 'F': 'aromatic', 'W': 'aromatic',
                'R': 'positive', 'H': 'positive', 'K': 'positive',
                'D': 'negative', 'E': 'negative',
                'S': 'polar', 'T': 'polar', 'Q': 'polar', 'N': 'polar',
                'A': 'hydrophobic', 'V': 'hydrophobic', 'I': 'hydrophobic',
                'L': 'hydrophobic', 'M': 'hydrophobic', 'C': 'hydrophobic',
                'G': 'hydrophobic', 'P': 'hydrophobic'}

RESIDUE_COLOR = {'aromatic': 'yellow', 'positive': 'blue', 'negative': 'red',
                 'polar': 'black', 'hydrophobic': 'orange'}

# Arc colour and curvature per interaction level. Positive levels arc above the sequence,
# negative levels below it.
INTERACTION_STYLE = {2: ('green', 'arc3,rad=-0.5'),
                     1: ('lightgreen', 'arc3,rad=-0.5'),
                     -1: ('orange', 'arc3,rad=0.5'),
                     -2: ('red', 'arc3,rad=0.5')}


# Positions of the residues along the sequence track.
# This used to build a networkx MultiDiGraph whose nodes carried a 'pos' attribute. The
# nodes sat on a straight horizontal line and no edge was ever drawn between them, so the
# graph was only ever used as a dictionary of coordinates. It is now exactly that.
def create_network(seq, length, size=10):
    '''
    Residue node positions, keyed by 1-based sequence position.

    :param seq: protein sequence
    :param length: sequence length
    :param size: height of the sequence track
    :return: dict mapping residue number to (x, y)
    '''
    layout = {}
    # Add residue node to network graphics
    i = 1
    while i < (length + 1):
        layout[i] = (i, size)
        i += 1
    return layout


# To help us identify different amino acid. We color code negative charged residues , positive
# charged residues and aromatic residues with different colors.
def seq_color(seq):
    '''

    :param seq:
    :return:
    '''
    # Create empty list
    negacharged = []
    posicharged = []
    aromatic = []
    # Identify special residues and append into corresponding list
    for index, i in enumerate(seq):
        if i in ['D', 'E']:
            negacharged.append(index)
        elif i in ['R', 'K', 'H']:
            posicharged.append(index)
        elif i in ['F', 'Y', 'W']:
            aromatic.append(index)
    return negacharged, posicharged, aromatic


# Produce position matrix and layout.
def create_position(layout):
    '''
    Kept because callers on the other branches still unpack two values. The two dicts were
    always identical: one came from nx.get_node_attributes, the other from reading the same
    attribute back off each node.

    :param layout: the dict returned by create_network
    :return: the same dict twice
    '''
    return layout, layout


# Plot function of sequence color coding
def color_text(pos, index, colorselec, seq, ax, residue_style='marker'):
    '''

    :param pos:
    :param index: 0-based position along the sequence
    :param colorselec: colour used by the letter style
    :param seq:
    :param ax:
    :param residue_style: 'marker' for coloured dots, 'letter' for coloured one letter codes
    :return:
    '''
    # Get the position of each residue
    (x, y) = pos[index + 1]
    # Get the name of each residue
    label = seq[index]  # this makes "1" and 1 labeled the same
    if residue_style == 'letter':
        # Color coding the sequence
        ax.text(
            x - 0.2,
            y,
            label,
            color=colorselec,
            transform=ax.transData,
            clip_on=True, fontfamily='monospace'
        )
    elif residue_style == 'marker':
        if label not in RESIDUE_TYPE:
            raise ValueError('Unknown residue {0!r} at sequence position {1}'.format(label, index + 1))
        ax.plot(
            x - 0.2,
            y,
            marker='o',
            color=RESIDUE_COLOR[RESIDUE_TYPE[label]], ms=7, markeredgecolor='black'
        )
    else:
        raise ValueError('Unknown residue_style {0!r}'.format(residue_style))


# Color coding the residue based on the position matrix and sequence color code.
def create_color_coding(seq, layout, pos, negacharged, posicharged, aromatic, figuresize=(10, 10), nodesize=0.1,
                        residue_style='marker'):
    '''

    :param seq:
    :param layout: residue positions from create_network
    :param pos: the same dict, kept for call compatibility
    :param negacharged:
    :param posicharged:
    :param aromatic:
    :param figuresize:
    :param nodesize:
    :param residue_style: 'marker' or 'letter'
    :return:
    '''
    fig = plt.figure(figsize=figuresize)
    ax = fig.add_subplot(111)
    # This stood in for nx.draw(graphg, pos, with_labels=False, node_size=nodesize, ax=ax).
    # Labels were switched off and the graph had no edges, so all it did was scatter the
    # node markers, silence the ticks and hide the axes. Reproduced here so the figure is
    # unchanged, with networkx no longer needed.
    fig.set_facecolor('w')
    nodes = list(layout)
    xy = np.asarray([layout[n] for n in nodes])
    node_collection = ax.scatter(xy[:, 0], xy[:, 1], s=nodesize, c=NODE_COLOR, marker='o',
                                 cmap=None, vmin=None, vmax=None, alpha=None,
                                 linewidths=None, edgecolors=None, label=None)
    ax.tick_params(axis='both', which='both', bottom=False, left=False,
                   labelbottom=False, labelleft=False)
    node_collection.set_zorder(2)
    ax.set_axis_off()
    # Set residue color based on their types
    for index, node in enumerate(nodes):
        if index in negacharged:
            color_text(pos, index, 'r', seq, ax, residue_style)
        elif index in posicharged:
            color_text(pos, index, 'b', seq, ax, residue_style)
        elif index in aromatic:
            color_text(pos, index, 'orange', seq, ax, residue_style)
        else:
            color_text(pos, index, 'black', seq, ax, residue_style)
    return fig, ax


# Plot selected interaction line
def interaction_plotting(interaction, layout, ax, inter_type, profile=None):
    '''

    :param interaction:
    :param layout:
    :param ax:
    :param inter_type: interaction level to draw, one of 2, 1, -1, -2
    :param profile: analysis profile supplying min_sep and the line width scale
    :return:
    '''
    profile = profiles.get_profile(profile)
    if inter_type not in INTERACTION_STYLE:
        raise ValueError('No style defined for interaction level {0!r}'.format(inter_type))
    colorset, connect = INTERACTION_STYLE[inter_type]
    # During the debug process, we only consider the condition when plot_c>0.
    if inter_type > 0:
        overall_strength = 1
    else:
        overall_strength = -1
    # r_1 / r_2 already hold 1-based sequence positions, which is how the network nodes are
    # numbered, so they index the layout directly.
    select_df = interaction[(interaction['plot_value'] == inter_type) &
                            (interaction['distance'] >= profile.min_sep)]
    for index, data in select_df.iterrows():
        strength = data['cont_prob']
        rela_strength = data['relative_strength']
        r1 = int(data['r_1'])
        r2 = int(data['r_2'])
        if inter_type > 0:
            linewidth = profile.line_width_scale * strength * (1 + rela_strength)
        else:
            linewidth = -profile.line_width_scale * strength * (rela_strength - 1)
        # Adjust location to improve visualization effect
        a = layout[r1][0] - 0.2
        b = layout[r1][1]
        c = layout[r2][0] + 0.2
        d = layout[r2][1]
        # Plot Interaction between pairs
        ax.annotate("",
                    xy=(a, b),
                    xytext=(c, d),
                    arrowprops=dict(arrowstyle="-", color=colorset,
                                    shrinkA=10, shrinkB=10, lw=linewidth,
                                    patchA=None, patchB=None,
                                    connectionstyle=connect,
                                    ), )
    return overall_strength


def interaction_map(seq, length, interaction, figname, profile=None, show=True, save=True):
    """
    :param seq: protein sequence
    :param length: sequence length
    :param interaction: normalized interaction dataframe
    :param figname: output file name, without extension
    :param profile: analysis profile (all_atom / cg)
    :param show: call plt.show(). Turn this off for batch runs.
    :param save: write <figname>.png and <figname>.svg
    :return: matplotlib figure and axis
    """
    profile = profiles.get_profile(profile)
    # Create network object
    graphg = create_network(seq, length)
    # Obtain position list
    (pos, layout) = create_position(graphg)
    # Color coding different residues
    (negacharged, posicharged, aromatic) = seq_color(seq)
    (fig, ax) = create_color_coding(seq, graphg, pos, negacharged, posicharged, aromatic,
                                    residue_style=profile.residue_style)
    # Plot interaction between each residue
    att1 = interaction_plotting(interaction, layout, ax, 1, profile=profile)
    att2 = interaction_plotting(interaction, layout, ax, 2, profile=profile)
    rep1 = interaction_plotting(interaction, layout, ax, -1, profile=profile)
    rep2 = interaction_plotting(interaction, layout, ax, -2, profile=profile)
    # Save the plot to png file
    if save:
        plt.savefig(figname + '.png')
        plt.savefig(figname + '.svg')
    # Show the plot
    if show:
        plt.show()
    return fig, ax
