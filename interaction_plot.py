# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 15:00:05 2021

@author: ShaharGroup-fyu
"""

import matplotlib.pyplot as plt
import numpy as np

import networkx as nx

# Create the networkx graph
#### problem need to be solved: width of the graph
##Check this with the latest graphx library
# Create a  networkx graph object
def create_network(seq, length, size=10):
    '''

    :param seq:
    :param length:
    :param size:
    :return:
    '''
    graphg = nx.MultiDiGraph()
    # Add residue node to network graphics
    i = 1
    while i < (length + 1):
        graphg.add_node(i, residue=seq[i - 1], pos=(i, size))
        i += 1
    return graphg


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
def create_position(graphg):
    '''

    :param graphg:
    :return:
    '''
    # Get the position and the layout for the networkx nodes
    pos = nx.get_node_attributes(graphg, 'pos')
    layout = dict((n, graphg.node[n]["pos"]) for n in graphg.nodes())
    return pos, layout


# Plot function of sequence color coding
def color_text(pos, index, colorselec, seq, ax):
    '''

    :param pos:
    :param index:
    :param colorselec:
    :param seq:
    :param ax:
    :return:
    '''
    # Get the position of each residue
    (x, y) = pos[index + 1]
    # Get the name of each residue
    label = seq[index]  # this makes "1" and 1 labeled the same
    # Color coding the sequence
    ax.text(
        x - 0.2,
        y,
        label,
        color=colorselec,
        transform=ax.transData,
        clip_on=True, fontfamily='monospace'
    )


# Color coding the residue based on the position matrix and sequence color code.
def create_color_coding(seq, graphg, pos, negacharged, posicharged, aromatic, figuresize=(10, 10), nodesize=0.1):
    '''

    :param seq:
    :param graphg:
    :param pos:
    :param negacharged:
    :param posicharged:
    :param aromatic:
    :param figuresize:
    :param nodesize:
    :return:
    '''
    fig = plt.figure(figsize=figuresize)
    ax = fig.add_subplot(111)
    seqdict = {}
    for index, i in enumerate(seq):
        seqdict[index + 1] = i
    nx.draw(graphg, pos, labels=seqdict, with_labels=False, node_size=nodesize, ax=ax)
    # Set residue color based on their types
    for index, node in enumerate(graphg):
        if index in negacharged:
            color_text(pos, index, 'r', seq, ax)
        elif index in posicharged:
            color_text(pos, index, 'b', seq, ax)
        elif index in aromatic:
            color_text(pos, index, 'orange', seq, ax)
        else:
            color_text(pos, index, 'black', seq, ax)
    return fig, ax


# Plot selected interaction line
def interaction_plotting(interaction, layout, ax, inter_type):
    '''

    :param interaction:
    :param layout:
    :param ax:
    :param intertype:
    :param plot_c:
    :return:
    '''
    if inter_type == 2:
        colorset = 'green'
        connect = "arc3,rad=-0.5"
    if inter_type == 1:
        colorset = 'lightgreen'
        connect = "arc3,rad=-0.5"
    if inter_type == -2:
        colorset = 'red'
        connect = "arc3,rad=0.5"
    if inter_type == -1:
        colorset = 'orange'
        connect = "arc3,rad=0.5"
    # During the debug process, we only consider the condition when plot_c>0.
    if inter_type > 0:
        #full_strength, overall_strength = interaction_strength.calculate_overall_strength(interaction, inter_type)
        full_strength = []
        overall_strength = 1
    else:
        full_strength = []
        overall_strength = -1
    for index, data in interaction[interaction['plot_value']==inter_type].iterrows():
        strength=data['relative_strength']
        distance = data['distance']
        r1 = data['r_1']
        r2 = data['r_2']
        if distance <=2:
            continue
        if inter_type > 0:
            #    linewidth = raw_value_new[index]
            linewidth = 2 * data['relative_strength']
            distance = data['distance']
        else:
            #    linewidth = -1 * raw_value_new[index]
            linewidth = -2 * data['relative_strength'] * inter_type * 0.1
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

def interaction_map(seq, length, interaction, figname):
    # Create network object
    graphg = create_network(seq, length)
    # Obtain position list
    (pos, layout) = create_position(graphg)
    # Color coding different residues
    (negacharged, posicharged, aromatic) = seq_color(seq)
    (fig, ax) = create_color_coding(seq, graphg, pos, negacharged, posicharged, aromatic)
    # Plot interaction between each residue
    att1 = interaction_plotting(interaction, layout, ax,  1)
    att2 = interaction_plotting(interaction, layout, ax,  2)
    graphg.add_edge(1, 5)
    # Save the plot to png file
    plt.savefig(figname + '.png')
    plt.savefig(figname + '.svg')
    # Show the plot
    plt.show()
    # Return the interaction strength


