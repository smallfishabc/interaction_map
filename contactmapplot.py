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
def create_network(seq, length, size=10):
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
    negacharged = []
    posicharged = []
    aromatic = []
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
    pos = nx.get_node_attributes(graphg, 'pos')
    layout = dict((n, graphg.node[n]["pos"]) for n in graphg.nodes())
    return pos, layout


# Plot function of sequence color coding
def color_text(pos, index, colorselec, seq, ax):
    (x, y) = pos[index + 1]
    label = seq[index]  # this makes "1" and 1 labeled the same
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
def interaction_plotting(reaction, raw_value, pairs, layout, ax, intertype, strength, targetmap):
    inot = strength  # default
    if intertype == 'att':
        if inot == 2:
            colorset = 'green'
            connect = "arc3,rad=-0.5"
        if inot == 1:
            colorset = 'aquamarine'
            connect = "arc3,rad=-0.5"
    elif intertype == 'rep':
        if inot == -2:
            colorset = 'red'
            connect = "arc3,rad=0.5"
        if inot == -1:
            colorset = 'mistyrose'
            connect = "arc3,rad=0.5"
    # Remove unnecessary pairs
    indexlist = []
    for index, i in enumerate(reaction):
        if i != inot:
            indexlist.append(index)
    pairsnew = np.delete(pairs, indexlist, axis=0)
    reaction_value = np.delete(reaction, indexlist, axis=0)
    raw_value_new = np.delete(raw_value, indexlist, axis=0)
    target_new = np.delete(targetmap, indexlist, axis=0)
    np.savetxt(intertype + str(inot) + 'interaction_pair.csv', pairsnew, delimiter=',')
    for index, edge in enumerate(pairsnew):
        if inot > 0:
        #    linewidth = raw_value_new[index]
            linewidth = 2*target_new[index]*inot
            print(linewidth)
        else:
        #    linewidth = -1 * raw_value_new[index]
            linewidth = -10*target_new[index]*inot
        # Adjust location to improve visualization effect
        a = layout[edge[0]][0] - 0.2
        b = layout[edge[0]][1]
        c = layout[edge[1]][0] + 0.2
        d = layout[edge[1]][1]
        # Plot Interaction between pairs
        ax.annotate("",
                    xy=(a, b),
                    xytext=(c, d),
                    arrowprops=dict(arrowstyle="-", color=colorset,
                                    shrinkA=10, shrinkB=10, lw=linewidth,
                                    patchA=None, patchB=None,
                                    connectionstyle=connect,
                                    ), )


def interaction_map(seq, length, interaction, raw_value, pairs, figname,targetmap):
    # Create network object
    graphg = create_network(seq, length)
    # Obtain position list
    (pos, layout) = create_position(graphg)
    # Color coding different residues
    (negacharged, posicharged, aromatic) = seq_color(seq)
    (fig, ax) = create_color_coding(seq, graphg, pos, negacharged, posicharged, aromatic)
    # Plot interaction between each residue
    interaction_plotting(interaction, raw_value, pairs, layout, ax, 'att', 1,targetmap)
    interaction_plotting(interaction, raw_value, pairs, layout, ax, 'att', 2,targetmap)
    interaction_plotting(interaction, raw_value, pairs, layout, ax, 'rep', -1,targetmap)
    interaction_plotting(interaction, raw_value, pairs, layout, ax, 'rep', -2,targetmap)
    # Save and demonstrate the plot
    plt.savefig(figname + '.png')
    plt.show()
