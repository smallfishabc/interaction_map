# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 15:00:05 2021

@author: ShaharGroup-fyu
"""

import matplotlib.pyplot as plt
import numpy as np

import networkx as nx


def seqcolor(seq):
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


def createnetwork(seq, length, size=10):
    graphg = nx.MultiDiGraph()
    # Add residue node to network graphics
    i = 1
    while i < (length + 1):
        graphg.add_node(i, residue=seq[i - 1], pos=(i, size))
        i += 1
    return graphg


# Produce position matrix and layout
def createposition(graphg):
    pos = nx.get_node_attributes(graphg, 'pos')
    layout = dict((n, graphg.node[n]["pos"]) for n in graphg.nodes())
    return pos, layout


def colortext(pos, index, colorselec, seq, ax):
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


def createcolorcoding(seq, graphg, pos, negacharged, posicharged, aromatic, figuresize=(10, 10), nodesize=0.1):
    fig = plt.figure(figsize=figuresize)
    ax = fig.add_subplot(111)
    seqdict = {}
    for index, i in enumerate(seq):
        seqdict[index + 1] = i
    nx.draw(graphg, pos, labels=seqdict, with_labels=False, node_size=nodesize, ax=ax)
    # Set residue color based on their types
    for index, node in enumerate(graphg):
        if index in negacharged:
            colortext(pos, index, 'r', seq, ax)
        elif index in posicharged:
            colortext(pos, index, 'b', seq, ax)
        elif index in aromatic:
            colortext(pos, index, 'orange', seq, ax)
        else:
            colortext(pos, index, 'black', seq, ax)
    return fig, ax


# Remove attraction interaction to plot repulsive interaction
def interactionploting(reaction, pairs, layout, ax, intertype):
    inot = 1     # default
    if intertype == 'att':
        inot = 1
        colorset = 'aquamarine'
        connect = "arc3,rad=-0.5"
    elif intertype == 'rep':
        inot = -1
        colorset = 'mistyrose'
        connect = "arc3,rad=0.5"
    indexlist = []
    for index, i in enumerate(reaction):
        if i != inot:
            indexlist.append(index)
    pairsnew = np.delete(pairs, indexlist, axis=0)
    for index, edge in enumerate(pairsnew):
        #         Adjust location to improve visualization effect
        a = layout[edge[0]][0] - 0.2
        b = layout[edge[0]][1]
        c = layout[edge[1]][0] + 0.2
        d = layout[edge[1]][1]
        ax.annotate("",
                    xy=(a, b),
                    xytext=(c, d),
                    arrowprops=dict(arrowstyle="-", color=colorset,
                                    shrinkA=10, shrinkB=10, lw=1,
                                    patchA=None, patchB=None,
                                    connectionstyle=connect,
                                    ), )


def interactionmap(seq, length, interaction, pairs):
    graphg = createnetwork(seq, length)
    (pos, layout) = createposition(graphg)
    (negacharged, posicharged, aromatic) = seqcolor(seq)
    (fig, ax) = createcolorcoding(seq, graphg, pos, negacharged, posicharged, aromatic)
    interactionploting(interaction, pairs, layout, ax, 'att')
    interactionploting(interaction, pairs, layout, ax, 'rep')
    plt.show()
