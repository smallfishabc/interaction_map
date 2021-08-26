# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os

import numpy as np

import contactmapgeneration as cg
import contactmapplot as cm
import normalization as nl
import path


def averagecontMean(trajlist):
    nn = 0
    repeats=len(trajlist)
    framesum = 0
    for n in trajlist:
        traj, frame = cg.loadtraj(n)
        [pairs, contMean] = cg.computecontact(traj, 6, 0.8)
        if nn == 0:
            contMeansum = np.empty_like(contMean)
        contMeansum += np.multiply(frame, contMean)
        framesum += frame
        del traj
        frame = 0
        nn += 1
    contMean = np.true_divide(contMeansum, framesum)
    return contMean, pairs


def interactionmapmain(proteinpath, psi, residue):
    print(proteinpath, psi, residue)
    os.chdir(proteinpath)
    [seq, length] = path.readsequence()
    path.subdir(proteinpath, psi, residue)
    trajlist = [4]
    contMean, pairs = averagecontMean(trajlist)
    interaction = nl.normalization(contMean, pairs)
    print(interaction)
    cm.interactionmap(seq, length, interaction, pairs)

#interactionmapmain('F:\globus\simulation_contactmap_validation\GS22-summary','S_0','BB')