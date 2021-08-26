# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:23:10 2021

@author: ShaharGroup-fyu
"""
# Optional module
import mdtraj as md
import numpy as np


class fengtraj:
    def __init__(self, name):
        self.name = name

    def loadsingletraj(self):
        return

    def loadmultitraj(self):
        return

    def contactCAcutoff(self):
        return

    def contactALLcutoff(self):
        return

    def averagecontact(self):
        return

    def stdcontact(self):
        return

    def plotcontact(self):
        return


def loadtrajname(pdbname, trajname):
    t = md.load(trajname, top=pdbname)
    u = t.top.select('protein')
    r = t.atom_slice(u)
    jframe = t.n_frames
    return r, jframe


def loadtraj(nn):
#    t = md.load('__traj_' + str(nn) + '.xtc', top='__START_0.pdb')
    t = md.load('__traj_' + str(nn) + '.xtc', top='__END_0.pdb')
    u = t.top.select('protein')
    r = t.atom_slice(u)
    jframe = t.n_frames
    return r, jframe


def computecontact(t, cutoffn, cutoffr):
    indexlist = []
    [cont, pairs] = md.compute_contacts(t, contacts='all', scheme='CA', ignore_nonprotein=True)
    for index, i in enumerate(pairs):
        if (i[1] - i[0]) < cutoffn:
            indexlist.append(index)
    cont = np.delete(cont, indexlist, axis=1)
    pairs = np.delete(pairs, indexlist, axis=0)
    contact = (cont < cutoffr)
    contMean = contact.mean(0)
    return pairs, contMean


def average(frame, cont, trajlist):
    nn = 0
    ntraj=len(trajlist)
    jframesum = 0
    contMeansum = []
    if nn == 0:
        contMeansum = np.empty_like(cont[nn])
    contMeansum += np.multiply(frame[nn], cont[nn])
    jframesum += frame[nn]
    if nn == ntraj - 1:
        contMean = np.true_divide(contMeansum, jframesum)
    return contMean
