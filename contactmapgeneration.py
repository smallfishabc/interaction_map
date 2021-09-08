# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:23:10 2021

@author: ShaharGroup-fyu
"""
# Optional module
import os

import mdtraj as md
import numpy as np


# We need to first load a trajectory file to use Contact map class
# Containing contact map and the pairwise information. Will be saved to csv file
# Name can use to distinguish between different proteins and solution conditions
# Proteinpath can be obtained from my protein database csv file
# Traj is the trajectory file pass to the class
# Read from file can save time and memory when dealing with large protein sequences
# I have already changed working directory in the main script. We can set traj_path to alter the
# working directory if needed.
class Contactmap:
    def __init__(self, name, cutoff, traj, traj_path=0, readfromfile=0):
        self.name = name
        self.cutoff = cutoff
        self.traj = traj
        # self.traj_path = traj_path
        # os.chdir(traj_path)
        if readfromfile == 0:
            self.pair, self.contact = self.compute_contact()
            self.save_pair()
            self.save_contact()
        else:
            self.contact = self.read_contact()
            self.pair = self.read_pair()

    def compute_contact(self):
        index_list = []
        [cont, pairs] = md.compute_contacts(self.traj, contacts='all', scheme='CA', ignore_nonprotein=True)
        print(self.cutoff)
        for index, i in enumerate(pairs):
            if (i[1] - i[0]) < 3:
                index_list.append(index)
        cont = np.delete(cont, index_list, axis=1)
        pairs = np.delete(pairs, index_list, axis=0)
        contact = (cont < self.cutoff)
        np.savetxt("debug.csv", contact, delimiter=",")
        cont_mean = contact.mean(0)
        del contact
        return pairs, cont_mean

    def read_contact(self):
        a = np.loadtxt(str(self.name) + "contact.csv", delimiter=',')
        return a

    def read_pair(self):
        a = np.loadtxt(str(self.name) + "pair.csv", delimiter=',')
        return a

    def save_contact(self):
        a = self.contact
        np.savetxt(str(self.name) + "contact.csv", a, delimiter=",")
        return

    def save_pair(self):
        a = self.pair
        np.savetxt(str(self.name) + "pair.csv", a, delimiter=",")
        return


# Normally, we input repeats number into this function to get the full trajectory file of the contact map
# If switch == 1, we can input single traj file.
# We can also input a list to traj_selection indicating the index of individual repeats to import part of
# the entire simulation.
# We will also create stdoutput function to calculate the standard diviation of the traj by only return a list of traj.
def loadtraj(repeats, pdbtype='__START_0.pdb', stdoutput=0, traj_selection=0, switch=0, traj_name=0):
    # if stdoutput is 1:
    #     tlist=[]
    if switch is 1:
        t = md.load(traj_name, top=pdbtype)
    elif traj_selection is 0:
    # Need further modification
    #    for i in range(repeats):
        if repeats is 5:
            t = md.load({'__traj_0.xtc','__traj_1.xtc','__traj_2.xtc','__traj_3.xtc','__traj_4.xtc'},top=pdbtype)
        elif repeats is 3:
            t = md.load({'__traj_0.xtc','__traj_1.xtc','__traj_2.xtc'},top=pdbtype)
            # t = md.load('__traj_' + str(i) + '.xtc', top=pdbtype)
            # if stdoutput is 1:
            #     tlist.append(t)
    else:
        for i in traj_selection:
            t = md.load('__traj_' + str(i) + '.xtc', top=pdbtype)
    # if stdoutput is 1:
    #     a = stdcalculation()
    #     np.savetxt('stdoutput.csv', a, delimiter=",")
    u = t.top.select('protein')
    r = t.atom_slice(u)
    jframe = t.n_frames
    return r, jframe


def generate_contactmap(protein_name, repeats=5, cutoff=0.8, workpath=0):
    traj, frames = loadtraj(repeats)
    contact = Contactmap(protein_name, cutoff, traj)
    return contact
