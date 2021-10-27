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
    def __init__(self, name, cutoff, traj=0, traj_path=0, readfromfile=0):
        # Protein name
        self.name = name
        # Cutoff distance for determining the contact
        self.cutoff = cutoff
        # Location of the protein data file. Will be modified in next update
        # self.traj_path = traj_path
        # os.chdir(traj_path)
        # If we do not want to read saved data
        if readfromfile == 0:
            # Retrieve trajectory file
            self.traj = traj
            # Calculate contact probability and residue pairs
            self.pair, self.contact = self.compute_contact()
            # Save residue pairs
            self.save_pair()
            # Save contact probability profile
            self.save_contact()
            print('finish2')
        # If we want to read saved data
        elif readfromfile == 1:
            # We do not need to pass traj file
            self.traj = 'from_read'
            # Read contact map from saved file
            self.contact = self.read_contact()
            # Read pair array from saved file
            self.pair = self.read_pair()

    # Define a function to calculate contact probability
    def compute_contact(self):
        # Calculate contact map using the MDtraj library. Cont stored the distance between two residues.
        [cont, pairs] = md.compute_contacts(self.traj, contacts='all', scheme='CA', ignore_nonprotein=True)
        print('finish1')
        # Create a empty list for removing unnecessary pairs
        index_list = []
        for index, i in enumerate(pairs):
            if (i[1] - i[0]) < 3:
                # Add unnecessary pair index into the list
                index_list.append(index)
        # Remove the unnecessary pair
        cont = np.delete(cont, index_list, axis=1)
        pairs = np.delete(pairs, index_list, axis=0)
        # Determine whether the distance is smaller than the cutoff.
        contact = (cont < self.cutoff)
        # Calculate the average contact frequency of each residue pair
        cont_mean = contact.mean(0)
        # Free up the memory
        del contact
        return pairs, cont_mean

    # Define a function to store the average contact probability into file
    def read_contact(self):
        a = np.loadtxt(str(self.name) + "contact.csv", delimiter=',')
        return a

    # Define a function to store the residue pairs into file
    def read_pair(self):
        a = np.loadtxt(str(self.name) + "pair.csv", delimiter=',')
        return a

    # Define a function to read the average contact probability from file
    def save_contact(self):
        a = self.contact
        np.savetxt(str(self.name) + "contact.csv", a, delimiter=",")
        return

    # Define a function to read the residue pairs from file
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
    # If the user only want to input individual pdb or trajectory file. Skip all other steps
    if switch == 1:
        t = md.load(traj_name, top=pdbtype)
    # Else, load multiple repeats from my simulation dataset.
    elif traj_selection == 0:
        t=load_multiple_xtc(repeats, pdbtype)
    else:
        for i in traj_selection:
            t = md.load('__traj_' + str(i) + '.xtc', top=pdbtype)
    # Select the protein from the pdb file
    u = t.top.select('protein')
    r = t.atom_slice(u)
    # Get frame number of the trajectory file.
    jframe = t.n_frames
    print(jframe)
    return r, jframe


# Designed for my multi repeat protein simulation data
def load_multiple_xtc(repeats, pdbtype):
    # Load 3 or 5 individual trajectory repeats
    if repeats == 5:
        t = md.load({'__traj_0.xtc', '__traj_1.xtc', '__traj_2.xtc', '__traj_3.xtc', '__traj_4.xtc'}, top=pdbtype)
    elif repeats == 3:
        try:
            t = md.load({'__traj_0.xtc', '__traj_1.xtc', '__traj_2.xtc'}, top=pdbtype)
        except:
            t = md.load({'__traj_0.xtc', '__traj_1.xtc', '__traj_2.xtc'}, top='__END_0.pdb')
    return t


def generate_contactmap(protein_name, read_from_file=0, cutoff=0.8, workpath=0, repeats=5):
    # If we do not want to read data from saved contact map data file. Load trajectory and calculate contact map.
    if read_from_file == 0:
        traj, frames = loadtraj(repeats)
        contact = Contactmap(protein_name, cutoff, traj)
    # If we want to load saved data. Create contact map array from saved file.
    elif read_from_file == 1:
        contact = Contactmap(protein_name, cutoff, readfromfile=read_from_file)
    print('finished')
    return contact


# Generate contactmap for single trajectory file
def generate_contactmap_single_traj(path, pdb_name, xtc_name):
    traj, frames = loadtraj(0, pdbtype=pdb_name, stdoutput=0, traj_selection=0, switch=1, traj_name=xtc_name)
    contact = Contactmap('single_traj', cutoff=0.8, traj=traj)
    return contact
