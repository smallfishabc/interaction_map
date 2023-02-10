# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:23:10 2021

@author: ShaharGroup-fyu
"""
# Optional module

import mdtraj as md
import pandas as pd


# We need to first load a trajectory file to use Contact map class
class ContactProbData:
    """
    Containing contact probability and the pairwise information. Will be saved to a dataframe.(Can be saved to csv)
    Traj is the trajectory file pass to the class
    Read from file can save time and memory when dealing with large protein sequences

    Attributes:
        name  The protein and the contact map name
        cutoff The cutoff distance of contact_map unit:nm
        contact 4 columns Dataframes with contact probability and residue pair information
    """

    def __init__(self, name, cutoff, traj, read_from_file=0):
        # Contact name include different information
        self.name = name
        # Cutoff distance for determining the contact
        self.cutoff = cutoff

        if read_from_file:
            # Read contact map from saved file
            self.read_contact()
        else:
            # Calculate contact probability and residue pairs
            self.contact = self.compute_contact(traj)
            self.update_csv()

    # Define a function to calculate contact probability
    def compute_contact(self, traj):
        # Calculate contact map using the MDtraj library. Cont stored the distance between two residues.
        # The scheme 'CA' and ignore_nonprotein removed ions and ACE NME caps
        [cont, pairs] = md.compute_contacts(traj, contacts='all', scheme='CA', ignore_nonprotein=True)
        # Determine whether the distance is smaller than the cutoff.
        contact = (cont < self.cutoff)
        # Calculate the average contact frequency of each residue pair
        cont_mean = contact.mean(0)
        cont_df = pd.DataFrame({'r_1': pairs[:, 0], 'r_2': pairs[:, 1], 'cont_prob': cont_mean[:]})
        cont_df['distance'] = cont_df['r_2'] - cont_df['r_1']
        # Free up the memory
        del contact
        del pairs
        return cont_df

    # Define a function to read the average contact probability dataframe from file
    def read_contact(self):
        self.contact = pd.read_csv('_'.join([str(self.name), str(self.cutoff), "contact_df_1201.csv"]), index_col=0)

    # Define a function to save the average contact probability dataframe into file
    def save_contact(self):
        self.contact.to_csv('_'.join([str(self.name), str(self.cutoff), "contact_df_1201.csv"]))

    def update_csv(self):
        self.save_contact()


# Normally, we input repeats number into this function to get the full trajectory file of the contact map
# If switch == 1, we can input single traj file.
# We can also input a list to traj_selection indicating the index of individual repeats to import part of
# the entire simulation.
# We will also create stdoutput function to calculate the standard deviation of the traj by only return a list of traj.
def load_traj_protein(pdb_top='__START_0.pdb', xtc_input='__traj_0.xtc'):
    '''
    :param pdb_top: pdb file name as topology file
    :type: string
    :param xtc_input: xtc name list o/ repeat number (will auto generate name list)
    :type: list/int
    :return: mdtraj.trajectory, trajectory framenumber
    '''
    if isinstance(xtc_input, int):
        list_xtc = ["__traj_{}.xtc".format(i) for i in range(xtc_input)]
    else:
        list_xtc = xtc_input
    t = load_xtc(list_xtc, pdb_top)
    # Select the protein from the pdb file
    u = t.top.select('protein')
    r = t.atom_slice(u)
    # Get frame number of the trajectory file.
    return r, r.n_frames


# input a xtc list to load 1 or more xtc files
def load_xtc(xtc, pdb):
    """
    Load xtc/xtc file list with correct pdb topology.
    :param xtc: a list contains xtc file names
    :type xtc: list/string
    :param pdb: pdb file name
    :type pdb: string
    :return: mdtraj.Trajectory
    """
    return md.load(xtc, top=pdb)

# What if we input a single xtc file, this need to be discuss and changed
def generate_contact(protein_name, cutoff=0.8, repeats=5, read_from_file=0):
    """
    :param protein_name: input protein_name/ contact_map_name
    :type protein_name: str
    :param cutoff: cutoff of contact_map unit:nm
    :type cutoff: float
    :param repeats: repeats number
    :type repeats: int
    :param read_from_file: read or recalculate the contact map
    :type read_from_file: bool
    :return: a Contactmap object with the corresponding protein contact map
    """
    if read_from_file:
        contact = ContactProbData(protein_name, cutoff, read_from_file)
        return contact

    traj, frames = load_traj_protein(xtc_input=repeats)
    contact = ContactProbData(protein_name, cutoff, traj)

    return contact


""" 
    def compute_contact(self, traj):
        # Calculate contact map using the MDtraj library. Cont stored the distance between two residues.
        # The scheme 'CA' and ignore_nonprotein removed ions and ACE NME caps
        [cont, pairs] = md.compute_contacts(self.traj, contacts='all', scheme='CA', ignore_nonprotein=True)
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
"""
