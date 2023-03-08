# This file is designed to search interactions after identified the strongest interaction.
# 1. Search what interaction happens simultaneously with the interaction
# 2. Compare of radius of gyration and helicity distribution with/without the interaction
import mdtraj as md
import pandas as pd
import numpy as np
import statistics as st
import contact_map_generation as cg

# Get frame numbers with the target interaction
def inter_frames(r1,r2,traj,chunk=True,cutoff=0.8):
    if chunk:
        contact_pairs=[[x,y] for x in range(r1-1,r1+2) for y in range(r2-1,r2+2)]
    else:
        contact_pairs=[[r1, r2]]
    [cont, pairs] = md.compute_contacts(traj, contacts=contact_pairs, scheme='ca', ignore_nonprotein=True)
    contact = (cont < cutoff)
    contact_select = np.any(contact, axis=1)
    return contact_select

# Use r[contact_select] to call the filtered trajectory

# Search other interactions with the strongest interaction.
def new_interaction_map(name,traj,contact_select,cutoff=0.8):
    contact = cg.ContactProbData(name+'with', cutoff, traj[contact_select])
    contact_non = cg.ContactProbData(name+'without', cutoff, traj[~contact_select])
    return contact, contact_non


# Compare radius of gyration and helicity

# Calculate the Hydrogen-bond number using Mdtraj library
# Here the r is the trajectory that only containing the protein atoms
# j is the number of frames in the trajectory.
def calc_HB(r,j,length):
    #Calculate the Hydogen_bond
    hbo=md.wernet_nilsson(r)
    # This is a couter for the loop over frames
    op=0
    # A list to store the final result.
    HB=[]
    # looping over each frame
    while op<j:
        # hbo[op] is a sub numpy array stored donor and acceptor
        # information of the H-bonds. The length of the array is the number of
        # H-bonds
        if (len(hbo[op])):
            HB.append(len(hbo[op])/length)
        op+=1
    # Average H-bonds per residue per frame
    return(sum(HB)/j)

# Calculate the helical propensity fo the protein.
# Here the t is the full trajectory file(We can try to use r instead)
def calc_Heli(t,length):
    # Calculate the secondary structure using dssp
    dssp = md.compute_dssp(t)
    # Create an empty array with the row number equals to the number of residues
    # We may not need the number 1 here
    dssp_count = np.zeros((1, t.n_residues))
    # Loop over frames
    for i in range(t.n_frames):
        # Loop over residues
        for j in range(t.n_residues):
            # In dssp, letter 'H' means alpha-helix
            if dssp[i,j] in 'H':
                dssp_count[0,j] += 1
    # Get the average of the dssp_count
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary=dssp_prob.sum(axis=1)
    return(summary[0]/length)

#calc_beta is added 03/28/2022
# Same algorithm as alpha
def calc_Beta(t,length):
    dssp = md.compute_dssp(t)
    dssp_count = np.zeros((1, t.n_residues))
    for i in range(t.n_frames):
        for j in range(t.n_residues):
            if dssp[i,j] in 'E':
                dssp_count[0,j] += 1
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary = dssp_prob.sum(axis=1)
    return(summary[0]/length)

# Average radius of gyration of the ensemble
def calc_Rg(r):
    d = md.compute_rg(r)
    return(st.mean(d))

# Average end to end distance of the ensemble
# The defination of the end to end distance is the distance between
# the first alpha carbon and the last alpha carbon of the protein sequence.
def calc_Re(r,t):
    topology=t.topology
    rpology=topology.select_atom_indices(selection='alpha')
    d = md.compute_distances(r,[[rpology[0],rpology[-1]]])
    listtemp=[]
    for temp in d:
        listtemp.append(float(temp[0]))
    return(st.mean(listtemp))

def calculate_ensemble(traj, sliced_traj, length):
    return {'Rg':calc_Rg(sliced_traj),'Re':calc_Re(sliced_traj,traj),'Heli':calc_Heli(traj,length),'Beta':calc_Beta(traj,length),'HB':calc_HB(traj,traj.n_frames,length)}
def compare_ensemble(traj, sliced_traj, length, contact_select):
    return {'interaction':calculate_ensemble(traj[contact_select],sliced_traj[contact_select],length),
           'non_interaction':calculate_ensemble(traj[~contact_select],sliced_traj[~contact_select],length)}