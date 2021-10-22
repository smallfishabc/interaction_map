import os
import numpy as np
import mdtraj as md
import contactmapgeneration as cg

# This module will pick interested frames from my simulation trajectory.
# The standard may be based on interaction or residual structure.

# input will be a  contact map (It is important to release memory after the running)
# If we want to select the conformational ensemble based on interaction map, we will also need to pass interaction map
# If we want to select based on the residue structure, we will need to pass a helicity array.
# The select function will only output a list of frame number.
# We will have a follow up function to retrieve the frames from the trajectory and create a separate pdb file or
# mdtraj object for further comparison with the entire trajectory.

def select_conformation_index_interaction(target_pair,contact_map, pair_index):
    # find index of target_pair from pair_index array
    index_wanted = np.where(np.all(pair_index == target_pair, axis=1))[0][0]# Need to be tested.
    # Search for interaction in the interaction_array
    # Create for a list to store all desired frames
    index_list=[]
    # loop over the contact_map
    for index, i in enumerate(contact_map):
        # If the contact between target_pair exist
        print(i[index_wanted] == True)
        if i[index_wanted] == True:
            # Attach the contact to the final list
            index_list.append(index)
    return index_list


def frame_selection(index_list, trajectory_object):
    # We'd better pass a object using the load traj function
    target_traj=trajectory_object.slice(index_list)
    # Optional: Save the new traj to a separate pdb file
    target_traj.save_pdb('test.pdb')
    return target_traj

# Based on previous function. Systematically ser

if __name__ == "__main__":
    # test purpose only
    target_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma123/puma_wildfull-summary/BB/S_0'
    os.chdir(target_directory)
    traj = md.load('__traj_0.xtc', top='__START_0.pdb')
    # Select the protein from the pdb file
    u = traj.top.select('protein')
    r = traj.atom_slice(u)
    # Get frame number of the trajectory file.
    jframe = traj.n_frames
    cutoff = 0.8
    [cont, pairs] = md.compute_contacts(traj, contacts='all', scheme='CA', ignore_nonprotein=True)
    contact = (cont < cutoff)
    a=select_conformation_index_interaction([13,23],contact,pairs)
    target_traj=frame_selection(index_list=a,trajectory_object=traj)
    target_u = target_traj.top.select('protein')
    target_r = target_traj.atom_slice(u)
    rg=md.compute_rg(target_r).mean()
    rg_standard=md.compute_rg(r).mean()
    print(rg,rg_standard)
