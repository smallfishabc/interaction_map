import os
import numpy as np
import mdtraj as md
import contactmapgeneration as cg
import mdtraj_function as mf
import create_standard as cs
import pandas as pd
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
        #print(i[index_wanted] == True)
        if i[index_wanted] == True:
            # Attach the contact to the final list
            index_list.append(index)
    return index_list


def frame_selection(index_list, trajectory_object):
    # We'd better pass a object using the load traj function
    target_traj=trajectory_object.slice(index_list)
    # Optional: Save the new traj to a separate pdb file
    target_traj.save_xtc('test.xtc')
    return target_traj

def comparison(target_pair,traj,contact,pairs,limit):
    select_list=select_conformation_index_interaction(target_pair,contact,pairs)
    if len(select_list)<limit:
        return (0,0,0)
    target_traj=frame_selection(index_list=select_list,trajectory_object=traj)
    m_dict=cs.compute_property(target_traj)
    frame_number=m_dict['Frame_number']
    compare_result = cs.compare_s_file(m_dict)
    ree_change=compare_result.loc[compare_result['name']=='End_to_end_distance(nm)','value_change'].iloc[0]
    helicity_change=compare_result.loc[compare_result['name']=='Helicity','value_change'].iloc[0]
    return frame_number,ree_change,helicity_change

# Based on previous function. Systematically quantify how interaction influence the ensemble property

if __name__ == "__main__":
    target_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum/puma_wildfull-summary/BB/S_0'
    os.chdir(target_directory)
    repeats=5
    # This part will be needed if the second part is removed
    read_from_file=1
    name='puma_wildfull'

    traj = cg.load_multiple_xtc(repeats, pdbtype='__START_0.pdb')
    s_dict=cs.compute_property(traj)
    cs.write_s_file(s_dict)
    # Can be saved to a individual function and removed
    u = traj.top.select('protein')
    r = traj.atom_slice(u)
    # Get frame number of the trajectory file.
    jframe = traj.n_frames
    print('end')
    cutoff = 0.8
    [cont, pairs] = md.compute_contacts(traj, contacts='all', scheme='CA', ignore_nonprotein=True)
    contact = (cont < cutoff)
    # Remove end
    df=pd.read_csv('att1interaction_pair.csv', sep=',',header=None)
    #df2=pd.read_csv('att2interaction_pair.csv',sep=',',header=None)
    #test=df.values.tolist()+df2.values.tolist()
    test = df.values.tolist()
    if not cs.check_directory(target_directory):
        raise Exception('Wrong_location')
    result_value=pd.DataFrame(columns=['interaction_pair','frame_number','ree_change','helicity_change','ree_weighted','helicity_weighted'])
    for i in test:
        if i[1]-i[0]<=5:
            continue
        frame_number,ree_change,helicity_change=comparison(i,traj,contact,pairs,limit=1400)
        if frame_number==0:
            continue
        ree_weighted=frame_number*ree_change
        helicity_weighted=frame_number*helicity_change
        result_value=result_value.append({'interaction_pair':str(i[0])+'_'+str(i[1]),
                             'frame_number':frame_number,
                             'ree_change':ree_change,
                             'helicity_change':helicity_change,
                             'ree_weighted':ree_weighted,
                             'helicity_weighted':helicity_weighted},ignore_index=True)
    result_value.to_csv('mutation.csv',index=False)

    # select_list=select_conformation_index_interaction([13,23],contact,pairs)
    # target_traj=frame_selection(index_list=select_list,trajectory_object=traj)
    # m_dict=cs.compute_property(target_traj)
    # compare_result=cs.compare_s_file(m_dict)
    # if cs.check_directory(target_directory):
    #     pass
    # print(compare_result)
'''
    # test purpose only
    target_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma123/puma_wildfull-summary/BB/S_0'
    os.chdir(target_directory)
    print('start')
    traj = md.load({'__traj_0.xtc', '__traj_1.xtc', '__traj_2.xtc', '__traj_3.xtc', '__traj_4.xtc'}, top='__START_0.pdb')
    # Select the protein from the pdb file
    u = traj.top.select('protein')
    r = traj.atom_slice(u)
    # Get frame number of the trajectory file.
    jframe = traj.n_frames
    print('end')
    cutoff = 0.8
    [cont, pairs] = md.compute_contacts(traj, contacts='all', scheme='CA', ignore_nonprotein=True)
    contact = (cont < cutoff)
    a=select_conformation_index_interaction([13,23],contact,pairs)
    target_traj=frame_selection(index_list=a,trajectory_object=traj)
    target_u = target_traj.top.select('protein')
    target_r = target_traj.atom_slice(u)
    target_jframe = target_traj.n_frames
    rg=md.compute_rg(target_r).mean()
    rg_standard=md.compute_rg(r).mean()
    re=mf.calc_Re(target_r,target_traj)
    re_standard=mf.calc_Re(r,traj)
    helicity=mf.calc_Heli(target_traj)
    helicity_standard=mf.calc_Heli(traj)
    HB=mf.calc_HB(target_r,target_jframe)
    HB_standard=mf.calc_HB(r,jframe)
    print('Number of frames of selected trajectory: ',target_jframe)
    print('Number of frames of full trajectory: ', jframe)
    print('Radius of Gyration(nm) of selected trajectory:   ',rg)
    print('Radius of Gyration(nm) of full trajectory:   ',rg_standard)
    print('End to End Distance(nm) of selected trajectory:   ',re)
    print('End to End Distance(nm) of full trajectory:   ',re_standard)
    print('Helicity of selected trajectory:   ',helicity)
    print('Helicity of full trajectory:   ',helicity_standard)
    print('Hydrogen bonds per residue of selected trajectory:   ',HB)
    print('Hydrogen bonds per residue of full trajectory:   ',HB_standard)
'''
