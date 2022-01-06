# This file is to create a standard file for mutation experiment.
# This standard file will store the end to end distance, radius of gyration, and helical propensity of the 'wt' simulation.
# These value will be compared to the corresponding property of the mutant.

# IMPORTANT: we need to change the length of the sequence before actually start to run the analysis

import mdtraj as md
import os
import pandas as pd
import numpy as np
import contactmapgeneration as cg
import mdtraj_function as mf

def compute_property(traj):
    # Select the protein from the pdb file
    u = traj.top.select('protein')
    r = traj.atom_slice(u)
    name=[]
    value=[]
    # Get frame number of the trajectory file.
    name.append('Frame_number')
    jframe = traj.n_frames
    value.append(jframe)
    # Calculate the radius of gyration of the trajectory
    name.append('Radius_of_Gyration(nm)')
    rg_standard = md.compute_rg(r).mean()
    value.append(rg_standard)
    # Calculate the end to end distance of the trajectory
    name.append('End_to_end_distance(nm)')
    re_standard = mf.calc_Re(r, traj)
    value.append(re_standard)
    # Calculate the Helicity of the trajectory
    name.append('Helicity')
    heli_standard = mf.calc_Heli(traj)
    value.append(heli_standard)
    storage_dict=dict(zip(name,value))
    return storage_dict

# We should pass storage_dict to this function
def write_s_file(s_dict):
    write=[]
    for key,value in s_dict.items():
        write.append(key+','+str(value)+'\n')
    with open('map_std_wt.csv', 'w') as f:
        f.writelines(write)
    f.close()

# Check whether the current directory matched with the desired directory.
def check_directory(dir_wanted):
    now=os.getcwd()
    if now==dir_wanted:
        return True
    return False

# Define a funciton to read and compare between mutant and the WT.
def compare_s_file(m_dict):
    # Read
    standard_df=pd.read_csv('map_std_wt.csv',names=['name','value'])
    standard_df['change']=np.nan
    for key,value in m_dict.items():
        value_s=standard_df.loc[standard_df['name']==key,'value']
        v_change=abs(value-value_s)
        p_change=abs(value-value_s)/value_s
        standard_df.loc[standard_df['name']==key,'percentage_change']=p_change
        standard_df.loc[standard_df['name'] == key, 'value_change'] = v_change
    standard_df.drop(columns=['value'])
    return standard_df


