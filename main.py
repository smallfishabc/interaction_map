# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 14:37:47 2021

@author: ShaharGroup-fyu
"""
import argparse
import os
import readpath
import showoff
import default_function


def multi_trajectory_test(args):
    test = 1
    if args.protein_directory:
        path = args.protein_directory
    elif test == 1:
        # path = 'F:\globus\simulation_sticker_spacer\F1_GS_40-summary'
        # path='F:\globus\simulation_contactmap_validation\GS18-summary'
        # path = 'F:\DATA_F\GSlinker\GS56-summary'
        # path = 'F:\DATA_F\PDBsumreal_0814\p53_new-summary'
        #path = 'F:\DATA_F\contact_mutation\p53_W53G-summary'
        #path = r'F:\DATA_F\UGDH0.5_simplified_interaction_map'
        #path = r'F:\DATA_F\puma_scrammble_sum\puma_scramble_20-summary'
        path = r'F:\DATA_F\E1A_pat-summary'
        #path = 'F:\DATA_F\puma_scramble_new\puma_scrammble_sum\puma_wildfull-summary'
        #path='/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma123/puma_wildfull-summary'
    if args.name:
        name = args.name
    elif single_traj != 1:
        try:
            name = path.split("\\")[-1].split("-")[0]
        except:
            print('err')
    if args.psi:
        psi = args.psi
    else:
        psi = '0'
    if args.restype:
        residue = args.restype
    else:
        residue = 'BB'
    if args.repeat:
        repeat=args.repeat
    else:
        repeat=5
    if single_traj != 1:
        traj_p, seq= default_function.multi_traj_pre(name, path, psi, residue)
        default_function.interaction_map_pairwise(name, traj_p, seq, traj_p, xtc_input=repeat, read_from_file=False)

if __name__ == "__main__":
    showoff.welcome()
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb","-p", help="Input PDB file")
    parser.add_argument("--xtc","-x", help="Input XTC file")
    parser.add_argument("--protein_directory", "-dir", help="protein directory")
    parser.add_argument("--single_traj", "-straj", default=1,help="Determining whether the input is a single pdb. "
                                                                  "Default value is 1."
                                                                  "Do not need to change.(test_function)")
    parser.add_argument("--repeat", "-repeat", help="Number of repeats(test_function)")
    parser.add_argument("--name", "-name", help="Name of the protein and the contact map(test_function)")
    parser.add_argument("--restype", "-restype", help="Name of Residue group, free energy altered(test_function)")
    parser.add_argument("--psi", "-psi", help="Psi value of transfer free energy(test_function)")
    args = parser.parse_args()
    single_traj = args.single_traj
    single_traj = 0
    if single_traj:
        map_name=args.name
        path = args.protein_directory
        pdb_name=args.pdb
        xtc_name=args.xtc
        sequence = readpath.readsequence_single(path)
        default_function.interaction_map_pairwise(map_name,path,sequence,path,pdb_name,xtc_name)
    else:
        multi_trajectory_test(args)
    # For multi trajectory analysis. Internal test only


