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
    elif test is 1:
        # path = 'F:\globus\simulation_sticker_spacer\F1_GS_40-summary'
        # path='F:\globus\simulation_contactmap_validation\GS44-summary'
        # path = 'F:\DATA_F\GSlinker\GS56-summary'
        # path = 'F:\DATA_F\puma_scramble_new\puma123\puma_scramble_3-summary'
        path='/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma123/puma_wildfull-summary'
    if args.name:
        name = args.name
    elif single_traj is not 1:
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
    if single_traj is not 1:
        # if test == 1:
        #     default_function.test_function(name,path,residue)
        # else:
        default_function.interactionmap_pairwise(name, path, psi, residue)
    # elif single_traj is 1:
    #     default_function.single_traj_contactmap()

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
    single_traj = 0
    if single_traj==1 :
        path = args.protein_directory
        pdb_name=args.pdb
        xtc_name=args.xtc
        sequence=readpath.readsequence_single()
        os.chdir(path)
        default_function.interactionmap_pairwise_single_traj(path,pdb_name,xtc_name,sequence)
    else:
        multi_trajectory_test(args)
    # For multi trajectory analysis. Internal test only


