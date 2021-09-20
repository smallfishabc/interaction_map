# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 14:37:47 2021

@author: ShaharGroup-fyu
"""
import argparse
import showoff
import default_function

if __name__ == "__main__":
    showoff.welcome()
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat","-repeat", help="Number of repeats")
    parser.add_argument("--name", "-name", help="Name of the protein and the contact map")
    parser.add_argument("--restype", "-restype", help="Name of Residue group, free energy altered")
    parser.add_argument("--protein_directory", "-o", help="protein directory")
    parser.add_argument("--psi", "-psi", help="Psi value of transfer free energy")
    parser.add_argument("--single_traj", "-straj", help ="If you want to only process specific protein. Set this value "
                                                         "to 1. The program will try to find the traj file and pdb file under "
                                                         "the protein directory" )
    args = parser.parse_args()
    single_traj=args.single_traj
    test=1
    if args.protein_directory:
        path = args.protein_directory
    elif test is 1:
        # path = 'F:\globus\simulation_sticker_spacer\F1_GS_40-summary'
        # path='F:\globus\simulation_contactmap_validation\GS44-summary'
        # path = 'F:\DATA_F\GSlinker\GS56-summary'
        path = 'F:\DATA_F\puma_scramble_new\puma123\puma_scramble_3-summary'
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
        psi = 'S_0'
    if args.restype:
        residue = args.restype
    else:
        residue = 'BB'
    if single_traj is not 1:
        default_function.interactionmap_pairwise(name, path, psi, residue)
    elif single_traj is 1:
        default_function.single_traj_contactmap()


