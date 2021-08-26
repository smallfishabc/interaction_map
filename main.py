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
    #    parser.add_argument("--repeat","-repeat", help="Number of repeats")
    parser.add_argument("--restype", "-restype", help="Name of Residue group, free energy altered")
    parser.add_argument("--protein_directory", "-o", help="protein directory")
    parser.add_argument("--psi", "-psi", help="Psi value of transfer free energy")
    parser.add_argument("--database", "-db", help="Protein information file")
    parser.add_argument("--protein", "-pro", help="Protein name (if use database)")
    #   parser.set_defaults(psi=False)
    args = parser.parse_args()

    if args.protein_directory:
        path = args.protein_directory
    else:
        #path = 'F:\globus\simulation_sticker_spacer\F1_GS_40-summary'
        path='F:\globus\simulation_contactmap_validation\ChargeP_30-summary'
        #path = 'F:\DATA_F\GSlinker\GS64-summary'
    if args.psi:
        psi = args.psi
    else:
        psi = 'S_0'
    if args.restype:
        residue = args.restype
    else:
        residue = 'BB'
    default_function.interactionmapmain(path, psi, residue)
