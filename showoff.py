# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 14:03:38 2021

@author: ShaharGroup-fyu
"""


##Get from Alex
def welcome():
    print("")
    print("#################################################################")
    print("")
    print(".................................................................")
    print("")
    print(".....................IDP Interaction Map.........................")
    print("")
    print(".................................................................")
    print("")
    print("#################################################################")
    print("Some function are labeled with test_function")
    print("These function are designed for developer's simulation database")
    print("For single pdb/xtc analysis.")
    print("Here is an example code")
    print("python3 script --pdb pdb_file.pdb --xtc xtc_file.xtc -dir F:\DATA")
    print("--pdb or -p is the name of the pdb file")
    print("--xtc or -x is the name of the xtc file")
    print("--protein_directory option or -dir is used to identify the location of data file")
    print("Please put your sequence as the first line in seq.txt text file under the same protein_directory")
    print("Interaction map will be saved as pdb_file.png under protein_directory")

# Reference code for using parser
# if __name__=="__main__":
#    import argparse 
#    import sys
#    parser = argparse.ArgumentParser()
###Show
#    welcome()
#
#    parser.add_argument("--pdb","-pdb", help="Input PDB file") 
#    parser.add_argument("--xtc","-xtc", help="Input XTC file") 
#    parser.add_argument("--output_directory","-o", help="Output directory") 
#    parser.add_argument("--verbose","-v", help="Be loud and obnoxious", action='store_true')
#    parser.add_argument("--stride", help="Number of frames to extract [D=1]")
