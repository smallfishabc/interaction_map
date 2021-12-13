# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:40:59 2021

@author: ShaharGroup-fyu
"""
import os
import re

# A good way to get the absolute path of the file
def getcurrentpath():
    pwd = os.path.dirname(os.path.realpath(__file__))
    return pwd

# Read the fasta file , clean the sequence, and measure the sequence length
def readsequence():
    # Read protein sequences
    seqopen = open('seq.fasta', 'r')
    seq = seqopen.read()
    # Clean the sequence
    while re.search('\s', seq[-1]) is not None:
        seq = seq[:-1]
    # Measure sequence length
    x = len(seq)
    return seq, x

def readsequence_single():
    # Read protein sequences
    print(os.getcwd())
    seqopen = open('seq.txt', 'r')
    seq = seqopen.read()
    # Clean the sequence
    while re.search('\s', seq[-1]) is not None:
        seq = seq[:-1]
    return seq

# Automatically change folder and do analysis on every sub-folder
def subdir(pwd, psi, h):
    # Can be altered based on different operation system and different trajectory file structure.
    p='S_'+str(psi)
    print(p)
    print(h)
    string = str(pwd) + '/' + h + '/' + p
    os.chdir(string)
    print(string)
    return string

