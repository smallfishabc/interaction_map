# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:40:59 2021

@author: ShaharGroup-fyu
"""
import os
import re


def getcurrentpath():
    pwd = os.path.dirname(os.path.realpath(__file__))
    return pwd


# Automaticlly change folder and do analysis on every subfolder
def subdir(pwd, p, h):
    print(p)
    print(h)
    string = str(pwd) + '/' + h + '/' + p
    os.chdir(string)
    print(string)
    return string


def readsequence():
    # Read protein sequences
    seqopen = open('seq.fasta', 'r')
    seq = seqopen.read()
    # Clean the sequences
    while re.search('\s', seq[-1]) is not None:
        seq = seq[:-1]
    # Measure sequence length
    x = len(seq)
    return seq, x


def detect_residue(type, seq):
    target = []
    nega_charged = ['D', 'E']
    posi_charged = ['R', 'K', 'H']
    aromatic = ['F', 'Y', 'W']
    polar = ['S', 'T', 'N', 'Q']
    hydrophobic = ['A', 'V', 'I', 'L', 'M']
    hydrophilic = ['Q', 'N', 'S', 'T', 'H']
    if type == 'nega_charged':
        standard = nega_charged
    elif type == 'posi_charged':
        standard = posi_charged
    elif type == 'aromatic':
        standard = aromatic
    elif type == 'polar':
        standard = polar
    elif type == 'hydrophobic':
        standard = hydrophobic
    elif type == 'hydrophilic':
        standard = hydrophilic
    else:
        raise NameError('Wrong type')
    for index, i in enumerate(seq):
        if i in standard:
            target.append(index)
    return target

# Actuall function start here
