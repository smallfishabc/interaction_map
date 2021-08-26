# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:40:59 2021

@author: ShaharGroup-fyu
"""
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#Optional module
import mdtraj as md
import statistics as st
import re
import networkx as nx

def getcurrentpath():
    pwd=os.path.dirname(os.path.realpath(__file__))
    return(pwd)
def readsequence():
#Read protein sequences
    seqopen=open('seq.fasta','r')
    seq=seqopen.read()
#Clean the sequences
    while (re.search('\s',seq[-1]) != None):
        seq=seq[:-1]
#Measure sequence length
    x=len(seq)
    return(seq,x)
# Automaticlly change folder and do analysis on every subfolder
def subdir(pwd,p,h='BB'):
    string = str(pwd)+'/'+h+'/'+p
    os.chdir(string)

#Actuall function start here