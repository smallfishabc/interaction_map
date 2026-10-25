# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:40:59 2021

@author: ShaharGroup-fyu
"""
import os


# A good way to get the absolute path of the file
def getcurrentpath():
    pwd = os.path.dirname(os.path.realpath(__file__))
    return pwd


# Strip trailing whitespace from a sequence read off disk. The original did this one
# character at a time with a regex; rstrip is the same operation, and does not raise on a
# file that holds nothing but whitespace.
def _clean(seq):
    return seq.rstrip()


# Read the fasta file , clean the sequence, and measure the sequence length
def readsequence(_traj_path):
    # Read protein sequences
    with open(os.path.join(_traj_path, 'seq.fasta'), 'r') as seqopen:
        seq = seqopen.read()
    # Clean the sequence
    # Measure sequence length
    return _clean(seq)


def readsequence_single(path):
    # Read protein sequences
    with open(os.path.join(path, 'seq.txt'), 'r') as seqopen:
        seq = seqopen.read()
    # Clean the sequence
    return _clean(seq)


# Read a sequence from a plain text file or a FASTA file. Header lines are dropped and
# wrapped sequences joined, so both forms work.
def read_sequence_file(path):
    """
    :param path: file holding the sequence
    :return: the sequence, as one string
    """
    with open(path, 'r') as handle:
        lines = [line.strip() for line in handle if not line.startswith('>')]
    return ''.join(lines)


# Where to look for the sequence of a single trajectory run, in order. The simulation
# database uses seq.fasta, the single trajectory path has always used seq.txt, and the
# coarse-grained sets name the file after the construct.
def find_sequence(directory, name=None, explicit=None):
    """
    :param directory: folder holding the trajectory
    :param name: construct name, used to look for <name>.fasta
    :param explicit: a file named on the command line, tried first
    :return: (sequence, path it came from)
    :raises IOError: when none of the candidates exist
    """
    candidates = []
    if explicit:
        candidates.append(explicit if os.path.isabs(explicit)
                          else os.path.join(directory or '.', explicit))
    if name:
        candidates.append(os.path.join(directory or '.', name + '.fasta'))
    candidates.append(os.path.join(directory or '.', 'seq.txt'))
    candidates.append(os.path.join(directory or '.', 'seq.fasta'))

    for candidate in candidates:
        if os.path.exists(candidate):
            return read_sequence_file(candidate), candidate
    raise IOError('No sequence file found. Looked for: {0}'.format(', '.join(candidates)))
