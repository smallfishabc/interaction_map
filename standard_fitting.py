import statistics

import numpy as np
import scipy
import os
import pandas
import chunk
import contactmapgeneration
import normalization
import statistics


# def pairwise_fitting():
#    return
def directory_location():
    os.chdir('F:\globus\simulation_contactmap_validation')
    df = pandas.read_csv('database_entry.csv')
    location = df['Directory'].tolist()
    return location


def chunk_fitting(seq_length, y):
    cutoff = 0
    x = [i for i in range(cutoff, seq_length)]
    popt, pcov = scipy.optimize.curve_fit(normalization.fitting_function, x, y)
    a = popt[0]
    b = popt[1]
    return a, b


# Generate contact probability with different distance and chunksize
def generate_chunk_contact_probability(distance, contactmap, contact_probability, chunknumber,chunksize=5):
    # Create two empty list. Every time we find a new
    chunkobject, needed = chunk.mapping_chunk(contactmap.name, distance, contactmap.contact, contactmap.pair)
    chunknumber += len(needed)
    contact_probability.append(needed)
    return contact_probability, chunknumber


# Get the path of selected protein
def get_contact_map_list(location, repeats=5):
    contactmaplist = []
    for i in location:
        protein_name = i.split('\\')[-1].split('-')[0]
        fullpath = os.path.join(i, 'BB', 'S_0')
        b = contactmapgeneration.generate_contactmap(protein_name, workpath=fullpath)
        print(protein_name)
        yield b


def full_probability_set(limit, location):
    contactmaplist = get_contact_map_list(location)
    ii=4
    i=4
    full_set = [[] for x in range(limit-ii)]
    chunk_set = [0 for x in range(limit-ii)]
    for k in contactmaplist:
        while i < limit:
            contactmaplist, average = generate_chunk_contact_probability(i, k,contact_probability=full_set[i-ii],chunknumber=chunk_set[i-ii])
            i+=1
    return full_set


if __name__ == "__main__":
    location = directory_location()
    full_probability_set(6,location)
