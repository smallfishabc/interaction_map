# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:36:23 2021

@author: ShaharGroup-fyu
"""

import os

import numpy as np

import contactmapgeneration
import contactmapgeneration as cg
import contactmapplot as cm
import normalization as nl
import path
import chunk


# def average_contmean(trajlist):
#     nn = 0
#     repeats=len(trajlist)
#     contactlist=[]
#     framelist=[]
#     framesum = 0
#     print(os.getcwd())
#     for n in trajlist:
#         traj, frame = cg.loadtraj(n)
#         [pairs, contMean] = cg.computecontact(traj, 6, 0.8)
#         np.savetxt("test" + str(n) + "contact.csv", contMean, delimiter=",")
#         np.savetxt("test" + str(n) + "pairs.csv", pairs, delimiter=",")
#         if nn == 0:
#             contMeansum = np.empty_like(contMean)
#         frametest = np.multiply(frame, contMean)
#         framelist.append(frame)
#         contactlist.append(frametest)
#         np.savetxt(str(n) + "frametest.csv",frametest, delimiter=",")
#         contMeansum += np.multiply(frame, contMean)
#         framesum += frame
#         del traj
#         frame = 0
#         nn += 1
#     np.savetxt("sumcontact.csv", contMeansum, delimiter=",")
#     contMean = np.true_divide(contMeansum, framesum)
#     np.savetxt("averagecontact.csv", contMean, delimiter=",")
#     contactall = sum(contactlist)
#     frameall = sum(framelist)
#     contMean2 = []
#     for i in contactall:
#         a=i/frameall
#         contMean2.append(a)
#     np.savetxt("newaveragecontact.csv",contMean2,delimiter=",")
#     return contMean2, pairs


def interactionmap_pairwise(name,proteinpath, psi, residue):
    print(proteinpath, psi, residue)
    os.chdir(proteinpath)
    [seq, length] = path.readsequence()
    workpath = path.subdir(proteinpath, psi, residue)
    #trajlist = [0, 1, 2, 3, 4]
    contact=contactmapgeneration.generate_contactmap(name,readfromfile=0)
    location,value =chunk.mapping_chunk(name,7,contact.contact,contact.pair)
    interaction, raw_value = nl.normalization(value,location)
    #contMean, pairs = average_contmean(trajlist)
    #testchunk=chunk.Chunk(15,10,'15_10',5,contact.contact,contact.pair)
    # testchunklist=normalization.chunklist(contMean,pairs,10,5)
    # interaction = normalization.chunk_normalization(testchunklist)
    # interaction, raw_value = nl.normalization(contact.contact, contact.pair)
    # np.savetxt("interaction.csv", interaction, delimiter=",")
    # np.savetxt("raw_value.csv", raw_value, delimiter=",")
    # print(interaction)
    cm.interaction_map(seq, length, interaction, raw_value, location, 'contact_S_0', value)
# interactionmapmain('F:\globus\simulation_contactmap_validation\GS22-summary','S_0','BB')
