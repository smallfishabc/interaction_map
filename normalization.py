# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 18:10:07 2021

@author: ShaharGroup-fyu
"""
import numpy as np


# class chunk:
#     def __init__(self,residue,distance,name):
#         self.name=name
#         self.centerpairs=searchpairs(residue,distance)
#         self.pair_list=pair_list(pairs)
#         self.interact_list=calcinteract(self,fullpairmap,contactmapraw,plist)
#         self.interact=suminteract(self,standard)
#     def searchpairs(residue,distance):
#         return ([residue,residue+distance])
#     def pair_list(centerpairs):
#         plist=[]
#         for i in [-2,-1,0,1,2]:
#             pair=[centerpairs[0]+i,centerpairs[1]+i]
#             plist.append(pair)
#         return plist
#     def calcinteract(self,fullpairmap,contactmapraw,plist):
#         interlist=[]
#         for i in plist:
#             identify i from fullpairmap
#         retrive corresponding contactmap
#         return(pariscontact)
#     def suminteract(self,standard):
#         average=np.log(np.average(self.interact)/standard)
#         return average
def function(x, a, b):
#    if x>40:
#        x=40
    return a * x ** b


def normalization(targetmap, pairs, a1=1.65 , a2=0.4, b1=-1.4, b2=-1.16):
    interaction = np.zeros(pairs.shape[0])
    for index, i in enumerate(targetmap):
        a=a1
        b=b1
        distance=(pairs[index][1] - pairs[index][0])
        # if distance>40:
        #      a=a2
        #      b=b2
        adjustment = function(distance, a, b)
        value = 0
        if i != 0:
            value = np.log(i / adjustment)
        else:
            value = 0
        if value > 1 and i > 0.005:
            interaction[index] = 1
            print(index,distance,value,i,adjustment)
        elif value < -1 and i > 0.005:
            interaction[index] = -1
    return (interaction)

# def standardfunction(x,a,b,c):
#    return (x)
# def newnormalization(targetmap,pairs,a,b,c):
#     interaction=np.zeros(pairs.shape[0])
#     for index,i in enumerate(targetmap):
#        adjustment=standardfunction((pairs[index][1]-pairs[index][0]),a,b,c)
#        print(index)
#        value=0
#        if i!=0:
#            value=determine(i,adjustment)
#        else:
#            value=0
#        if value>1:
#            interaction[index]=1
#        elif value<-1 and i>0.005:
#            interaction[index]=-1
#     return(interaction)
# def determine(i,adjustment):
#    return(i)
