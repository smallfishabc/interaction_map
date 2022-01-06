import mdtraj as md
import statistics as st
import numpy as np
# Several customized function for conformational property calculation.
# Here the x is the length of the sequence.
# Here the r is the selected trajectory, t is the raw trajectory and j is the frame number.
x=34
def calc_HB(r,j):
    hbo=md.wernet_nilsson(r)
    op=0
    HB=[]
    while op<j:
        if (len(hbo[op])):
            HB.append(len(hbo[op])/x)
        op+=1
    return(st.mean(HB))

def calc_Heli(t):
    dssp = md.compute_dssp(t)
    dssp_count = np.zeros((1, t.n_residues))
    for i in range(t.n_frames):
        for j in range(t.n_residues):
            if dssp[i,j] in 'H':
                dssp_count[0,j] += 1
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary=dssp_prob.sum(axis=1)
    # Here we should have a number as sequence length
    return(summary[0]/x)

def calc_Rg(r):
    d = md.compute_rg(r)
    return(st.mean(d))

def calc_Re(r,t):
    topology=t.topology
    rpology=topology.select_atom_indices(selection='alpha')
    d = md.compute_distances(r,[[rpology[0],rpology[-1]]])
    listtemp=[]
    for temp in d:
        listtemp.append(float(temp[0]))
    return(st.mean(listtemp))