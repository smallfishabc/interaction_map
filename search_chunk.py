# This file is designed to search interactions after identified the strongest interaction.
# 1. Search what interaction happens simultaneously with the interaction
# 2. Compare of radius of gyration and helicity distribution with/without the interaction
import mdtraj as md
import pandas as pd
import numpy as np
import contact_map_generation as cg

# Get frame numbers with the target interaction
def inter_frames(r1,r2,traj,chunk=True,cutoff=0.8):
    if chunk:
        contact_pairs=[[x,y] for x in range(r1-1,r1+2) for y in range(r2-1,r2+2)]
    else:
        contact_pairs=[[r1, r2]]
    [cont, pairs] = md.compute_contacts(traj, contacts=contact_pairs, scheme='ca', ignore_nonprotein=True)
    contact = (cont < cutoff)
    contact_select = np.any(contact, axis=1)
    return contact_select

# Use r[contact_select] to call the filtered trajectory

# Search other interactions with the strongest interaction.


# Compare radius of gyration and helicity

