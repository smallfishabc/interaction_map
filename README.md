# Intrinsically Disordered Protein(IDP) Interaction Map for Coarse-Grained simulation

Author: Feng Yu


This branch is designed based on the single_publish version.


In this branch, we focus on analyzing the intramolecular interaction of IDPs with the CALVADOS coarse-grained simulation engine designed by Kresten Lindorff-Larsen.
We will analyze the simulation trajectory and plot the strong interaction between protein residues which may impact the overall structural preference of IDPs


**A Google Colab notebook is created for this repository to help biologists design IDRs without a deep understanding of the code.**


Program structure:
default_function.py describes protocols to analyze, design, and plot IDR intramolecular interactions.
main.py provides an interface for users to separate the backend and the frontend by calling functions in the default_function.py. 


The functions in default_function.py are calling the function in the backend to generate an interaction map 
or calculate the interaction strength.

