# inteactionmap_single_publish
Author:Feng Yu
This branch is designed based on single_publish version.

In this branch, we have two types of functions.

The first type of functions is designed for analyzing single pdb/xtc file. 

The second type of functions is designed for multi-protein analysis based my simulation data structure.

We will focus on explaining the first type of functions in this readme file.

Basic usage of the database:
##copy from the main

File structure:
Our main script will usually called the function in default_function.py. Therefore, the main function will
not directly interact with the backend.

The functions in default_function.py are calling the function in the backend to generate a interaction map 
or calculte the interaction strength.

The other files are storing functions having different purposes. We will explain some of them in the
following paragraph.

Function explaination
