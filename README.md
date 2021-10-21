# inteactionmap_single_publish
Author:Feng Yu
This branch is designed based on single_publish version.

In this branch, we have two types of functions.

The first type of functions is designed for analyzing single pdb/xtc file. 

The second type of functions is designed for multi-protein analysis based my simulation data structure.

We will focus on explaining the first type of functions in this readme file.

Basic usage of the database:
python3 main.py --pdb pdb_file.pdb --xtc xtc_file.xtc -dir F:\DATA
--pdb or -p is the name of the pdb file
--xtc or -x is the name of the xtc file
--protein_directory option or -dir is used to identify the location of data file
Please put your sequence as the first line in seq.txt text file under the protein_directory
Interaction map will be saved as pdb_file.png under protein_directory

File structure:
Our main script will usually called the function in default_function.py. Therefore, the main function will
not directly interact with the backend.

The functions in default_function.py are calling the function in the backend to generate a interaction map 
or calculte the interaction strength.

The other files are storing functions having different purposes. We will explain some of them in the
following paragraph.

class Contactmap
For a contact map calculated by mdtraj, it will give out a pair array and a contact probability array.
During the calculation process, we will do some modification or remove unnecessary pairs from the array.
In addtion, we will also generate corresponding interaction map array. To prevent data mismatch, we created 
Contactmap. We will further intergrate the operation and function into this class to make this program easy to use.

Function explaination

generate_contactmap_single_traj

We will load the trajectory file, generate contact map data and store it into the Contactmap class in this function.

normalization

We will use GS-linker standard curve to calculate the interaction in this function

interaction_map

In this function, we will calculate the interaction strength and plot the interaction based on pyplot and networkx library.
