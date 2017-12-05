# PyRAMA
Python3 implementation of the Ramachandran plot

Usage:

PyRAMA.py my_pdb_file.pdb

Note: The script is able to read in multiple PDB files, however all of the torsion angles will be displayed on the same plot, with the same color!

Dependencies:

Running PyRAMA requires *matplotlib* and *biopython*

To install these on a standard Linux system:

    apt install python3-matplotlib
    apt install python3-biopython

For the standard PSI and PHI preferences see:

Lovell *et al*. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. 2003
DOI: 10.1002/prot.10286