# pycell
PyMOL plugin that can draw supercells on each frame from a pdb trajectory. Currently it has been tested you PDB files created with GROMACS and Materials Studio. The cell is drawn following the convention, A along X, B in YX plane.
Pycell requires that each PDB model have the CRYST1 field. Also the NumPy package for python is needed.

Tested with 
python 2.7.6
GROMACS 5.1.1
PyMOL 1.7.4.5 Edu
Materials Studio 6.0




