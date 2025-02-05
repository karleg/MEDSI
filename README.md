# MEDSI
A python implementation of the method described in https://engrxiv.org/preprint/view/2566

Requirements:  The Gurobi Optimizer (the academic license is free): https://www.gurobi.com/
                Some python modules: gurobipy, scipy, numpy, scikit-learn 

In order to execute the analysis, run the file main.py:

python main.py DATA_DIR PREFIX NETWORK

where DATA_DIR is the directory where the data files (trajectories and steady states) are stored, PREFIX is the prefix of every data file name and NETWORK is the file that contain edges between regulators and targes that will be considred in the analysis.

Examples for data files and a network files are provided in the directory data and the file network.txt
The data file contain either steady states or trajectories, where the data has already been binarized (converted to 0s and 1s)
Each column corresponds to a gene (node in the network) and each row is a time point (steady states have just one row).
A file that contains a trajectory will start with the line >trajectory , and a file that contains a steady state will start with the line >steady state .  The columns in the data file must correspond to the genes that are listed in the first line of the network file.
The network file contains all the edges that will be considered in the reconstruction, in the format: regulator   target.  Each regulator-target pair are given in a separate line.  The first line contains a list of all the genes.  It is recommended to keep the number of possible regulators of each gene as small as possible, so if a gene with completely unknown regulation is modeled, it is recommended to prioritize potential regulators first and then choose only the top ones for the network file.

To run the example provided with this repository, run:
python main.py data experiment network.txt

The file main_asyn.py and heuristic_async.py are the asynchronous equivalent of main.py and heuristic.py.  They are used in the following work: https://doi.org/10.31224/4269

for questions:  gkarleba@fitchburgstate.edu





