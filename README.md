# MoleculeTopologyTracking

This is a Python script that tracks the topological variation of the molecules in a system.

# Prerequisite

    Numpy
    ASE             https://wiki.fysik.dtu.dk/ase/

# Usage

    python topology <trajectory_dir> <output_dir>

where <trajectory_dir> is a folder containing a list of *.extxyz trajectory files and <output_dir> is the folder that stores the ouput trajectories remaining each *.extxyz files with *-cluster.extxyz. Each output trajectory contains one additional per-atom value, "cluster", which is 1 if the molecule to which the atom belongs has a topology change between the initial and the final snapshot and 0 otherwise.
