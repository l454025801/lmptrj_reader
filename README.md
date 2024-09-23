Python scripts used to analyze lammpstrj with fix bond/create. 
This scripts require lammp dump file containing bonds and output new bonds formed on each timestep.
The dump file should have the following format

    ITEM: TIMESTEP
    THE ACTUAL TIMESTEP
    ITEM: BOX BOUNDS pp pp pp
    x_low x_high
    y_low y_high
    z_low z_high
    ITEM: ENTRIES index c_bond[1] c_bond[2] c_bond[3]
    index_1 bond_atom_1 bond_atom_2 bond_type
    index_2 bond_atom_1 bond_atom_2 bond_type
    ...

This file can be created with lammp input command 
  dump dump_name all local file_name index c_bond[*]
    
The scripts support grand canonical simulation (changing number of atoms).

Here the scripts are used to analyze a lammps simulation of fiber-antibody crosslinks. You need to adapt the script for other analysis.
