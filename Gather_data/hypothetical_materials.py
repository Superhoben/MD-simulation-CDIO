"""This is for making hypothetical materials."""
from ase import Atoms
import os
import os.path
import sys
from ase.io.trajectory import Trajectory
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')


def mix_materials(traj_name: str, add_element: str, interval=10):
    """Mix an element into existing structure.

    This function mixes add_element evenly into the atoms object in the
    trajectory file by replacing atoms in it. It creates several different
    trajectory files for different amounts of atoms added. Interval specifies
    how many atoms to add each step; high interval means few trajectory files
    will be made, low interval means many trajectory files will be made.

    When using this, the structure in the traj-file should be similar to the
    structure of the material to add, and they should work with the same
    potential, otherwise this is very unreasonable to do. Also the lattice
    constants need to be similar.

    TODO: Rescale each step to decrease the internal pressure that will be
          created from the mixing.

    Args:
        traj_name(str): Name of trajectory file of atoms object to mix
                        with add_element.
        add_element(str): Name of element to mix.
        interval(int): How many of add_element to add at each step.

    Returns:
        None

    """
    traj_path = os.path.dirname(os.path.abspath(__file__)) + \
        '/../Input_trajectory_files/'
    traj = Trajectory(traj_path+traj_name)
    atoms = traj[0]
    original_positions = atoms.get_chemical_symbols()
    len_atoms = len(atoms)
    traj_name = traj_path + traj_name.split(".")[0] + "_" + add_element + "_"
    for x in range(1, len_atoms, interval):
        atoms_positions = original_positions.copy()
        interval = len_atoms/x
        for y in range(x):
            atoms_positions[int(y*interval)] = add_element
        atoms.set_chemical_symbols(atoms_positions)
        # TODO: add scaling = optimize_scaling(atoms) here
        traj = Trajectory(traj_name+str(x)+".traj", "w")
        traj.write(atoms)
