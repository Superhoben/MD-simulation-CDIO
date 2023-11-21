"""This is for making hypothetical materials."""
import os
import sys
import random
import os.path
from ase import Atoms
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
    original_symbols = atoms.get_chemical_symbols()
    nr_of_atoms = len(atoms)
    traj_name = traj_path + traj_name.split(".")[0] + "_" + add_element + "_"
    for x in range(1, nr_of_atoms, interval):
        atoms_positions = original_symbols.copy()
        random_mix(atoms_positions, add_element, x)
        atoms.set_chemical_symbols(atoms_positions)
        # TODO: add scaling = optimize_scaling(atoms) here
        traj = Trajectory(traj_name+str(x)+".traj", "w")
        traj.write(atoms)


def random_mix(chemical_symbols: list, add_element: str, amount: int):
    """Switch a given number of instances from a given list to a given string.

    This can be used to create 'alloys': use atoms.get_chemical_symbols to
    get the chemical symbols of the atoms object, then send this into this
    function together with the element to mix with, and how many to add.
    Then use atoms.set_chemical_symbols to set the new chemical symbols.

    For this to make sense, the added element needs to be similar to the
    host material in terms of structure and size.

    Args:
        chemical_symbols(list[str]): List of chemical symbols.
        add_element(str): The element to add to the list.
        amount(int): The amount of elements to add to the list.

    Returns:
        chemical_symbols(list[str]): The new list with randomly inserted
            elements.
    """
    insert_indexes = random.sample(range(0, len(chemical_symbols)), amount)
    for index in insert_indexes:
        chemical_symbols[index] = add_element

    return chemical_symbols


def create_mix_from_concentration(traj_name: str, add_element: str, concentration):
    """Mix an element into existing structure from concentration.
    
    This creates one trajectory file by mixing add_element into the structure
    in the input file with the given concentration.
    
    Args:
        traj_name(str): Name of trajectory file of atoms object to mix
                        with add_element.
        add_element(str): Name of element to mix.
        interval(int): The desired concentration of the add_element.
    """
    traj_path = os.path.dirname(os.path.abspath(__file__)) + \
        '/../Input_trajectory_files/'
    traj = Trajectory(traj_path+traj_name)
    atoms = traj[0]
    original_symbols = atoms.get_chemical_symbols().copy()
    atoms.set_chemical_symbols(random_mix(original_symbols, add_element,
                               int(concentration*len(atoms))))
    # TODO: optimize_scaling here?
    traj_name = traj_path + traj_name.split(".")[0] + "_" + add_element + "_" + str(concentration) + ".traj"
    traj = Trajectory(traj_name, "w")
    traj.write(atoms)
