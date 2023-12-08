"""This is for making hypothetical materials."""
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
import shutil
import random
from ase import Atoms
from asap3 import EMT
from Simulation.lattice_constant import optimize_scaling
from ase.io.trajectory import Trajectory


def mix_materials(traj_name: str, add_element: str, target_concs, mix_dir_name=False):
    """Mix an element into an existing structure at diffrent concentrations

    This function mixes add_element randomly into the atoms object in the trajectory file 
    by replacing atoms in it. It creates several different trajectory files for different 
    amounts of atoms added. 

    When using this, the structure in the traj-file should be similar to the
    structure of the material to add, and they should work with the same
    potential, otherwise this is very unreasonable to do. Also the lattice
    constants should be similar if the original structure is supposed to be preserved.

    Args:
        traj_name(str): Name of trajectory file of atoms object to mix with add_element.
        add_element(str): Name of element to mix.
        target_concs(list(float)): The concentration which will be aimed at for each mixing.
        mix_dir_name(str or bool): Should be a string with the directory in which we want to
            save the mixed materials, if we don't want to save them it can be set to false.

    Returns:
        mixed_atoms_list(list(Atoms)): A list with the mixed materials created.
        actual_concs(list(float)): A list with the final concetration of the mixed in element.
    """
    input_traj_dir_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
    traj = Trajectory(input_traj_dir_path+traj_name)
    atoms = traj[0]

    if mix_dir_name:
        alloy_dir = input_traj_dir_path+mix_dir_name
        if os.path.exists(alloy_dir):
            shutil.rmtree(alloy_dir)
        os.mkdir(alloy_dir)

    original_symbols = atoms.get_chemical_symbols()
    mixed_atoms_list = []
    actual_concs = []
    for c in target_concs:
        mixed_symbols = original_symbols.copy()
        number_of_symbols_to_switch = round(c*len(atoms))
        c = number_of_symbols_to_switch/len(atoms)
        actual_concs.append(c)
        random_mix(mixed_symbols, add_element, round(c*len(atoms)))
        atoms_copy = atoms.copy()
        atoms_copy.set_chemical_symbols(mixed_symbols)
        atoms_copy.calc = EMT()
        suggested_scaling = optimize_scaling(atoms_copy)
#        print(suggested_scaling)
        if suggested_scaling != 0:
            atoms_copy.set_cell(atoms_copy.cell*suggested_scaling, scale_atoms=True)
        if mix_dir_name:
            traj = Trajectory(alloy_dir+'/'+str(0.1*round(1000*c))+add_element+'_in_'+traj_name.split("/")[-1], "w")
            traj.write(atoms_copy)
        mixed_atoms_list.append(atoms_copy)

    return mixed_atoms_list, actual_concs


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
        concentration(float): The desired concentration of add_element.
                              Should be between 0 and 1.
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

if __name__ == "__main__":
    mix_list, actual_conc = mix_materials("1728_atoms_of_mp-30.traj", "Ag", [0.2, 0.5, 0.8], "Cu_Ag_mixes")
    print(actual_conc)
