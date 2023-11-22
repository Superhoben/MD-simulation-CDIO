import os
from math import floor, cbrt
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.visualize import view
from ase.build import molecule, nanotube, bulk, make_supercell
from ase.lattice.triclinic import Triclinic
from ase.lattice.monoclinic import SimpleMonoclinic, BaseCenteredMonoclinic
from ase.lattice.hexagonal import Hexagonal
from ase.lattice.tetragonal import SimpleTetragonal, CenteredTetragonal
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic, SimpleCubic
from ase.lattice.orthorhombic import SimpleOrthorhombic, BaseCenteredOrthorhombic, BodyCenteredOrthorhombic, FaceCenteredOrthorhombic


def save_any_atoms(file_name, target_number_of_atoms=300, *args, **kwargs):
    """Use the Atoms class initializer and save an atoms object to a .traj file

    Args:
        file_name(string): Name of the .traj file that will be created
        target_number_of_atoms(int): Determines size of the created super cell
        *args: An arbitary number of arguments which will be passed forward to the
            Atoms class initializer
        **kwargs: An arbitary number of keyword arguments which will be passed forward
            to the Atoms class initializer

    Return:
        None
    """
    primitive_cell = Atoms(args, kwargs)
    number_atoms_primitive = primitive_cell.get_number_of_atoms()
    n = floor(cbrt(target_number_of_atoms/number_atoms_primitive))
    M = [[n, 0, 0], [0, n, 0], [0, 0, n]]
    atoms = make_supercell(primitive_cell, M)
    size_descripition = str(atoms.get_number_of_atoms())+'_atoms_of_'
    path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
    location_and_name = path_to_traj_folder + size_descripition + file_name + '.traj'
    traj = Trajectory(location_and_name, "w")
    traj.write(atoms)


def create_view_and_save_crystal_guided():
    """Guide the user through creating a crystal .traj input file through the terminal
    """
    primitive_cell = create_crystal_guided()
    number_atoms_primitive = primitive_cell.get_number_of_atoms()
    print("This is the primitive cell you've created")
    view(primitive_cell, block=False)
    target_number_of_atoms = input("Input the target number of atoms in the supercell: ")
    input_error_message = "Must be an integer bigger above " + str(primitive_cell.get_number_of_atoms()) + ": "
    while True:
        if target_number_of_atoms.isdigit():
            target_number_of_atoms = int(target_number_of_atoms)
            if target_number_of_atoms > primitive_cell.get_number_of_atoms():
                break
        target_number_of_atoms = input(input_error_message)
    file_name = input("Input of the name .traj file that will contain the material: ")
    n = floor(cbrt(target_number_of_atoms/number_atoms_primitive))
    M = [[n, 0, 0], [0, n, 0], [0, 0, n]]
    atoms = make_supercell(primitive_cell, M)
    path_to_traj_folder = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
    location_and_name = path_to_traj_folder + file_name + '.traj'
    traj = Trajectory(location_and_name, "w")
    traj.write(atoms)


def create_crystal_guided():
    """Provides guided creation of primitive cell by the ase.lattice crystal builders
    """
    print("Choose crystal system. System names and description of the lattice vectors follows below. \n",
          "1 - Triclinic, arbitrary lengths and directions \n",
          "2 - Monoclinic, arbitrary lengths, two vectors are orthagonal \n",
          "3 - Orthorhombic, arbitrary lengths, all vectors are orthagonal \n",
          "4 - Tetragonal, two equal lenghts, all vectors are orthagonol \n",
          "5 - Trigonal, not yet implemented \n",
          "6 - Hexagonal, not yet implemented \n",
          "7 - Cubic, equal lenghts and all vectors are orthagonal")
    crystal_system = input("Input number of choosen crystal system: ")
    while not (crystal_system in ['0', '1', '2', '3', '4', '7']):
        crystal_system = input("Invalid input, enter an integer between 1 and 7 except 5 and 6: ")
    crystal_system = int(crystal_system)
    if crystal_system == 1:
        lengths_and_angles = {}
        lengths_and_angles['a'] = float(input("Input length in of the first lattice translation vector in Å: "))
        lengths_and_angles['b'] = float(input("Input length of the second lattice translation vector in Å: "))
        lengths_and_angles['c'] = float(input("Input length of the third lattice translation vector in Å: "))
        lengths_and_angles['alpha'] = float(input("Input angle between the first and second translation vector in degrees: "))
        lengths_and_angles['beta'] = float(input("Input angle between the first and third translation vector in degrees: "))
        lengths_and_angles['gamma'] = float(input("Input angle between the second and third translation vector in degrees: "))
        element = input("Input the chemical symbol for the base atom: ")
        return Triclinic(symbol=element, latticeconstant=lengths_and_angles)

    if crystal_system == 2:
        print("Choose subsystem of the monoclinic cell\n",
              "1 - Simple monoclinic (also known as primitive monoclinic)\n",
              "2 - Base centered monoclinic")
        monoclinic_system = ("Input number of choosen monoclinic subsystem: ")
        while not (int(crystal_system) in range(1, 3)):
            cubic_system = input("Invalid input, enter an integer between 1 and 2: ")
        if monoclinic_system == 1:
            lengths_and_angles = {}
            lengths_and_angles['a'] = float(input("Input length of the first lattice translation vector in Å: "))
            lengths_and_angles['b'] = float(input("Input length of the second lattice translation vector in Å: "))
            lengths_and_angles['c'] = float(input("Input length of the third lattice translation vector in Å: "))
            print("First and second translation vectors are orthagonal")
            lengths_and_angles['beta'] = float(input("Input angle between the first and third translation vector: "))
            lengths_and_angles['gamma'] = float(input("Input angle between the second and third translation vector: "))
            element = input("Input one element in chemical notation: ")
            return SimpleMonoclinic(symbol=element, latticeconstant=lengths_and_angles)

        if monoclinic_system == 2:
            lengths_and_angles = {}
            lengths_and_angles['a'] = float(input("Input lengths of the first lattice translation vector: "))
            lengths_and_angles['b'] = float(input("Input lengths of the second lattice translation vector: "))
            lengths_and_angles['c'] = float(input("Input lengths of the third lattice translation vector: "))
            print("First and second translation vectors are orthagonal")
            lengths_and_angles['beta'] = float(input("Input angle between the first and third translation vector: "))
            lengths_and_angles['gamma'] = float(input("Input angle between the second and third translation vector: "))
            elements = input("Input the chemical symbols for the two base atoms separated by space: ")
            elements = elements.split()
            atoms = BaseCenteredMonoclinic(symbol=elements[0], latticeconstant=lengths_and_angles)
            atoms.set_chemical_symbols(elements)
            return atoms

    if crystal_system == 3:
        print("Choose subsystem of the orthorhombic cell\n",
              "1 - Simple orthorhomic (also known as primitive orthorhombic)\n",
              "2 - Base centered orthorhomic\n"
              "3 - Body centered orthorhomic\n",
              "4 - Face centered orthorhomic")
        cubic_system = input("Input number of choosen cubic subsystem: ")
        while (cubic_system in ['1', '2', '3', '4']):
            cubic_system = input("Invalid input, enter an integer between 1 and 4: ")
        cubic_system = int(cubic_system)
        if cubic_system == 1:
            lattice_constants = {}
            lattice_constants['a'] = float(input("Input x-direction lattice constant in Ångström: "))
            lattice_constants['b'] = float(input("Input y-direction lattice constant in Ångström: "))
            lattice_constants['c'] = float(input("Input z-direction lattice constant in Ångström: "))
            element = input("Input the chemical symbol for the base atom: ")
            return SimpleOrthorhombic(symbol=element, latticeconstant=lattice_constants)
        if cubic_system == 2:
            lattice_constants = {}
            lattice_constants['a'] = float(input("Input x-direction lattice constant in Ångström: "))
            lattice_constants['b'] = float(input("Input y-direction lattice constant in Ångström: "))
            lattice_constants['c'] = float(input("Input z-direction lattice constant in Ångström: "))
            elements = input("Input the chemical symbols for the two base atoms separated by space: ")
            elements = elements.split()
            atoms = BaseCenteredOrthorhombic(symbol=elements[0], latticeconstant=lattice_constants)
            atoms.set_chemical_symbols(elements)
            return atoms
        if cubic_system == 3:
            lattice_constants = {}
            lattice_constants['a'] = float(input("Input x-direction lattice constant in Ångström: "))
            lattice_constants['b'] = float(input("Input y-direction lattice constant in Ångström: "))
            lattice_constants['c'] = float(input("Input z-direction lattice constant in Ångström: "))
            elements = input("Input the chemical symbols for the two base atoms separated by space: ")
            elements = elements.split()
            atoms = BodyCenteredOrthorhombic(symbol=elements[0], latticeconstant=lattice_constants)
            atoms.set_chemical_symbols(elements)
            return atoms
        if cubic_system == 4:
            lattice_constants={}
            lattice_constants['a'] = float(input("Input x-direction lattice constant in Ångström: "))
            lattice_constants['b'] = float(input("Input y-direction lattice constant in Ångström: "))
            lattice_constants['c'] = float(input("Input z-direction lattice constant in Ångström: "))
            element = input("Input the chemical symbol for the base atom: ")
            atoms = FaceCenteredOrthorhombic(symbol=element, latticeconstant=lattice_constants)
            return atoms

    if crystal_system == 4:
        print("Choose subsystem of the tetragonal cell\n",
              "1 - Simple tetragonal\n",
              "2 - Centered tetragonal (also known as body centered tetragonal)")
        tetragonal_system = input("Input number of choosen tetragonal subsystem: ")
        while not (tetragonal_system in ['1', '2']):
            tetragonal_system = input("Invalid input, enter an integer between 1 and 2: ")
        tetragonal_system = int(tetragonal_system)
        if tetragonal_system == 1:
            lattice_constants = []
            lattice_constants.append(float(input("Input x-, y-direction lattice constant in Ångström: ")))
            lattice_constants.append(float(input("Input z-direction lattice constant in Ångström: ")))
            element = input("Input the chemical symbol for the base atom: ")
            return SimpleTetragonal(symbol=element, latticeconstant=lattice_constants)
        if tetragonal_system == 2:
            lattice_constants = []
            lattice_constants.append(float(input("Input x-, y-direction lattice constant in Ångström: ")))
            lattice_constants.append(float(input("Input z-direction lattice constant in Ångström: ")))
            elements = input("Input the chemical symbols for the two base atoms separated by space: ")
            elements = elements.split()
            atoms = CenteredTetragonal(symbol=elements[0], latticeconstant=lattice_constants)
            atoms.set_chemical_symbols(elements)
            return atoms

    if crystal_system == 7:
        print("Choose subsystem of the cubic cell \n",
              "1 - Simple Cubic\n",
              "2 - Body centered cubic\n",
              "3 - Face centered cubic")
        cubic_system = input("Input number of choosen cubic subsystem: ")
        while not (cubic_system in ['1', '2', '3']):
            cubic_system = input("Invalid input, must be integer between 1 and 3: ")
        cubic_system = int(cubic_system)
        if cubic_system == 1:
            lattice_constant = float(input("Input lattice constant in Ångström: "))
            element = input("Input one element in element notation: ")
            return SimpleCubic(symbol=element, latticeconstant=lattice_constant)
        if cubic_system == 2:
            lattice_constant = float(input("Input lattice constant in Ångström: "))
            elements = input("Input the chemical symbols for the two base atoms separated by space: ")
            elements = elements.split()
            atoms = BodyCenteredCubic(symbol=element, latticeconstant=lattice_constant)
            atoms.set_chemical_symbols(elements)
            return atoms
        if cubic_system == 3:
            lattice_constant = float(input("Input lattice constant in Ångström: "))
            elements = input("Input the chemical symbols for the four base atoms separated by space: ")
            elements = elements.split()
            atoms = FaceCenteredCubic(symbol=elements[0], latticeconstant=lattice_constant)
            atoms.set_chemical_symbols(elements)
            return atoms


if __name__ == "__main__":
    create_view_and_save_crystal_guided()
