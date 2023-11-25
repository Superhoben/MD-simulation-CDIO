from ase import Atoms
import numpy as np
from configparser import ConfigParser
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')


def calc_temp(atoms: Atoms, output_dict={'temperature': []}):
    """ Calculates temperature of atoms object

    Args:
        atoms(ase atom object): the system to calculate the temperature for

    Returns:
        (float): the calculated temperature
    """

    temperature = atoms.get_temperature()
    output_dict['temperature'].append(temperature)
    return temperature


def calc_energy(atoms: Atoms, output_dict={'total_energy': [], 'kinetic_energy': [], 'potential_energy': []}):
    """ Calculates total, kinetic, and potential energy of atoms object

    Args:
        atoms(ase atom object): the system to calculate the energy for
        output_dict(dict): dictionary to append the result to

    Returns:
        (float): the calculated total energy
    """
    total_energy = atoms.get_total_energy()
    output_dict['total_energy'].append(total_energy)
    output_dict['kinetic_energy'].append(atoms.get_kinetic_energy())
    output_dict['potential_energy'].append(atoms.get_potential_energy())
    return total_energy


def calc_pressure(atoms: Atoms, output_dict={'pressure': []}, external_field=None):
    """Calculate pressure of atoms object with or without an external field.

    The formula used is P=1/3V*(2*E_kin(t)+sum_over_all_atoms{r_i*f_i}) where
    r_i and f_i is the position of and force on atom i and V referes to the
    volume of the unitcell.

    Args:
        atoms(ase atom object): the system to calculate the pressure for
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]

    Returns:
        (float): the calculated pressure in GPa (giga pascal)
    """
    forces = np.array(atoms.get_forces(apply_constraint=False, md=True))
    positions = np.array(atoms.get_positions())

    volume = atoms.get_volume()        # in Å^3
    ekin = atoms.get_kinetic_energy()  # in eV

    if external_field == None:
        pressure_in_eV_per_Å3 = (2*ekin+np.sum(np.multiply(forces, positions)))/(3*volume)
        output_dict['pressure'].append(pressure_in_eV_per_Å3*160.21766208)
        return pressure_in_eV_per_Å3*160.21766208
    else:
        forces += external_field(atoms)
        # When having an external field point of origin needs to be centered
        centered_positions = positions-positions.mean(axis=0)
        pressure_in_eV_per_Å3 = (2*ekin+np.sum(np.multiply(forces, centered_positions)))/(3*volume)
        output_dict['pressure'].append(pressure_in_eV_per_Å3*160.21766208)
        return pressure_in_eV_per_Å3*160.21766208


def calculate_specific_heat(atoms, config_file, output_dict):
    """Calculate specific heat capacity of atoms object.
    The two formulas used are from lecture 4 slide 49.

    Args:
        atoms(ase atom object): the system to calculate the specific heat capacity for.
        config_file(str): Name of the file with user's parameters.
        output_dict(dict): dictionary to append the result to.


    Returns:
        (float): the calculated specific heat capacity in Joule/Kilogram * Kelvin (J/Kg*K)
    """

    # Define the Boltzmann constant (kB)
    kB = 1.380649e-23  # Boltzmann constant in J/K

    # Constants for unit conversion Physics handbook CU-2.4
    eV_to_Joules = 1.602177e-19  # Conversion factor from eV to Joules
    AtomicMass_to_Kg = 1.660539e-27  # Conversion factor from atomic mass units (u) to kilograms

    # Read the user config data
    config_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    config_data = ConfigParser()
    config_data.read(config_path+config_file)

    # Get the initial temperature
    initial_temperature = float(config_data['SimulationSettings']['temperature'])

    if config_data['SimulationSettings']['ensemble'] == "NVE":
        kinetic_energies = output_dict["kinetic_energy"]
        kinetic_energies_array = np.array(kinetic_energies)
        # Unit change from eV to J
        kinetic_energies_array_joule = kinetic_energies_array * eV_to_Joules
        # Calculate the variance of the kinetic energies
        var_kinetic_energies = np.var(kinetic_energies_array_joule)

        print("var of kinetic energy " + str(var_kinetic_energies))

        # Formula from lecture 4 slide 49
        first_term = (3 * len(atoms) * kB) / 2
        second_term = 1 / (1 - (2 * var_kinetic_energies / (3 * kB**2 * initial_temperature**2)))
        heat_capacity = first_term * second_term

        print("Heat Capacity: ", heat_capacity)
        # Formula from physics handbook F-2.3 c = C/m where c: specific heat capacity and C: heat capacity
        specific_heat_capacity = heat_capacity / (sum(atoms.get_masses()) * AtomicMass_to_Kg)
        print("sum of atomic mass: "+str(sum(atoms.get_masses())))
        print("specific_heat_capacity: " + str(specific_heat_capacity))

    elif config_data['SimulationSettings']['ensemble'] == "NVT":
        total_energies = output_dict["total_energy"]
        total_energies_array = np.array(total_energies)
        # Unit change from eV to J
        total_energies_array_joule = total_energies_array * eV_to_Joules
        # Calculate the variance of the total energies
        var_total_energies = np.var(total_energies_array_joule)

        # Formula from lecture 4 slide 49
        heat_capacity = var_total_energies / (kB * initial_temperature**2)

        print("Heat Capacity:", heat_capacity)
        # Formula from physics handbook F-2.3 c = C/m where c: specific heat capacity and C: heat capacity
        specific_heat_capacity = heat_capacity / sum(atoms.get_masses() * AtomicMass_to_Kg)
        print("Specific Heat Capacity:", specific_heat_capacity)

    return specific_heat_capacity
