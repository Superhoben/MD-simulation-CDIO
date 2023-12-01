"""This script calculated the intesive and extensive properties."""
from ase import Atoms
import numpy as np
from configparser import ConfigParser
import os
import sys
from numpy import linalg as LA
from ase import units
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')


def calc_temp(atoms: Atoms, output_dict={'temperature': []}):
    """Calculate the temperature of atoms object.

    Args:
        atoms(ase atom object): The system to calculate the temperature for.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        (float): The calculated temperature.
    """
    temperature = atoms.get_temperature()
    output_dict['temperature'].append(temperature)
    return temperature


def calc_energy(atoms: Atoms, output_dict={'total_energy': [], 'kinetic_energy': [], 'potential_energy': []}):
    """Calculate the total, kinetic, and potential energy of atoms object.

    Args:
        atoms(ase atom object): The system to calculate the energy for.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        (float): The calculated total energy.
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
        atoms(ase atom object): The system to calculate the pressure for.
        output_dict(dict): Dictionary to append the result to.
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]

    Returns:
        (float): The calculated pressure in GPa (giga pascal).
    """
    forces = np.array(atoms.get_forces(apply_constraint=False, md=True))
    positions = np.array(atoms.get_positions())

    volume = atoms.get_volume()        # in Å^3
    ekin = atoms.get_kinetic_energy()  # in eV

    if external_field is None:
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

    Two formulas are used here, for NVE or NVT ensemble:
    1. For NVE: heat capacity = 3*N*kB / 2*(1-(2*variance of the kinetic energy/Number of atoms^2)/3*kB^2*T^2)
    2. For NVT: heat capacity = (variance of the total energy/Number of atoms^2) /kB * T^2
    Finally the specific heat capacity = Heat capactiy / sum of the atoms masses

    Args:
        atoms(ase atom object): The system to calculate the specific heat capacity for.
        config_file(str): Name of the file with user's parameters.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        Specific heat capacity(float): The calculated specific heat capacity in Joule/Kilogram * Kelvin (J/Kg*K)
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

    # Get the average of the instantaneous temperatures
    # Ignoring the first 10 unstable values at the begining of the simulation
    average_temperature = np.mean(output_dict["temperature"][10:])

    if config_data['SimulationSettings']['ensemble'] == "NVE":
        # Ignoring the first 10 unstable values at the begining of the simulation
        kinetic_energies = output_dict["kinetic_energy"][10:]
        kinetic_energies_array = np.array(kinetic_energies)
        # Unit change from eV to J
        kinetic_energies_array_joule = kinetic_energies_array * eV_to_Joules
        # Calculate the variance of the kinetic energies
        var_kinetic_energies = np.var(kinetic_energies_array_joule)

        # Formula from lecture 4 slide 49
        first_term = (3 * len(atoms) * kB) / 2
        second_term = 1 / (1 - ((2 * var_kinetic_energies / len(atoms)**2) / (3 * kB**2 * average_temperature**2))) * 2
        heat_capacity = first_term * second_term

        # print("Heat Capacity: ", heat_capacity)
        # Formula from physics handbook F-2.3 c = C/m where c: specific heat capacity and C: heat capacity
        specific_heat_capacity = heat_capacity / (sum(atoms.get_masses()) * AtomicMass_to_Kg)
        # print("specific_heat_capacity: " + str(specific_heat_capacity))

    elif config_data['SimulationSettings']['ensemble'] == "NVT":
        # Ignoring the first 10 unstable values at the begining of the simulation
        total_energies = output_dict["total_energy"][10:]
        total_energies_array = np.array(total_energies)
        # Unit change from eV to J
        total_energies_array_joule = total_energies_array * eV_to_Joules
        # Calculate the variance of the total energies
        var_total_energies = np.var(total_energies_array_joule)

        # Formula from lecture 4 slide 49
        heat_capacity = (var_total_energies / len(atoms)**2) / (kB * average_temperature**2)
        # Formula from physics handbook F-2.3 c = C/m where c: specific heat capacity and C: heat capacity
        specific_heat_capacity = heat_capacity / (sum(atoms.get_masses()) * AtomicMass_to_Kg)

    return specific_heat_capacity


def calc_mean_square_displacement(atoms: Atoms, output_dict={'mean_square_displacement': []}):
    """Calculate the mean square displacement of atoms object.

    The formula used is MSD=1/N*sum{(r_i(t_n)-r_i(t_0))^2 where r_i is the
    position of atom i at time t_n and N is the number of atoms

    Args:
        atoms(ase atom object): the system to calculate the mean square displacement for
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]

    Returns:
        (float): the calculated mean square displacement
    """
    positions = np.array(atoms.get_positions())
    MSD = 0

    if output_dict['mean_square_displacement'] == []:
        # Append initial position array which is used for all future iterations
        output_dict['mean_square_displacement'].append(positions)
    else:
        atom_pos_diffs = positions-output_dict['mean_square_displacement'][0]
        MSD_sum = 0
        for atom in atom_pos_diffs:
            MSD_sum += LA.norm(atom)**2
        MSD = MSD_sum/len(positions)
        output_dict['mean_square_displacement'].append(MSD)

    return MSD


def lindemann_criterion(atoms: Atoms, output_dict={'lindemann_criterion': []}, d = 1):
    """Calculate the Lindemann criterion of atoms object.

    The formula used is L = 1/d*(MSD)^(1/2).

    Args:
        atoms(ase atom object): the system to calculate the lindemann criterion for
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]
        d(int): the nearest neighbour distance between the atoms

    Returns:
        (float): the calculated Lindemann criterion
    """
    lindemann = 0

    if output_dict['lindemann_criterion'] != []:
        lindemann = np.sqrt(output_dict['mean_square_displacement'][-1])/d

    output_dict['lindemann_criterion'].append(lindemann)

    return lindemann


def self_diffusion_coefficent(atoms: Atoms, output_dict={'lindemann_criterion': []}, time_elapsed_per_interval = 1):
    """Calculate the self-diffusion coefficient of atoms object.

    The formula used is D = 1/(6*t)*MSD where t is time elapsed at a certian iteration

    Args:
        atoms(ase atom object): the system to calculate the self-diffusion coefficient for
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]
        time_elapsed_per_interval(s): number of seconds that have passed per iteration

    Returns:
        (float): the calculated self-diffusion coefficient
    """
    self_diffusion_coefficient = 0

    if output_dict['self_diffusion_coefficient'] != []:
        # t is calculated by using how many calculations of calc_self_diffusion_coefficient have
        # been performed, and multiplying this by how long the intervals are.
        self_diffusion_coefficient = output_dict['mean_square_displacement'][-1]/(6*time_elapsed_per_interval*(len(output_dict['mean_square_displacement']) - 1)*units.fs)

    output_dict['self_diffusion_coefficient'].append(self_diffusion_coefficient)

    return self_diffusion_coefficient
