"""This script calculated the intesive and extensive properties."""
from ase import Atoms
import numpy as np
from configparser import ConfigParser
import os
import sys
from numpy import linalg as LA
from ase import units
import numpy as np
from scipy.spatial.distance import cdist
from heapq import nsmallest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')


def approx_lattice_constant(atoms):
    """Calculate the approximate lattice constant

    This is done by finding the distances to the four nearest neighbours for 
    every atom, summing the results and dividing by the amount of results

    Args:
        atoms(ase atom object): The system to calculate the lattice constant for.

    Returns:
        (float): The approximate lattice constant.
    """
    positions = np.array(atoms.get_positions())
    distances_between_atoms = cdist(positions, positions)

    lattice_contributions = 0
    for element in distances_between_atoms:
        # Since 0 is always present for each atom, we calculate the 5 nearest distances
        lattice_contributions += sum(nsmallest(5,element))
    
    return lattice_contributions/(4 * len(positions))


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
    """Calculate specific heat capacity of atoms object, Bear in mind, the system MUST reach equilibrium to get correct answers.

    Two formulas are used here, for NVE or NVT ensemble:
    1. For NVE: heat capacity = 3*N*kB / 2*(1-(2* (variance of the kinetic energy / (Nr of atoms)^2))/3*kB^2*T^2)
    2. For NVT: heat capacity = (variance of the total energy) /kB * T^2
    Finally the specific heat capacity = Heat capactiy / atoms mass

    Args:
        atoms(ase atom object): The system to calculate the specific heat capacity for.
        config_file(str): Name of the file with user's parameters.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        Specific heat capacity(float): The calculated specific heat capacity in Joule/Kilogram * Kelvin (J/Kg*K)
    """
    # Constants for unit conversion Physics handbook CU-2.4
    eV_to_Joules = 1.602177e-19  # Conversion factor from eV to Joules
    avogadro_number = 6.02214076e23  # atoms/mol

    # Read the user config data
    config_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    config_data = ConfigParser()
    config_data.read(config_path+config_file)

    # Skipping the first 80% of the unstable values 
    total_nr_values = len(output_dict["total_energy"])
    skippable_values = int(0.8 * total_nr_values)

    # Get the average of the instantaneous temperatures
    # Ignoring the first 80% of the unstable values at the beginning of the simulation
    average_temperature = np.mean(output_dict["temperature"][skippable_values:])

    if config_data['SimulationSettings']['ensemble'] == "NVE":
        # Ignoring the first 80% of the unstable values at the beginning of the simulation
        kinetic_energies = output_dict["kinetic_energy"][skippable_values:]
        kinetic_energies_array = np.array(kinetic_energies)

        # Calculate the variance of the kinetic energies
        var_kinetic_energies = np.var(kinetic_energies_array)

        # Formula from lecture 4 slide 49
        first_term = (3 * len(atoms) * units.kB) / 2
        # Added factor of 2 becuase we suspect the lecture's formula is wrong and it lacks factor of 2
        second_term = 2 * (1 / (1 - ((2 * var_kinetic_energies / len(atoms)**2) / (3 * (units.kB**2) * average_temperature**2))))
        heat_capacity = first_term * second_term

        # Unit change from eV to J
        heat_capacity_joule = heat_capacity * eV_to_Joules

        # Convert molar mass from g/mol to g/atom and then calculat this for the whole system
        molar_mass_atom = atoms.get_masses()[0] / avogadro_number
        total_molar_mass_atom = len(atoms) * molar_mass_atom

        # Formula from physics handbook F-2.3 c = C/m where c: specific heat capacity and C: heat capacity
        specific_heat_capacity = heat_capacity_joule / total_molar_mass_atom

        # Converting from J/K*g to J/K*Kg
        specific_heat_capacity = specific_heat_capacity * 1e3
        #print("specific_heat_capacity : ", specific_heat_capacity, "J/Kg*K")

    elif config_data['SimulationSettings']['ensemble'] == "NVT":
        # Ignoring the first 80% of the unstable values at the beginning of the simulation
        total_energies = output_dict["total_energy"][skippable_values:]
        total_energies_array = np.array(total_energies)

        # Calculate the variance of the total energies
        var_total_energies = np.var(total_energies_array)

        # Formula from lecture 4 slide 49
        heat_capacity = var_total_energies / (units.kB * average_temperature**2)

        # Unit change from eV to J
        heat_capacity_joule = heat_capacity * eV_to_Joules

        # Convert molar mass from g/mol to g/atom and then calculat this for the whole system
        molar_mass_atom = atoms.get_masses()[0] / avogadro_number
        total_molar_mass_atom = len(atoms) * molar_mass_atom

        # Formula from physics handbook F-2.3 c = C/m where c: specific heat capacity and C: heat capacity
        specific_heat_capacity = heat_capacity_joule / total_molar_mass_atom

        # Converting from J/K*g to J/K*Kg
        specific_heat_capacity = specific_heat_capacity * 1e3
        #print("specific_heat_capacity : ", specific_heat_capacity, "J/Kg*K")

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


def time_average_of_debye_temperature(atoms: Atoms, output_dict={'debye_temperature': []}):
    """Calculate the time average of debye temperature of an atoms object.

    The formula which is used: Debye temperature = (Planck constant * Debye frequency) / Boltzmann constant
    where Debye frequency = velocity of sound * ((6 * pi^2 * N) / volume)^(1/3)).

    Args:
        atoms(ase atom object): The system to calculate the time average of debye temperature for.
        output_dict(dict): Dictionary to append the result to.
    Returns:
        (float): The calculated time average of debye temperature
    """
    # From physics handbook CU-1.1
    hbar = 6.5821196e-16  # Planck constant in eV*s

    # Constant for unit conversion Physics handbook CU-2.4
    AtomicMass_to_Kg = 1.660539e-27     # Conversion factor from atomic

    volume_angstrom = atoms.get_volume()    # in Å^3
    volume = volume_angstrom * 1e-30    # in m^3
    num_atoms = len(atoms)

    # Density in Kg/m^3
    density = (sum(atoms.get_masses()) * AtomicMass_to_Kg) / volume

    # Calculating the velocity of sound in m/s
    # Velocity of sound values can be found in Physics handbook T-4.1
    velocity_of_sound = np.sqrt(((output_dict["bulk_modulus"][-1]) * 1e9) / density)  # in m/s

    # Calculate Debye frequency w_D, the formula can be found in three different places:
    # 1. In "Introduction to Solid State Physics" by Charles Kittel page 112
    # 2. Physics handbook F-10.4
    # 3.wiki: https://en.wikipedia.org/wiki/Debye_model#Debye_frequency
    # Debye temperature values can be found in Table 1 Page 116 in the same book above.
    w_D = velocity_of_sound * ((6 * np.pi**2 * num_atoms) / volume)**(1/3)

    # Calculate Debye temperature, formula from wiki debye temp = h_bar/kB * debye_frequency
    debye_temperature = (hbar * w_D) / units.kB

    # Append the calculated Debye temperature to the output dictionary
    output_dict['debye_temperature'].append(debye_temperature)

    # Calculate the time average, usless line of code
    # time_average_of_debye_temperature = np.mean(output_dict['debye_temperature'])
    return debye_temperature
