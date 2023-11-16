from ase import Atoms
import numpy as np
from numpy import linalg as LA
from scipy.spatial.distance import cdist


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
    

def calc_mean_square_displacement(atoms: Atoms, output_dict={'mean_square_displacement': []}, external_field=None):
    """Calculate pressure of atoms object with or without an external field.

    The formula used is MSQ=1/N*sum{(r_i(t_n)*r_i(t_0))^2 where r_i is the 
    position of atom i

    Args:
        atoms(ase atom object): the system to calculate the pressure for
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]

    Returns:
        (float): the calculated pressure in GPa (giga pascal)
    """
    positions = np.array(atoms.get_positions()) 
    MSD = 0

    if output_dict['mean_square_displacement'] == []:
        output_dict['mean_square_displacement'].append(positions)
    else:
        atom_pos_diffs = positions-output_dict['mean_square_displacement'][0]
        MSD_sum = 0
        for atom in atom_pos_diffs:   
            MSD_sum += LA.norm(atom)**2
        MSD = MSD_sum/len(positions)
        output_dict['mean_square_displacement'].append(MSD)

    return MSD


def lindemann_criterion(atoms: Atoms, output_dict={'lindemann_criterion': []}, external_field=None):
    """Calculate pressure of atoms object with or without an external field.

    The formula used is MSQ=1/N*sum{(r_i(t_n)*r_i(t_0))^2 where r_i is the 
    position of atom i

    Args:
        atoms(ase atom object): the system to calculate the pressure for
        external_field(function(Atoms)->np.array): A function which takes an ase Atom object
            and returns the force on each atom as an array in the same format as Atoms.force
            function would but converted to an np.array. For N atoms in 3 dimensions:
            [[f_x_atom1, f_y_atom1, f_z_atom1], ..., [f_x_atomN, f_y_atomN, f_z_atomN]]

    Returns:
        (float): the calculated pressure in GPa (giga pascal)
    """
    positions = np.array(atoms.get_positions())
    lindemann = 0
    
    if output_dict['lindemann_criterion'] == []:
        output_dict['lindemann_criterion'].append(lindemann)
    else:
        MSD = output_dict['mean_square_displacement'][-1]
        distance_between_atoms = cdist(positions, positions)
        i = 0
        while i < len(positions):
            distance_between_atoms

        print(cdist(positions, positions))

        #lindemann = MSD/np.min(cdist(positions, positions))
        
        """atom_pos_diffs = positions-output_dict['mean_square_displacement'][0]
        MSD_sum = 0
        for atom in atom_pos_diffs:
            MSD_sum += LA.norm(atom)
        MSD = MSD_sum/len(positions)"""


        output_dict['lindemann_criterion'].append(lindemann)

    return lindemann