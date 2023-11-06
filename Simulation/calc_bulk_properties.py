"""The file will create the trajectory file with different configuration.

It also will read from the created trajectory file and calculate the wishes properties.
"""
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from ase.io import read
from ase.eos import EquationOfState
import numpy as np
from ase.units import kJ


def calc_bulk_modulus(atoms: Atoms, output_dict={'bulk_modulus': []}):
    """Calculate the equilibrium bulk modulus B for solids using the equation of state.

    It is directly connnected to the second derivative given by the equation
    B = V * (d^2E/dV^2), bulk modulus (incompressibility constant) is a measure of
    substance's resistance to changes in volume when subject to compressive force from
    all directions.

    Args:
        atoms(ase atom object): the system to calculate the bulk modulus for
        traj_file: list of configuration files which provide us with volumes and potential energies

    Returns:
        (float): the optimal bulk modulus in eV/Ångstrom^3
    """
    cell = atoms.get_cell()
    volumes = []
    energies = []
    # Create for example only 10 differnet configuration for different lattice constant
    for element in np.linspace(0.999, 1.001, 10):
        atoms.set_cell(cell * element, scale_atoms=True)
        volumes.append(atoms.get_volume())
        energies.append(atoms.get_potential_energy())
    # Equation of state
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    B = B / kJ * 1.0e24
    # print(B, "GPa")
    output_dict['bulk_modulus'].append(B)
    return B


def calculate_cohesive_energy(isolated_atoms, bulk_atoms):
    """Calculate the cohesive energy of an object, Atoms, cluster or Bulk.

    To calculate the cohesive energy of a crystal(cluster/Bulk), we need to find the energy required
    to separate its components into neutral free atoms at rest and at infinite separation,
    the formula: Cohesive energy = (energy of free atoms - crystal atoms energy) / nr of crystal (cluster/bulk) atoms

    Args:
        isolated_atoms (ase atoms object): the isolated atoms object
        bulk_atoms (ase atoms object): cluster or bulk atoms object

    Returns:
        (float): the Cohesive energy in eV
    """
    # Get all potential energy in a list for each atom exist in our object, molecule, soild etc.
    # Bear in mind, here the atoms are static so total energy = potentail energy only
    isolated_atoms_potential_energies = isolated_atoms.get_potential_energies()
    # Get the total potential energy for all our atoms
    total_isolated_atoms_potential_energy = 0
    for atom_potential_energy in isolated_atoms_potential_energies:
        total_isolated_atoms_potential_energy += atom_potential_energy
    bulk_atoms_potential_energy = bulk_atoms.get_potential_energy()
    bulk_atoms_kinetic_energy = bulk_atoms.get_kinetic_energy()
    bulk_total_energy = bulk_atoms_potential_energy + bulk_atoms_kinetic_energy
    num_atoms = len(bulk_atoms)
    cohesive_energy = (num_atoms*total_isolated_atoms_potential_energy - bulk_total_energy)/num_atoms
    return cohesive_energy

    # This will be used later probably so i am keeping this comment
    # This uses Issa Nseir's personal API-key to access the database
    # with MPRester("t4XwMQ3LLvLcnugLQQCCCII6BG85APG8") as mpr:
    # some_material = mpr.materials.search(material_ids=[material_id])
    # print(some_material)
