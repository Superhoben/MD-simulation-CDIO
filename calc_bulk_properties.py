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


def create_traj_file(atoms: Atoms, lattice_constant: float):
    """Create our trajectory file from the best given lattice constant.

    Args:
        atoms(ase atom object): the system to create the trajectory file for
        lattice constant (float): best lattice constant

    Returns:
       None
    """
    # approximate_pos = lattice_constant /2
    cell = atoms.get_cell()
    traj = Trajectory("atoms.traj", "w")   # Creating the traj. file
    # Create for example only 10 differnet configuration for different lattice constant
    for element in np.linspace(0.999, 1.001, 10):
        atoms.set_cell(cell * element, scale_atoms=True)
        atoms.get_potential_energy()
        traj.write(atoms)


def calc_bulk_modulus(traj_file):
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
    traj_file = read(traj_file)
    # Extract volumes and energies
    volumes = []
    energies = []
    for element in traj_file:
        volume = element.get_volume()   # in Å^3
        volumes.append(volume)
        energy = element.get_potential_energy()  # in eV
        energies.append(energy)
    # Equation of state
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    B = B / kJ * 1.0e24
    # print(B, "GPa")
    return B


def calculate_cohesive_energy(atoms, bulk):
    """Calculate the cohesive energy of an object, Atoms,cluster or Bulk.

    To calculate the cohesive energy of a crystal(cluster/Bulk), we need to find the energy required
    to separate its components into neutral free atoms at rest and at infinite separation, the formula is
    Cohesive energy = (energy of free atoms - crystal atoms energy) / nr of crystal (cluster/bulk) atoms

    Args:
        atoms(ase atom object): the isolated atom object
        bulk:  cluster or bulk object

    Returns:
        (float): the optimal Cohesive energy in eV
    """
    # Assigning better names for the isolated atom object and the bulk object
    isolated_atom = atoms
    isolated_atom_energy = isolated_atom.get_potential_energy()
    bulk_atoms = bulk
    bulk_atoms_energy = bulk_atoms.get_potential_energy()
    num_atoms = len(bulk_atoms)
    cohesive_energy = (bulk_atoms_energy - isolated_atom_energy)/num_atoms
    return cohesive_energy

    # This will be used later probably so i am keeping this comment
    # This uses Issa Nseir's personal API-key to access the database
    # with MPRester("t4XwMQ3LLvLcnugLQQCCCII6BG85APG8") as mpr:
    # some_material = mpr.materials.search(material_ids=[material_id])
    # print(some_material)
