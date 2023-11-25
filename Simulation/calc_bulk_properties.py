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
from elastic import get_pressure, BMEOS, get_strain
from elastic import get_elementary_deformations, scan_volumes
from elastic import get_BM_EOS, get_elastic_tensor
from elastic.elastic import get_cij_order
import ase.units as units


def calc_bulk_modulus(atoms: Atoms, output_dict={'bulk_modulus': []}):
    """Calculate the equilibrium bulk modulus B for solids using the equation of state.

    It is directly connnected to the second derivative given by the equation
    B = V * (d^2E/dV^2), bulk modulus (incompressibility constant) is a measure of
    substance's resistance to changes in volume when subject to compressive force from
    all directions.

    Args:
        atoms(ase atom object): The system to calculate the bulk modulus for
        traj_file: List of configuration files which provide us with volumes and potential energies

    Returns:
        (float): The optimal bulk modulus in eV/Ångstrom^3
    """
    cell = atoms.get_cell()
    volumes = []
    energies = []
    atoms_copy = atoms.copy()
    atoms_copy.calc = atoms.get_calculator()
    # Create for example only 10 differnet configuration for different lattice constant
    for element in np.linspace(0.999, 1.001, 10):
        atoms_copy.set_cell(cell * element, scale_atoms=True)
        volumes.append(atoms_copy.get_volume())
        energies.append(atoms_copy.get_potential_energy())
    # Equation of state
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    B = B / kJ * 1.0e24
    # print(B, "GPa")
    output_dict['bulk_modulus'].append(B)
    return B


def calculate_cohesive_energy(atoms, output_dict={'cohesive_energy': []}):
    """Calculate the cohesive energy of an Atoms object, molecule, cluster or Bulk).

    Keep in mind that the attached calculator should support the atoms object your are trying to create
    To calculate the cohesive energy of an Atoms object, crystal(cluster/Bulk), etc.
    we need first to find the energy required to separate its components
    into neutral free atoms at rest and at infinite separation,
    Formula: Cohesive energy=(isolated atoms potential energies-total atoms potential energy)/nr of atoms

    Args:
        atoms (ASE atoms object): The ASE atoms object for which to calculate the cohesive energy.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        (float): The Cohesive energy in eV.
    """
    # Check if the atoms object is a single atom to begin with
    if len(atoms) == 1:
        cohesive_energy = atoms.get_potential_energy()
        return cohesive_energy
    else:
        # Get the current calculator for the atoms objectmolecule
        original_calculator = atoms.get_calculator()

        # Get the potential energy for the entire atoms object as whole:
        total_atoms_potential_energy = atoms.get_potential_energy()
        # print("Total atoms potential energy:", total_atoms_potential_energy)

        # Loop for calculating the potential energy of each atom in the atoms object
        isolated_atoms_potential_energies = 0
        for i in range(len(atoms)):
            isolated_atom = Atoms([atoms[i]])
            # Set the  original calculator for the isolated atom
            isolated_atom.set_calculator(original_calculator)
            isolated_atom_potential_energy = isolated_atom.get_potential_energy()
            isolated_atoms_potential_energies += isolated_atom_potential_energy
        # print("isolated_atoms_potential_energies "+ str(isolated_atoms_potential_energies))

        # Calculate cohesive energy
        cohesive_energy = (isolated_atoms_potential_energies - total_atoms_potential_energy) / len(atoms)
        return cohesive_energy


def calc_elastic(atoms: Atoms, output_dict={'elastic_tensor': []}):
    """Calculate the elastic tensor C11.

    This is somewhat computationally heavy for large systems. Computational
    cost increases exponentially with the size of the system.

    Args:
        atoms(ase atoms object): Atoms object to calculate the tensor for.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        (float): Elastic tensor C11 in GPa.

    """
    systems = get_elementary_deformations(atoms, n=5, d=0.33)

    # Here we should run stress calculations on each deformed system

    Cij, Bij = get_elastic_tensor(atoms, systems)
    cij_order = get_cij_order(atoms)

    # After the program is restructured we can look at the following properties with the help of the elastic tensor:
    # (For the numbers to make sense, we might want to convert them to GPa by calculating 'number/units.GPa')
    #bulk_modulus = (Cij[cij_order.index('C_11')] + 2*Cij[cij_order.index('C_12')])/3
    #shear_modulus = (3*Cij[cij_order.index('C_44')]+Cij[cij_order.index('C_11')]-Cij[cij_order.index('C_12')])/5
    #youngs_modulus = 9*bulk_modulus*shear_modulus/(3*bulk_modulus+shear_modulus)
    #poisson_ratio = 1/2-youngs_modulus/(6*bulk_modulus)

    # We have compared the bulk modulus here with the one implemented the 'normal' way, and we got similar numbers at the start of 
    # the simulation, however when the energy freaks out in the NVE simulations, the bulk modulus here also freaks out, but the
    # one Issa implemented still looks good when it goes crazy. Or maybe not after further investigation.

    output_dict['elastic_tensor'].append(Cij[0]/units.GPa)
    return Cij/units.GPa