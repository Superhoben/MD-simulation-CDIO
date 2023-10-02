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
    for element in np.linspace(0.95, 1.05, 10):
        atoms.set_cell(cell * element, scale_atoms=True)
        atoms.get_potential_energy()
        traj.write(atoms)


def calc_bulk_modulus(traj_file):
    """Calculate the equilibrium bulk modulus B for solids using the equation of state.

    it is directly connnected to the second derivative given by the equation
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
    print(B, "GPa")
    return B


if __name__ == "__main__":
    # Testing the code
    lattice_constant = 4
    approximate_pos = lattice_constant/2
    atoms = Atoms('Ni',
                  cell=[(0, approximate_pos, approximate_pos),
                        (approximate_pos, 0, approximate_pos),
                        (approximate_pos, approximate_pos, 0)],
                  pbc=1,
                  calculator=EMT())
    create_traj_file(atoms, lattice_constant)
    # Read the traj file from the first atom to the 10th atom
    calc_bulk_modulus("atoms.traj@0:9")
