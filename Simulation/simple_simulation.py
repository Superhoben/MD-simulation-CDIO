"""This is intended for relaxing structures during property calculations."""
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3.md.verlet import VelocityVerlet
from ase import units
from ase import Atoms
from asap3 import EMT, LennardJones
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import Stationary
import numpy as np
from scipy.spatial.distance import cdist


def run_simple_md_simulation(atoms: Atoms, simulation_settings: dict, time_average_energy=False):
    """Run md simulation for an atoms object, with parameters from a dict.

    Args:
        atoms(ase atom object): The system to bully with simulation.
        simulation_settings(dict): Dictionary of simulation settings.

    Returns:
        atoms(ase atoms object): The ase atoms object after simulation.

    """
    # Set potential for simulation
    potential = simulation_settings['potential']
    if potential == "EMT":
        # EMT needs parameter values as input arguments for all materials except
        # Ni, Cu, Pd, Ag, Pt and Au. Note that we don't know if it works with mixed
        # materials or if it will make incorrect assumptions
        atoms.calc = EMT()
    elif potential == "LennardJones":
        # Lennard Jones is generally valid for gases and liquid but rarely solids
        # and not metals as far as I understand it //Gustav
        atoms.calc = LennardJones()
    else:
        # TODO: implement running with other potentials, e.g.,:
        # atoms.calc = OtherPotential()
        raise Exception("Running calculations with potential '" + potential + "' is not implemented yet.")

    # Set dynamics module depending on simulation type
    ensemble = simulation_settings['ensemble']
    if ensemble == "NVE":
        MaxwellBoltzmannDistribution(atoms, temperature_K=2*int(simulation_settings['temperature']))
        Stationary(atoms)
        dyn = VelocityVerlet(atoms, int(simulation_settings['time_step'])*units.fs)
    elif ensemble == "NVT":
        dyn = Langevin(atoms, timestep=int(simulation_settings['time_step'])*units.fs,
                       temperature_K=int(simulation_settings['temperature']),
                       friction=(float(simulation_settings['friction']) or 0.005))
        MaxwellBoltzmannDistribution(atoms, temperature_K=2*int(simulation_settings['temperature']))
        Stationary(atoms)
    else:
        # TODO: implement other ensembles
        raise Exception("Running calculations with ensemble '" + ensemble + "' is not implemented yet.")

    if time_average_energy:
        energies = {}
        energies['total_energy'] = []
        dyn.attach(calc_total_energy, 5, atoms, energies)

    dyn.run(int(simulation_settings['step_number']))

    if time_average_energy:
        return atoms, np.average(energies['total_energy'])/len(atoms)

    return atoms


def calc_total_energy(atoms: Atoms, output_dict={'total_energy': []}):
    """Calculate the total energy of atoms object.

    Args:
        atoms(ase atom object): The system to calculate the energy for.
        output_dict(dict): Dictionary to append the result to.

    Returns:
        (float): The calculated total energy.
    """
    total_energy = atoms.get_total_energy()
    output_dict['total_energy'].append(total_energy)
    return total_energy
