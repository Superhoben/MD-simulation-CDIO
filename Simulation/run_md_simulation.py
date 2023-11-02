"""This runs MD simulations."""
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
import configparser

def run_md_simulation(config_file_name: str, trajectory_file: str):

    """Skeleton for the MD simulation program.

    Currently it is written mostly in pseudo code.
    It is not supposed to work yet.

    Args:
        config_file_name(str): Name of file with parameters
        to use in simulation.

    Returns:
        atoms(ase atoms object): The ase atoms object after simulation.

    """

    # Parse the config file to get dictionary of data
    config_data = parse_config(config_file_name)

    # Create atoms object for simulation
    atoms = create_atoms_object(config_data)

    # Call appropriate function depending on simulation type
    simulation_type = config_data['ensemble']
    if simulation_type == "NVE" or simulation_type == "NVT":
        atoms = run_NVE_NVT(atoms, config_data, simulation_type)
    elif simulation_type == "something_else":
        # TODO: implement run_other_simulation, e.g.,:
        # atoms = run_other_simulation(atoms, config_data)
        raise Exception("Running calculations with 'run_other_simulation' " +
                        "is not implemented yet.")
    # Return atoms object after simulation
    return atoms


def run_NVE_NVT(atoms, config_data, simulation_type):
    """Run NVE or NVT simulation.

    Args:
        atoms(ase atoms object): The ase atoms object to simulate.
        config_data(dict): Dictionary of parameters for simulation.
        simulation_type(string): NVE or NVT

    Returns:
        atoms(ase atoms object): The ase atoms object after simulation.

    """
    # Set potential for simulation
    potential = config_data['potential']
    if potential == "EMT":
        atoms.calc = EMT()
    elif potential == "LennardJones":
        atoms.calc = LennardJones()
    elif potential == "other_potential":
        # TODO: implement running with other potentials, e.g.,:
        # atoms.calc = OtherPotential()
        raise Exception(
            "Running calculations with 'other_potential'" +
            "is not implemented yet."
        )

    # Set dynamics module depending on simulation type
    if simulation_type == "NVE":
        MaxwellBoltzmannDistribution(atoms,
                                     temperature_K=int(config_data['temperature']))
        dyn = VelocityVerlet(atoms, int(config_data['time_step'])*units.fs)
    elif simulation_type == "NVT":
        dyn = Langevin(
            atoms,
            timestep=config_data['time_step']*units.fs,
            temperature_K=config_data['temperature'],
            friction=(config_data['friction'] or 0.005),
        )
        MaxwellBoltzmannDistribution(atoms,
                                     temperature_K=config_data['temperature'])

    if config_data['show_properties'] is True:
        dyn.attach(show_properties(atoms, config_data),
                   interval=config_data['interval'])

    dyn.run(5000)

    return atoms


def show_properties(atoms, config_data):
    """Show properties in GUI duing simulations.

    Note that this is yet to be written.

    Args:
        atoms(ase atoms object): The ase atoms object to simulate.
        config_data(dict): Dictionary of parameters for simulation.
    """
    # TODO: implement print_in_gui and calc_temp:
    # print_in_gui(calc_temp())
    raise Exception("run_md_simulation.py: print_in_gui and calc_temp not " +
                    "implemented yet.")
    return


def parse_config(config_file_name):
    """Parse config file into a dictionary.

    Note that this is yet to be written.

    Args:
        config_file_name(str): Parameters to use
        in simulation.

    Returns:
        config_data(dict): Dictionary of parameters for simulation.

    """
    config = configparser.ConfigParser()
    config.read('../User_interface/' + config_file_name)
    config_data = config['config1']
    return config_data


def create_atoms_object(traj_file):
    """Create a list which contains our atom/atoms from the trajectory file.

    Args:
        traj_file: The trajectory file created from the material ID

    Returns:
        atom/atoms(list): The ase atom/atoms object.

    """
    traj = Trajectory(traj_file)
    atoms_list = []
    for atom in traj:
        atoms_list.append(atom)

    return atoms_list



if __name__ == "__main__":
    run_md_simulation("config1.ini")
