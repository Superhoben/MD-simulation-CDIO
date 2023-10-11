"""This runs md simulations."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from ase.calculators.lj import LennardJones
from asap3.md.langevin import Langevin


def run_md_simulation(config):
    """Skeleton for the MD simulation program.

    Currently it is written mostly in pseudo code.
    It is not supposed to work yet.

    Args:
        config(either the file path or the whole file): Parameters
        to use in simulation.

    Returns:
        atoms(ase atoms object): The ase atoms object after simulation.

    """
    # Parse the config file to get dictionary of data
    config_data = parse_config(config)

    # Create atoms object for simulation
    atoms = create_atoms_object(config_data)

    # Set potential for simulation
    if config_data.potential == 'EMT':
        atoms.calc = EMT()
    elif config_data.potential == 'LennardJones':
        atoms.calc = LennardJones()
    elif config_data.potential == 'other_potential':
        # TODO: implement running with other potentials, e.g.,:
        # atoms.calc = OtherPotential()
        raise Exception("Running calculations with 'other_potential'
                        is not implemented yet.")

    # Call appropriate function depending on simulation type
    simulation_type = config_data.simulation_type
    if simulation_type == 'NVE' or simulation_type == 'NVT':
        atoms = run_NVE_NVT(atoms, config_data, simulation_type)
    elif simulation_type == 'something_else':
        # TODO: implement run_other_simulation, e.g.,:
        # atoms = run_other_simulation(atoms, config_data)
        raise Exception("Running calculations with 'run_other_simulation'
                        is not implemented yet.")

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
    MaxwellBoltzmannDistribution(atoms,
                                 temperature_K=config_data.temperature)

    if simulation_type == "NVE":
        dyn = VelocityVerlet(atoms, config_data.time_step)
        # time_step could be entered in the unit fs instead, like this:
        # dyn = VelocityVerlet(atoms, config_data.time_step * units.fs)
    elif simulation_type == "NVT":
        dyn = Langevin(atoms, config_data.time_step,
                       temperature_K=config_data.temperature,
                       friction=(config_data.friction or 0.005))

    dyn.attach(show_properties(atoms, config_data),
               interval=config_data.interval)
    dyn.run(config_data.iterations)

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
    raise Exception("run_md_simulation.py: print_in_gui and calc_temp not
                    implemented yet.")
    return


def parse_config(config):
    """Parse config file into a dictionary.

    Note that this is yet to be written.

    Args:
        config(either the file path or the whole file): Parameters to use
        in simulation.

    Returns:
        config_data(dict): Dictionary of parameters for simulation.

    """
    # TODO: implement parse_config to actually parse config data
    # return config_data
    raise Exception("run_md_simulation.py: parse_config not implemented yet.")


def create_atoms_object(config_data):
    """Create atoms object from the config data.

    Note that this is yet to be written.

    Args:
        config_data(dict): Dictionary of parameters for simulation.

    Returns:
        atoms(ase atoms object): The ase atoms object to simulate.

    """
    # TODO: implement parse_config to actually parse config data
    # return atoms
    raise Exception("run_md_simulation.py: create_atoms_object not
                    implemented yet.")
