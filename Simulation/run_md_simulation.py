"""This runs MD simulations."""
import os
import sys
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT, LennardJones
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from configparser import ConfigParser
import json
from Simulation import calc_properties, calc_bulk_properties, lattice_constant
from parcalc import ParCalculate
from elastic import get_elastic_tensor
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')


def run_single_md_simulation(config_file: str, traj_file: str, output_name: str):
    """Run md simulation for a single trajectory file, with parameters specified in config

    Args:
        config_file(str): Name of file with parameters to use in simulation.
        traj_file(str): Name of trajectory file with the atoms object to use
                        in simulation
        output_name(str): Name of file to write results to

    Returns:
        atoms(ase atoms object): The ase atoms object after simulation.

    """
    # Parse the config file to get dictionary of data
    config_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    config_data = ConfigParser()
    config_data.read(config_path+config_file)

    # Create separate dicts for easier access of data
    simulation_settings = config_data['SimulationSettings']
    recording_intervals = config_data['RecordingIntervals']

    # Create atoms object for simulation
    traj_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
    traj = Trajectory(traj_path+traj_file)
    atoms = traj[0]

    # Set potential for simulation
    potential = simulation_settings['potential']
    if potential == "EMT":
        atoms.calc = EMT()
    elif potential == "LennardJones":
        atoms.calc = LennardJones()
    else:
        # TODO: implement running with other potentials, e.g.,:
        # atoms.calc = OtherPotential()
        raise Exception("Running calculations with potential '" + potential + "' is not implemented yet.")

    # Set dynamics module depending on simulation type
    ensemble = simulation_settings['ensemble']
    if ensemble == "NVE":
        MaxwellBoltzmannDistribution(atoms, temperature_K=int(simulation_settings['temperature']))
        dyn = VelocityVerlet(atoms, int(simulation_settings['time_step'])*units.fs)
    elif ensemble == "NVT":
        dyn = Langevin(atoms, timestep=int(simulation_settings['time_step'])*units.fs,
                       temperature_K=int(simulation_settings['temperature']),
                       friction=(float(simulation_settings['friction']) or 0.005))
        MaxwellBoltzmannDistribution(atoms, temperature_K=int(simulation_settings['temperature']))
    else:
        # TODO: implement other ensembles
        raise Exception("Running calculations with ensemble '" + ensemble + "' is not implemented yet.")

    # This dict will contain output data of the simulation to be written into the output text file
    output_dict = {}
    output_dict["config_file"] = config_file
    # Attach recorders that calculate a certain property and store in the output_dict
    interval_to_record_energy = int(recording_intervals['record_energy'])
    if interval_to_record_energy:
        output_dict['total_energy'] = []
        output_dict['kinetic_energy'] = []
        output_dict['potential_energy'] = []
        dyn.attach(calc_properties.calc_energy, interval_to_record_energy, atoms, output_dict)

    interval_to_record_temperature = int(recording_intervals['record_temperature'])
    if interval_to_record_temperature:
        output_dict['temperature'] = []
        dyn.attach(calc_properties.calc_temp, interval_to_record_temperature, atoms, output_dict)

    interval_to_record_pressure = int(recording_intervals['record_pressure'])
    if interval_to_record_pressure:
        output_dict['pressure'] = []
        dyn.attach(calc_properties.calc_pressure, interval_to_record_pressure, atoms, output_dict)

    interval_to_record_bulk_modulus = int(recording_intervals['record_bulk_modulus'])
    if interval_to_record_bulk_modulus:
        output_dict['bulk_modulus'] = []
        dyn.attach(calc_bulk_properties.calc_bulk_modulus, interval_to_record_bulk_modulus, atoms, output_dict)

    interval_to_record_optimal_scaling = int(recording_intervals['record_optimal_scaling'])
    if interval_to_record_optimal_scaling:
        output_dict['optimal_scaling'] = []
        output_dict['iterations_to_find_scaling'] = []
        dyn.attach(lattice_constant.optimize_scaling, interval_to_record_optimal_scaling, atoms, output_dict)

    interval_to_record_configuration = int(config_data['RecordingIntervals']['record_configuration'])
    if interval_to_record_configuration:
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files/'
        traj = Trajectory(path+output_name+'.traj', "w", atoms)
        dyn.attach(traj.write, interval_to_record_configuration)

    interval_to_record_elastic_properties = int(recording_intervals['record_elastic'])
    if interval_to_record_elastic_properties:
        output_dict['elastic_tensor'] = []
        dyn.attach(calc_bulk_properties.calc_elastic, interval_to_record_elastic_properties, atoms, output_dict)

    # Run simulation with the attached recorders
    dyn.run(int(simulation_settings['step_number']))

    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/'
    with open(path + output_name + '.txt', 'w') as file:
        file.seek(0)
        json.dump(output_dict, file)

    return atoms


def run_md_simulation(config_file_list, trajectory_file_list):
    """Run md simulations for multiple configs and trajectory files.

    Args:
        config_file_list(list[str]): List of names of files with parameters to
                                     use in the simulations.
        traj_file(str): List of name of trajectory files with the atoms object
                        to use in the simulations

    Returns:
        None

    """
    i = 0
    for traj_file in trajectory_file_list:
        i = i+1
        run_single_md_simulation(config_file_list[0], traj_file, 'output'+str(i))


if __name__ == "__main__":
    run_single_md_simulation("example_config.ini", 'mp-30.traj', 'out_put1')
