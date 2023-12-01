"""This runs MD simulations."""
import os
import sys
import shutil
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from ase import Atoms
from asap3 import EMT, LennardJones
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from configparser import ConfigParser
import json
from Simulation import calc_properties, calc_bulk_properties, lattice_constant
from ase.md.velocitydistribution import Stationary
import numpy as np
from scipy.spatial.distance import cdist
from multiprocessing import Process
from ase.lattice.cubic import FaceCenteredCubic

def print_and_increase_progress(progress, sim_number):
    if sim_number:
        print("Simulation ", sim_number, ": Progress "+str(progress[0])+"%")
    print("Progress "+str(progress[0])+"%")
    progress[0] += 10

def run_single_md_simulation(config_file: str, traj_file: str, output_name: str, sim_number=0):
    """Run md simulation for a single trajectory file, with parameters specified in config.

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

    # This dict will contain output data of the simulation to be written into the output text file

    # Approximate lattice constant
    d = calc_properties.approx_lattice_constant(atoms)

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

    interval_to_record_configuration = int(recording_intervals['record_configuration'])
    if interval_to_record_configuration:
        path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files/'
        traj = Trajectory(path+output_name+'.traj', "w", atoms)
        dyn.attach(traj.write, interval_to_record_configuration)

    interval_to_record_elastic_properties = int(recording_intervals['record_elastic'])
    if interval_to_record_elastic_properties:
        output_dict['elastic_tensor_c11'] = []
        output_dict['bulk_modulus_from_tensor'] = []
        output_dict['shear_modulus'] = []
        output_dict['youngs_modulus'] = []
        output_dict['poisson_ratio'] = []
        dyn.attach(calc_bulk_properties.calc_elastic, interval_to_record_elastic_properties, atoms, output_dict)

    interval_to_record_mean_square_displacement = int(recording_intervals['record_mean_square_displacement'])
    if interval_to_record_mean_square_displacement:
        output_dict['mean_square_displacement'] = []
        dyn.attach(calc_properties.calc_mean_square_displacement, interval_to_record_mean_square_displacement, atoms, output_dict)

    interval_to_record_lindemann_criterion = int(recording_intervals['record_lindemann_criterion'])
    if interval_to_record_lindemann_criterion:
        output_dict['lindemann_criterion'] = []
        dyn.attach(calc_properties.lindemann_criterion, interval_to_record_lindemann_criterion, atoms, output_dict, d)

    interval_to_record_self_diffusion_coefficient = int(recording_intervals['record_self_diffusion_coefficient'])
    if interval_to_record_lindemann_criterion:
        output_dict['self_diffusion_coefficient'] = []
        dyn.attach(calc_properties.self_diffusion_coefficent, interval_to_record_self_diffusion_coefficient, atoms, output_dict, interval_to_record_self_diffusion_coefficient*int(simulation_settings['time_step']))

    progress = [0]
    ten_percent_interval = int(0.1 * float(simulation_settings['step_number']))
    dyn.attach(print_and_increase_progress, ten_percent_interval, progress, sim_number)

    # Run simulation with the attached recorders
    if sim_number:
        print("Simulation ", sim_number, " started")
    else:
        print("Simulation started")

    dyn.run(int(simulation_settings['step_number']))

    if interval_to_record_energy:
        output_dict['specific_heat_capacity'] = [calc_properties.calculate_specific_heat(atoms, config_file, output_dict)]

    if int(recording_intervals['record_mean_square_displacement']):
        output_dict['mean_square_displacement'][0] = 0
        time_avg_MSD = 0
        for MSD in output_dict['mean_square_displacement']:
            time_avg_MSD = time_avg_MSD + MSD
        output_dict["avg_MSD"] = time_avg_MSD/int(simulation_settings['step_number'])

    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/'
    with open(path + output_name + '.txt', 'w') as file:
        file.seek(0)
        json.dump(output_dict, file)

    if sim_number:
        print("Simulation ", sim_number, " finished!")
    else:
        print("Simulation finished!")
    return atoms


def run_md_simulations(config_file_name, trajectory_file_dir, output_dir_name):
    """Run md simulations for multiple configs and trajectory files.

    Args:
        config_file_list(list[str]): List of names of files with parameters to
                                     use in the simulations.
        traj_file(str): List of name of trajectory files with the atoms object
                        to use in the simulations

    Returns:
        None
    """
    base_path = os.path.dirname(os.path.abspath(__file__)) + '/../'
    path_new_output_txt_dir = base_path + 'Output_text_files/' + output_dir_name
    path_new_output_traj_dir = base_path + 'Output_trajectory_files/' + output_dir_name
    if os.path.exists(path_new_output_txt_dir):
        shutil.rmtree(path_new_output_txt_dir)
    if os.path.exists(path_new_output_traj_dir):
        shutil.rmtree(path_new_output_traj_dir)
    os.mkdir(base_path + 'Output_text_files/' + output_dir_name)
    os.mkdir(base_path + 'Output_trajectory_files/' + output_dir_name)
    traj_files_full_path = base_path + 'Input_trajectory_files/' + trajectory_file_dir
    traj_file_names = os.listdir(traj_files_full_path)
    process_list = []
    sim_number = 1
    for traj_file_name in traj_file_names:
        os.path.dirname(os.path.abspath(__file__))
        dir_and_traj_file_name = trajectory_file_dir + '/' + traj_file_name
        output_dir_and_file_name = output_dir_name + traj_file_name
        simulation_args = [config_file_name, dir_and_traj_file_name, output_dir_and_file_name, sim_number]
        process_list.append(Process(target=run_single_md_simulation, args=simulation_args))
        sim_number += 1

    for process in process_list:
        process.start()
    for process in process_list:
        process.join()

    print("All simulations finished!")


if __name__ == "__main__":
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [1, 1, 1]],
                              size=(2, 2, 3), symbol='Cu', pbc=(1, 1, 0))
    atoms.calc = EMT()
    #optimize_scaling(atoms, {'optimal_scaling': [], 'iterations_to_find_scaling': []})
    print(lattice_constant.optimize_scaling_using_simulation(atoms, {'potential': 'EMT', 'time_step': 5, 'temperature': 300, 'ensemble': 'NVE', 'step_number': 250}))

    #run_md_simulations("example_config.ini", 'Demo_multi_sim', 'Demo_multi_sim')
