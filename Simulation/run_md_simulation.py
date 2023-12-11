"""This runs MD simulations."""
import os
import sys
import shutil
from math import ceil
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.visualize import view
from ase.md.verlet import VelocityVerlet
from ase.md.npt import NPT
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
from Gather_data.hypothetical_materials import mix_materials

def print_and_increase_progress(progress, sim_number):
    """Prints to the terminal how far a simulations has come"""
    if sim_number:
        print("Simulation ", sim_number, ": Progress "+str(progress[0])+"%")
    else:
        print("Progress "+str(progress[0])+"%")
    progress[0] += 10


def run_single_md_simulation(config_file: str, traj_file: str, output_name: str, sim_number=0, enable_prints=True):
    """Run md simulation for a single trajectory file, with parameters specified in config.

    Snapshoots of the atom state during the simulation will be saved in Output_trajectory_files/output_name.traj
    and calculated properties in Output_text_files/output_name.txt

    Args:
        config_file(str): Name of file with parameters to use in simulation.
        traj_file(str): Name of trajectory file with the atoms object to use in simulation
        output_name(str): Name of file to write results to
        sim_number(int): Used by the function print_and_increase_progress specify the simulation

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
#    elif ensemble == "NPT":
#        dyn = NPT(atoms, timestep=int(simulation_settings['time_step'])*units.fs,
#                  temperature_K=int(simulation_settings['temperature']),
#                  ttime=float(simulation_settings['ttime']), pfactor=float(simulation_settings['pfactor']))
#        MaxwellBoltzmannDistribution(atoms, temperature_K=2*int(simulation_settings['temperature']))
#        Stationary(atoms)
    else:
        # TODO: implement other ensembles
        raise Exception("Running calculations with ensemble '" + ensemble + "' is not implemented yet.")

    # This dict will contain output data of the simulation to be written into the output text file

    # Approximate lattice constant
    d = calc_properties.approx_lattice_constant(atoms)

    output_dict = {}
    output_dict["config_file"] = [config_file]

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
        output_dict["debye_temperature"] = []
        dyn.attach(calc_bulk_properties.calc_bulk_modulus, interval_to_record_bulk_modulus, atoms, output_dict)
        dyn.attach(calc_properties.time_average_of_debye_temperature, interval_to_record_bulk_modulus, atoms, output_dict)

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
        dyn.attach(calc_properties.lindemann_criterion, interval_to_record_lindemann_criterion, output_dict, d)

    interval_to_record_self_diffusion_coefficient = int(recording_intervals['record_self_diffusion_coefficient'])
    if interval_to_record_lindemann_criterion:
        output_dict['self_diffusion_coefficient'] = []
        dyn.attach(calc_properties.self_diffusion_coefficent, interval_to_record_self_diffusion_coefficient, output_dict, interval_to_record_self_diffusion_coefficient*int(simulation_settings['time_step']))

    if enable_prints:
        progress = [0]
        ten_percent_interval = int(0.1 * float(simulation_settings['step_number']))
        dyn.attach(print_and_increase_progress, ten_percent_interval, progress, sim_number)

    # Run simulation with the attached recorders
    if sim_number:
        print("Simulation ", sim_number, " started")
    else:
        print("Simulation started")
    dyn.run(int(simulation_settings['step_number']))

    # The simulation is finished, calculate certain properties based on the collected data from the simulation
    if interval_to_record_energy:
        output_dict['specific_heat_capacity'] = [calc_properties.calculate_specific_heat(atoms, config_data, output_dict)]

    if int(recording_intervals['record_mean_square_displacement']):
        output_dict['mean_square_displacement'][0] = 0

    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/'
    with open(path + output_name + '.txt', 'w') as file:
        file.seek(0)
        json.dump(output_dict, file)

    if sim_number:
        print("Simulation ", sim_number, " finished!")
    else:
        print("Simulation finished!")
    return atoms


def queue_simulations(config_file_name, trajectory_file_dir, output_dir_name, 
                      process_list=[], execute=True, sim_number=[1], enable_prints=True):
    """For every material in trajectory_file_dir, append a simulation function object to process list

    Args:
        config_file_name (str): The name of the file with the simulations setting, it will be used for all simulations
        trajectory_file_dir (str): The directory containing the trajectory files which give the initial atoms objects
        output_dir_name (str): The name of the directories where the results will be saved, the trajectory files will 
        be saved in Output_trajectory_files/output_dir_name and the text files in Output_text_fils/output_dir_name
        process_list (list, optional): The list in which simulation function objects will be appended
        execute (bool, optional): Runs all simulation function objects from the process list in parallel
        sim_number (list, optional): The number of the simulation, is counted up for every new simulation function object
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

    for traj_file_name in traj_file_names:
        os.path.dirname(os.path.abspath(__file__))
        dir_and_traj_file_name = trajectory_file_dir + '/' + traj_file_name
        output_dir_and_name = output_dir_name + '/' + traj_file_name[:-5]
        sim_args = [config_file_name, dir_and_traj_file_name, output_dir_and_name, sim_number[0], enable_prints]
        process_list.append(Process(target=run_single_md_simulation, args=sim_args))
        sim_number[0] += 1

    if execute:
        for process in process_list:
            process.start()
        for process in process_list:
            process.join()
        if enable_prints:
            print("All simulations finished!")


def is_float(element: any) -> bool:
    if element is None: 
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False


def get_starting_num_string(file_name):
    """Get the number at the start of a string, intended to be used with certain file names

    Args:
        file_name(str): Using high_throughput_mix_and_simulate the start of a files name 
            will be the concentration of the mixed in material.

    Return:
        concentration(int): Given in percentage of the mixed in material
    """
    i = 0
    while is_float(file_name[:i+1]):
        i += 1
    return file_name[:i]


def summerize_text_files(dir_name):
    """Summarizes the main results of each file in dir_name and save it in the parent directory

    Args:
        dir_name(str): Should be given relative to the Output_text_files directory
    """
    dir_path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/' + dir_name
    dir_contents = os.listdir(dir_path)
    text_files = []
    for content in dir_contents:
        if os.path.isfile(dir_path+'/'+content) and content.endswith('.txt'):
            text_files.append(content)

    summary_dict = {'origin_files': text_files}
    for text_file in text_files:
        opened_file = open(dir_path+'/'+text_file, 'r')
        data_dict = json.load(opened_file)
        opened_file.close()
        for key in data_dict:
            if not (key in summary_dict):
                summary_dict[key] = []
            data = data_dict[key]
            number_of_measurements_to_include = ceil(0.2*len(data_dict[key]))
            if number_of_measurements_to_include == 1:
                summary_dict[key].append(data[0])
            else:
                summary_dict[key].append(np.mean(data[-number_of_measurements_to_include:]))

    with open(dir_path+'.txt', 'w') as file:
        file.seek(0)
        json.dump(summary_dict, file)


def summerize_traj_files(dir_name):
    """Put the last atom state of each trajectory file in the directory Output_trajectory_files/dir_name in one traj files

    Args:
        dir_name (str): The name of the directory with the trajectory files which should be given relative to
            Output_trajectory_files. The created trajectory file can be found as Output_trajectory_files/dir_name.traj
    """
    dir_path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files/' + dir_name
    dir_contents = os.listdir(dir_path)
    dir_contents.sort(key=get_starting_num_string)
#    parent_of_dir = '/'.join(dir_path.split('/')[:-1])
    traj_writer = Trajectory(dir_path+'.traj', 'w')
    for content in dir_contents:
        if os.path.isfile(dir_path+'/'+content) and content.endswith('.traj'):
            atoms = Trajectory(dir_path+'/'+content)[-1]
            traj_writer.write(atoms)


def high_throughput_mix_and_simulate(config_file, input_traj_dir, element_to_mix_in, mixing_concentrations,
                                     output_dir_name, enable_prints=True):
    """Mixes element_to_mix_in into every material in input_traj_dir then simulates all the mixes in parallell

    Args:
        config_file (str): Name of the config file which will be used to get the simulation settings
        input_traj_dir (str): The path to the directory with input trajectory files
        element_to_mix_in (): The chemical symbol of the elements the we want to mix in such as 'Cu', 'Ni', etc.
        mixing_concentrations (list(float)): List of the concentration with which we will mix in the element
        output_dir_name (str): Name for the main output directory that will be created in Output_trajectory_files
            and Output_text_files
    """
    base_path = os.path.dirname(os.path.abspath(__file__)) + '/../'
    path_new_output_txt_dir = base_path + 'Output_text_files/' + output_dir_name
    path_new_output_traj_dir = base_path + 'Output_trajectory_files/' + output_dir_name
    if os.path.exists(path_new_output_txt_dir):
        shutil.rmtree(path_new_output_txt_dir)
    if os.path.exists(path_new_output_traj_dir):
        shutil.rmtree(path_new_output_traj_dir)
    os.mkdir(path_new_output_txt_dir)
    os.mkdir(path_new_output_traj_dir)
    input_traj_dir_path = base_path + 'Input_trajectory_files/' + input_traj_dir
    input_traj_dir_contents = os.listdir(input_traj_dir_path)
    traj_file_names = []
    for content in input_traj_dir_contents:
        if os.path.isfile(input_traj_dir_path+'/'+content) and content.endswith('.traj'):
            traj_file_names.append(content)
    if len(traj_file_names) == 0:
        print("No traj files in the input traj dir, no alloys created, no simulations run")
        return

    config_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    config_data = ConfigParser()
    config_data.read(config_path+config_file)
    config_data['SimulationSettings']["record_configuration"] = '0'

    mix_dirs = []
    process_list = []
    sim_number = [1] # Simulation number is encapsulated in a list in order to become mutable
    for traj_file_name in traj_file_names:
        mix_dir_name =element_to_mix_in + '_mixed_into_' + traj_file_name[:-5]
        mix_materials(input_traj_dir+'/'+traj_file_name, element_to_mix_in, mixing_concentrations,  
                      input_traj_dir + '/' + mix_dir_name)
        mix_dirs.append(output_dir_name+'/'+mix_dir_name)
        queue_simulations(config_file, input_traj_dir+'/'+mix_dir_name, output_dir_name+'/'+mix_dir_name,
                          process_list, False, sim_number, enable_prints)

    for process in process_list:
        process.start()
    for process in process_list:
        process.join()

    for mix_dir in mix_dirs:
        summerize_traj_files(mix_dir)
        summerize_text_files(mix_dir)

    if enable_prints:
        print("All simulations finished!")


if __name__ == "__main__":
    mixing_concentrations = np.arange(0, 1, 0.01)
    high_throughput_mix_and_simulate("supercomputer_config.ini", 'Supercomputer_demo', 'Ag', mixing_concentrations,
                                     'Demo_multi_sim', False)
#    traj = Trajectory('/Users/gustavwassback/Documents/CDIO/MD-simulation-CDIO/Gather_data/../' +
#                      'Output_trajectory_files/Demo_multi_sim/Ni_mixed_into_1728_atoms_of_mp-30.traj', 'r')
#    view(traj)
