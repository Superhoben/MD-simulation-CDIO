"""This runs MD simulations."""
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
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


def save_configuration(atoms, output_file_name):
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_trajectory_files/'
    traj = Trajectory(path+output_file_name+'.traj', "w")
    traj.write(atoms)


def run_single_md_simulation(config_file: str, traj_file: str, output_name: str):
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
    config_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    config_data = ConfigParser()
    config_data.read(config_path+config_file)

    # Create atoms object for simulation
    traj_path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_trajectory_files/'
    traj = Trajectory(traj_path+traj_file)
    atoms = traj[0]

    # Set potential for simulation
    potential = config_data['SimulationSettings']['potential']
    if potential == "EMT":
        atoms.calc = EMT()
    elif potential == "LennardJones":
        atoms.calc = LennardJones()
    else:
        # TODO: implement running with other potentials, e.g.,:
        # atoms.calc = OtherPotential()
        raise Exception("Running calculations with potential '" + potential + "' is not implemented yet.")

    # Set dynamics module depending on simulation type
    ensemble = config_data['SimulationSettings']['ensemble']
    if ensemble == "NVE":
        MaxwellBoltzmannDistribution(atoms, temperature_K=int(config_data['SimulationSettings']['temperature']))
        dyn = VelocityVerlet(atoms, int(config_data['SimulationSettings']['time_step'])*units.fs)
    elif ensemble == "NVT":
        dyn = Langevin(atoms, timestep=config_data['SimulationSettings']['time_step']*units.fs,
                       temperature_K=config_data['SimulationSettings']['temperature'],
                       friction=(config_data['SimulationSettings']['friction'] or 0.005))
        MaxwellBoltzmannDistribution(atoms, temperature_K=config_data['SimulationSettings']['temperature'])
    else:
        # TODO: implement run_other_simulation, e.g.,:
        # atoms = run_other_simulation(atoms, config_data)
        raise Exception("Running calculations with ensemble '" + ensemble + "' is not implemented yet.")

    # This dict will contain output data of the simulation to be written into the output text file
    output_dict = {}
    # Attach recorders that calculate a certain property and store in the output_dict
    interval_to_record_temperature = int(config_data['RecordingIntervals']['record_temperature'])
    if interval_to_record_temperature:
        output_dict['temperature'] = []
        dyn.attach(calc_properties.calc_temp, interval_to_record_temperature, atoms, output_dict)

    interval_to_record_pressure = int(config_data['RecordingIntervals']['record_pressure'])
    if interval_to_record_pressure:
        output_dict['pressure'] = []
        dyn.attach(calc_properties.calc_pressure, interval_to_record_pressure, atoms, output_dict)

    interval_to_record_bulk_modulus = int(config_data['RecordingIntervals']['record_bulk_modulus'])
    if interval_to_record_bulk_modulus:
        output_dict['bulk_modulus'] = []
        dyn.attach(calc_bulk_properties.calc_bulk_modulus, interval_to_record_bulk_modulus, atoms, output_dict)

    interval_to_record_optimal_scaling = int(config_data['RecordingIntervals']['record_optimal_scaling'])
    if interval_to_record_optimal_scaling:
        output_dict['optimal_scaling'] = []
        output_dict['iterations_to_find_scaling'] = []
        dyn.attach(lattice_constant.optimize_scaling, interval_to_record_optimal_scaling, atoms, output_dict)

    interval_to_record_configuration = int(config_data['RecordingIntervals']['record_configuration'])
    if interval_to_record_configuration:
        dyn.attach(save_configuration, interval_to_record_configuration, atoms, output_name)

    # Run simulation with the attached recorders
    dyn.run(int(config_data['SimulationSettings']['step_number']))

    path = os.path.dirname(os.path.abspath(__file__)) + '/../Output_text_files/'
    with open(path + output_name + '.txt', 'w') as file:
        file.seek(0)
        json.dump(output_dict, file)

    return atoms


def run_md_simulation(config_file_list, trajectory_file_list):
    i=0
    for traj_file in trajectory_file_list:
        i=i+1
        run_single_md_simulation(config_file_list[0], traj_file, 'output'+str(i))


if __name__ == "__main__":
    run_single_md_simulation("config1.ini", 'mp-124.traj', 'out_put1')
