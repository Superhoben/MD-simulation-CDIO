"""Import the configparser library which creates the config file."""
import os
from configparser import ConfigParser


def config_file(file_name='default_config', ensemble='NVE', temperature=500, potential='EMT',
                step_number=5000, time_step=5, friction=0.005, record_energy = 0,
                record_cohesive_energy = 0, record_temperature = 0, record_pressure = 0, 
                record_configuration = 0, record_bulk_modulus = 0, record_optimal_scaling = 0,
                record_elastic = 0):
    """Create the configuration file

    Args:
    	file_name(string): Name of the file to create
        ensemble(string): Ensemble to use in simulation
        temperature(int): Initial temperature in simulation
        potential(string): Potential to use in simulation
        step_number(int): Number of steps to use in simulation
        time_step(int): Time step in fs to use in simulation
        friction(float): Friction for NVT simulation
        record_energy(int): Interval to record energy in simulation
        record_cohesive_energy(int): Interval to record cohesive energy in simulation
        record_temperature(int): Interval to record temperature in simulation
        record_pressure(int): Interval to record pressure in simulation
        record_configuration(int): Interval to record configuration in simulation
        record_bulk_modulus(int): Interval to record bulk modulus in simulation
        record_optimal_scaling(int): Interval to record optimal scaling in simulation
        record_elastic(int): Interval to record elastic properties in simulation

    Returns:
        None
    """
    if ensemble == "Select an Option":
        ensemble = "NVE"
    if potential == "Select an Option":
        potential = "EMT"
    
    config = ConfigParser()
    config['SimulationSettings'] = {'ensemble': ensemble, 'temperature': temperature,
                                    'potential': potential, 'step_number': step_number,
                                    'time_step': time_step, 'friction': friction}

    config['RecordingIntervals'] = {'record_energy': record_energy,
    				    'record_cohesive_energy': record_cohesive_energy,
                                    'record_temperature': record_temperature,
                                    'record_pressure': record_pressure,
                                    'record_configuration': record_configuration,
                                    'record_bulk_modulus': record_bulk_modulus,
                                    'record_optimal_scaling': record_optimal_scaling,
                                    'record_elastic': record_elastic}

    # Write to config file
    path = os.path.dirname(os.path.abspath(__file__)) + '/../Input_config_files/'
    with open(path+file_name+'.ini', 'w') as config_file:
        config.write(config_file)


def temperature_parameter():
    """Fetch the user's temperature parameter.

    Args:
        None

    Returns:
        (float): the input temperature
    """
    while True:
        temperature = input("Choose the Simulation temperature in Kelvin?  (hint: interval between [0k - 10000k]) \n")
        try:
            temperature = float(temperature)
            if 0 <= temperature <= 10000:
                return temperature
            else:
                print("Invalid temperature, try again")
        except ValueError:
            print("Invalid input, please enter a numeric value for temperature.")


def Steps_number():
    """Fetch the user's steps number parameter.

    Args:
        None

    Returns:
        (float): the input steps number
    """
    while True:
        Steps_number = input("Choose the Simulation steps number?  (hint: pick a step around ≈1 fsec, interval=[1-10000] (ex: 5 femtosecond)) \n")
        try:
            Steps_number = float(Steps_number)
            if 1 <= Steps_number <= 10000:
                return Steps_number
            else:
                print("Invalid steps number, try again")
        except ValueError:
            print("Invalid input, please enter a numeric value for Steps_number.")


def potential():
    """Fetch the user's potential.

    Args:
        None

    Returns:
        (string): the desired potential
    """
    available_potential_list = ["Lennard Jones", "EMT", "bla bla"]
    while True:
        print("Choose from the exsisting potentials: " + ", ".join(available_potential_list))
        potential = input("Enter the potential name: ").strip()
        if potential in available_potential_list:
            return potential
        else:
            print("Unavailable potential, try again")


def choose_ensemble():
    """Fetch the user's ensemble.

    Args:
    None

    Returns:
    (string): the desired ensemble
    """
    available_ensemble_list = ["NVE", "NVT", "NPT"]
    while True:
        print("Choose from the exsisting ensembles: " + ", ".join(available_ensemble_list))
        ensemble = input("Enter the ensemble name: ").strip()
        if ensemble in available_ensemble_list:
            return ensemble
        else:
            print("Unavailable ensemble, try again")


def material_id():
    """Fetch the user's material_id parameter.

    Args:
        None

    Returns:
        (string): the material id
    """
    material_id = input("What is the Material ID?  (hint: mp-123456 or cif file)   \n")
    return material_id


if __name__ == "__main__":
    config_file('config_parser_test')
