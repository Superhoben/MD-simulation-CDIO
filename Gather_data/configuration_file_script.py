"""Import the configparser library which creates the config file."""
from configparser import ConfigParser

config = ConfigParser()


def config_file(ensemble, materialID, temperature, stepsnumber, potential):
    """Create the configuration file.

    Args:
        ensemble(string): the ensemble used in the simulation
        materialID(string): specifies which structure is used
        temperature(string): the temperature of the simulation
        stepsnumber(string): the number of steps used in the simulation
        potential(string): the type of potential that is used

    Returns:
        None
    """

    config['User inputs - Ensembles'] = {'Ensemble': ensemble}
    config['User inputs - Material ID'] = {'Structure': materialID}
    config['User inputs - Simulation parameter'] = {'Temperature': temperature,
                                                    'StepsNumber': stepsnumber, 'Potential': potential}
    # Write to our config file
    with open('config1.ini', 'w') as config_file:
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
    # here it is better to call a file (cif file for example) which will be created from Gustav function some how
    # which then fetch the data and store it into a seperate cif file, it was Rickard advice :/
    material_id = input("What is the Material ID?  (hint: mp-123456 or cif file)   \n")
    return material_id


if __name__ == "__main__":
    config_file()
