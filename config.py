from configparser import ConfigParser

config = ConfigParser()

def temperature_parameter():
    """ Fetch the user's temperature parameter 

    Args:
        None

    Returns:
        (int): the input temperature
    """
    while True:
        temperature = input("Choose the Simulation temperature?  (hint: interval between [0k - 10000k]) \n")
        try:
            temperature = int(temperature)
            # if temperature >= 0 and temperature <= 10000:
            if 0 <= temperature <= 10000:
                return  temperature
            else:
                print("Invalid temperature, try again")
        except ValueError:
            print("Invalid input, please enter a numeric value for temperature.")
                
        
    


def config_file():
    
  """ Creat the configuration file

  Args:
      None

  Returns:
      config file with "ini" format
  """
  print("Welcome to MD simulation Software \nPlease provide the program with the parameters and conditions needed \n")
  print("What is the Material ID?  (hint: mp-123456)  ")
  config['User inputs - Material ID'] = {'Structure':input() }

  config['User inputs - Simulation parameter'] = {'Temperature (K)': temperature_parameter(),
                      'StepsNumber': 'StepsNumber',
                      'Potential': 'Potential'}


  # Write to our config file
  with open('config.ini','w') as config_file:
    config.write(config_file)
    






config_file()






