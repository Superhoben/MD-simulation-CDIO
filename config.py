from configparser import ConfigParser

def config_file():
  """ Create the configuration file

  Args:
      None

  Returns:
      None
  """
  config = ConfigParser()
  config['User inputs - Material ID'] = {'Structure':'Structure'}
  config['User inputs - Simulation parameter'] = {'Temperature (K)':'Temperature K','StepsNumber':'StepsNumber','Potential':'Potential'}

  # Write to our config file
  with open('config.ini','w') as config_file:
    config.write(config_file)

config_file()
