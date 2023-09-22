from configparser import ConfigParser

# get all functions from the module ConfigParser
config = ConfigParser()

# Read our config file
config.read("config.ini")
#user = input("Please provide a valid Material ID: ")
print(config.sections())
print(config['User inputs - Material ID']["structure"])