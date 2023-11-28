"""This handel the user's API key."""
import os
from configparser import ConfigParser
from tkinter import simpledialog
from tkinter import messagebox

# Configuration settings
config_file = "API_key.ini"
section = "User"

# Create and load the configuration file
config = ConfigParser()


# Function to prompt the user for the API key
def prompt_for_api_key():
    
    x = True
    while x:
        
        api_key = simpledialog.askstring(title="Input API key for materials project", 
                                    prompt="API key:") 
        if api_key is None:
            return
        
        if len(api_key) != 32:
            # askyesno returns true if yes and false if no
            if messagebox.askyesno("Warning", "The API-key is not 32 tokens long. Are you sure it is correct?"):
                break
        else:
            x = False
        
        
    config[section] = {"api_key": api_key}
    path = os.path.dirname(os.path.abspath(__file__)) + '/../API_key/'
    with open(path+config_file, "w") as configfile:
        config.write(configfile)
    return api_key
           


# Function to get user's API key
def get_api_key():
    # Check if the API key is already saved in the configuration file
    path = os.path.dirname(os.path.abspath(__file__)) + '/../API_key/'
    full_path = path+config_file
    # Config.read(config_file_path+config_file)
    if os.path.isfile(full_path):
        config.read(full_path)
        if section in config and "api_key" in config[section]:
            saved_api_key = config[section]["api_key"]
            return saved_api_key 
    return "API key doesnt exist"


if __name__ == "__main__":
    api_key = get_api_key()
    print(f"Using API key: {api_key}")
