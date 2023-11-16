"""This handel the user's API key."""
import os
from configparser import ConfigParser


# Configuration settings
config_file = "User_API_key.ini"
section = "User"

# Create and load the configuration file
config = ConfigParser()


# Function to prompt the user for the API key
def prompt_for_api_key():
    while True:
        api_key = input("Enter your API key (32 Character): ").strip()
        if len(api_key) == 32:
            # Save the API key in the configuration file for future use
            config[section] = {"api_key": api_key}
            path = os.path.dirname(os.path.abspath(__file__)) + '/../User_API_key/'
            with open(path+config_file, "w") as configfile:
                config.write(configfile)
            return api_key
        else:
            print("Vaild key is required.")


# Function to get user's API key
def get_api_key():
    # Check if the API key is already saved in the configuration file
    path = os.path.dirname(os.path.abspath(__file__)) + '/../User_API_key/'
    full_path = path+config_file

    # Config.read(config_file_path+config_file)
    if os.path.isfile(full_path):
        config.read(full_path)
        if section in config and "api_key" in config[section]:
            saved_api_key = config[section]["api_key"]
            confirmation = input(f"Is this your API key? (yes/no): {saved_api_key} ").strip().lower()
            if confirmation == "yes":
                print(f"Using saved API key: {saved_api_key}")
                return saved_api_key
            elif confirmation == "no":
                print("API key not confirmed.")
            else:
                print("Invalid option, try again")
    return prompt_for_api_key()


if __name__ == "__main__":
    api_key = get_api_key()
    print(f"Using API key: {api_key}")
