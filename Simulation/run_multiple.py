"""This runs multiple simulations."""
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from Simulation import run_md_simulation

if len(sys.argv) <= 3:
    raise Exception("Not enough arguments. Config file name, trajectory file directory and output directory name needs to be specified.")    

queue_simulations(sys.argv[1], sys.argv[2], sys.argv[3])
