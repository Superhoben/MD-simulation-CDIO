"""This creates alloys and runs simulations on them."""
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
from Simulation import run_md_simulation

if len(sys.argv) <= 5:
    raise Exception("Not enough arguments. Config file name, trajectory file directory, element to mix in, mixing concentrations and output directory name needs to be specified.")    

high_throughput_mix_and_simulate(sys.argv[1], sys.argv[2], sys.argv[3], [float(i) for i in sys.argv[4:-1]], sys.argv[-1]):
