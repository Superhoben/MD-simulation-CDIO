#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A LiU-compute-2023-32
#SBATCH --reservation devel
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --exclusive
#
export NSC_MODULE_SILENT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
source ~/MD-simulation-CDIO/activate_conda.sh
conda activate new_MD

time python3 ~/MD-simulation-CDIO/Simulation/run_md_simulation.py

echo "job completed"
