#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --account=def-aparamek
#SBATCH --time=00:30:00
module load python/3.8.5
time python OP_1nn_main.py 4
