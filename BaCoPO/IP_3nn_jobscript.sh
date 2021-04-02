#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --account=def-aparamek
#SBATCH --time=01:30:00
module load python/3.8.5
time python IP_3nn_main.py 1 4
time python IP_3nn_main.py 2 4
time python IP_3nn_main.py 3 4
