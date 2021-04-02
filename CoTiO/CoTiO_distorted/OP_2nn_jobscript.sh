#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --account=def-aparamek
#SBATCH --time=03:00:00
module load python/3.8.5
time python OP_2nn_main.py 1 4
time python OP_2nn_main.py 2 4
time python OP_2nn_main.py 3 4
time python OP_2nn_main.py 4 4
time python OP_2nn_main.py 5 4
time python OP_2nn_main.py 6 4
