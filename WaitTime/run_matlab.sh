#!/bin/bash

# Slurm sbatch options
#SBATCH -o matlab-test.log-%j
#SBATCH -n 12
#SBATCH --constraint=xeon-g6

# Initialize Modules
source /etc/profile


# Run the script
matlab2020a -nodisplay -r "eval(pRUN('param_sweep', 12, 'grid')); exit"
