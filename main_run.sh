#!/bin/bash

# Slurm sbatch options
#SBATCH -o main_run.log-%j
#SBATCH -c 40
#SBATCH -N 1

# Initialize Modules
source /etc/profile

#Load Julia module
export PATH=$HOME/julia-1.6.1/bin:$PATH

# Call your script as you would from the command line
julia --project=@. main.jl 
