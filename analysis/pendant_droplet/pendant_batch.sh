#!/bin/bash
#SBATCH -J pendant
#SBATCH -o pendant.o%j
#SBATCH -e pendant.e%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=rjnoh@utexas.edu

module list
pwd
date

julia -t 32 -- analysis_new/pendant_droplet/Main.jl

date
