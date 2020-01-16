#!/bin/bash
#bash script for submitting a job to the sharcnet Graham queue

#SBATCH --ntasks=48     # number of MPI processes

#SBATCH --time=00-36:00:00  # time (DD-HH:MM:SS)

#SBATCH --mem-per-cpu=4096M  # memory per cpu in megabytes

#SBATCH --error=error-mod01_t00.25p.ca

#SBATCH --output=output-mod01_t00.25p.ca

#SBATCH --mail-type=END

#SBATCH --mail-user=smrgho@gmail.com

#SBATCH --account=def-carol      # UWO Physics & Astronomy

srun ./hdustparv2.10.bc input=Proj_velocity/mod01/mod01_t00.25p.inp
