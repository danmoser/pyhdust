#!/bin/bash
#SBATCH --job-name Proj
#SBATCH --output Proj-%j.out

# Walltime
#SBATCH --time=23:59:59

# Number of MPI tasks
#SBATCH -n 60

# Number of tasks per node
#SBATCH --ntasks-per-node 20  # sur fdr
#SBATCH -p fdr

echo "Running on: $SLURM_NODELIST"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "SLURM_NTASKS_PER_NODE=$SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "SLURM_NNODES=$SLURM_NNODES"
echo "SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE"

module purge
module load intel-gnu8-runtime/19.1.2.254
module load impi/2019.8.254
module list

I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so srun --mpi=pmi2 ./hdustparv2.10.bc input=mod01.inp 
