#!/bin/bash -l

# Example jobscript to run an R MPI parallel job

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=12:0:0

# Request 1 gigabyte of RAM per process.
#$ -l mem=2G

# Request 15 gigabytes of TMPDIR space per node (default is 10 GB)
#$ -l tmpfs=10G

# Set the name of the job.
#$ -N DR_sim

# Select the MPI parallel environment with 32 processes
#$ -pe smp 32

# Set the working directory to somewhere in your scratch space.  This is
# necessary because the compute nodes cannot write to your $HOME
# NOTE: this directory must exist.
# Replace "<your_UCL_id>" with your UCL user ID
#$ -wd /home/zczlcga/Scratch/doubly_robust_PAIC/DR_sim


# Load the R module
module unload -f compilers mpi gcc-libs
module load r/new


# Run our MPI job. GERun is our wrapper for mpirun, which launches MPI jobs  
R --no-save < /home/zczlcga/Scratch/doubly_robust_PAIC/DR_sim/main.R > DR_sim_out
