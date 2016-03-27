#!/bin/bash

#SBATCH --job-name pyJ_mpiTest
#SBATCH -n 4 -N 4
#SBATCH --output mpi_test.out
#SBATCH --exclusive
#SBATCH -x node[1024,2001-2020]

## mpiexec python test.py
mpirun -np 4 /home/yalcin/JDFTx_latest/build/jdftx -i li.in
