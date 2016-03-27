#!/bin/bash

#SBATCH --job-name pyJ_mpiTest
#SBATCH -n 4 -N 4
#SBATCH --output mpi_testPython.out
#SBATCH --exclusive
#SBATCH -x node[1024,2001-2020]

## mpiexec python test.py
## mpirun -np 2 /home/yalcin/JDFTx_latest/build/jdftx -i test/li.in
mpirun -np 4 python ../test.py
## mpiexec -n 2 python test.py
