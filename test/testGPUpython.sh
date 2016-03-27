#!/bin/bash

#SBATCH --job-name pyJ_GPU
#SBATCH -n 1
#SBATCH --output testGPUpython.out
#SBATCH --exclusive
#SBATCH --gres=gpu

## mpiexec python test.py
## mpirun -np 2 /home/yalcin/JDFTx_latest/build/jdftx -i test/li.in
python ../test.py --gpu
## mpiexec -n 2 python test.py
