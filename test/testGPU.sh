#!/bin/bash

#SBATCH --job-name pyJ_mpiTest
#SBATCH --output gpu_test.out
#SBATCH --exclusive
#SBATCH --gres=gpu
## mpiexec python test.py
jdftx/build/jdftx_gpu -i test/li.in
