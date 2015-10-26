#!/bin/bash
#
#SBATCH -D /home/hpc/pr63so/lu32dec/workspace/lbm/
#SBATCH -o /home/hpc/pr63so/lu32dec/workspace/lbm/results/2.o.txt
#SBATCH -e /home/hpc/pr63so/lu32dec/workspace/lbm/results/2.e.txt
#SBATCH -J 1
#SBATCH --get-user-env
#
#SBATCH --partition=nvd
#Total number of tasks/ranks/processes:
#SBATCH --ntasks=2
#Number of tasks/ranks/processes per node:
#SBATCH --ntasks-per-node=2
#Number of threads per task/rank/process:
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#
#SBATCH --mail-type=END
#SBATCH --mail-user=riesinge@in.tum.de

mpirun -np 2 -ppn 2 ./bin/lbm configurations/mac-cluster.xml