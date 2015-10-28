#!/bin/bash
#
#SBATCH -D /home/hpc/pr63so/lu32dec/workspace/lbm/
#SBATCH -o /home/hpc/pr63so/lu32dec/workspace/lbm/results/4.o.txt
#SBATCH -e /home/hpc/pr63so/lu32dec/workspace/lbm/results/4.e.txt
#SBATCH -J 2
#SBATCH --get-user-env
#
#SBATCH --partition=nvd
#Total number of tasks/ranks/processes:
#SBATCH --ntasks=4
#Number of tasks/ranks/processes per node:
#SBATCH --ntasks-per-node=2
#Number of threads per task/rank/process:
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#
#SBATCH --mail-type=END
#SBATCH --mail-user=riesinge@in.tum.de

mpirun -np 4 -ppn 2 ./bin/lbm configurations/mac-cluster.xml
