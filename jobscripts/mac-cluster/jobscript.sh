#!/bin/bash
#
#SBATCH -D /home/hpc/pr63so/lu32dec/workspace/lbm/
#SBATCH -o /home/hpc/pr63so/lu32dec/workspace/lbm/results/8.o.txt
#SBATCH -e /home/hpc/pr63so/lu32dec/workspace/lbm/results/8.e.txt
#SBATCH -J lbm
#SBATCH --get-user-env
#
#SBATCH --partition=nvd
#Total number of tasks/ranks/processes:
#SBATCH --ntasks=8
#Number of tasks/ranks/processes per node:
#SBATCH --ntasks-per-node=2
#Number of threads per task/rank/process:
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#
#SBATCH --mail-type=end
#SBATCH --mail-user=riesinge@in.tum.de

mpiexec.hydra -genv OMP_NUM_THREADS 8 -np 8 -ppn 2 nvprof -o lbm_%p.nvprof ./bin/lbm configurations/mac-cluster_8.xml
