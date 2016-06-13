#!/usr/bin/python
import sys

code =  "#!/bin/bash\n";
code += "#\n";
code += "#SBATCH -D " + str(sys.argv[4]) + "/workspace/lbm/\n";
code += "#SBATCH -o " + str(sys.argv[4]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".o.txt\n";
code += "#SBATCH -e " + str(sys.argv[4]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".e.txt\n";
code += "#SBATCH -J lbm\n";
code += "#SBATCH --get-user-env\n";
code += "#\n";
code += "#SBATCH --partition=nvd\n";
code += "#Total number of tasks/ranks/processes:\n";
code += "#SBATCH --ntasks=" + str(int(sys.argv[1]) * int(sys.argv[2])) + "\n";
code += "#Number of tasks/ranks/processes per node:\n";
code += "#SBATCH --ntasks-per-node=" + str(sys.argv[2]) + "\n";
code += "#Number of threads per task/rank/process:\n";
code += "#SBATCH --cpus-per-task=" + str(sys.argv[3]) + "\n";
code += "#SBATCH --time=24:00:00\n";
code += "#\n";
code += "#SBATCH --mail-type=end\n";
code += "#SBATCH --mail-user=riesinge@in.tum.de\n";
code += "\n";
# code += "mpiexec.hydra -genv OMP_NUM_THREADS " + str(sys.argv[3]) + " -np " + str(int(sys.argv[1]) * int(sys.argv[2])) + " -ppn " + str(sys.argv[2]) + " ./bin/lbm configurations/mac-cluster_" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".xml\n";
code += "mpiexec.hydra -genv OMP_NUM_THREADS " + str(sys.argv[3]) + " -np " + str(int(sys.argv[1]) * int(sys.argv[2])) + " -ppn " + str(sys.argv[2]) + " nvprof -o lbm_%p.nvprof ./bin/lbm configurations/mac-cluster_" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".xml\n";
	
jobscript = open("jobscript.sh", "w")
jobscript.write(code)
jobscript.close()

