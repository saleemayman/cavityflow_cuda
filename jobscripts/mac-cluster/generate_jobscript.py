#!/usr/bin/python
import sys

code =  "#!/bin/bash\n";
code += "#\n";
code += "#SBATCH -D " + str(sys.argv[4]) + "/workspace/lbm/\n";
code += "#SBATCH -o " + str(sys.argv[4]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".o.txt\n";
code += "#SBATCH -e " + str(sys.argv[4]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".e.txt\n";
code += "#SBATCH -J " + str(sys.argv[1]) + "\n";
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
code += "#SBATCH --mail-type=END\n";
code += "#SBATCH --mail-user=riesinge@in.tum.de\n";
code += "\n";
# code += "module load cuda/6.5\n";
# code += "module load mpi.ompi/1.6\n";
# code += "\n";
code += "mpirun -np " + str(int(sys.argv[1]) * int(sys.argv[2])) + " -ppn " + str(sys.argv[2]) + " ./bin/lbm configurations/mac-cluster.xml\n";
	
jobscript = open("jobscript.sh", "w")
jobscript.write(code)
jobscript.close()

