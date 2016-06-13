#!/usr/bin/python
import sys

code =  "# @ shell=/bin/bash\n";
code += "#\n";
code += "# @ error            = " + str(sys.argv[4]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".e.txt\n";
code += "# @ output           = " + str(sys.argv[4]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".o.txt\n";
code += "# @ job_type         = parallel\n";
code += "# @ requirements     = (Feature==\"gpu\")\n";
code += "# @ node_usage       = not_shared\n";
code += "#\n";
code += "# Number of nodes:\n";
code += "# @ node             = " + str(sys.argv[1]) + "\n";
code += "# Number of tasks/ranks/processes per node:\n";
code += "# @ tasks_per_node   = " + str(sys.argv[2]) + "\n";
code += "# Number of threads per task/rank/process:\n";
code += "# @ resources        = ConsumableCpus(" + str(sys.argv[3]) + ")\n";
code += "# @ wall_clock_limit = 00:05:00\n";
code += "#\n";
code += "# @ network.MPI      = sn_all,not_shared,us\n";
code += "# @ notification     = complete\n";
code += "# @ notify_user      = $(user)@rzg.mpg.de\n";
code += "# @ queue\n";
code += "\n";
code += "module unload intel/14.0\n";
code += "module load intel/16.0\n";
code += "module load cuda/7.5\n";
code += "module unload mpi.ibm/1.4.0\n";
code += "module load mpi.intel/5.1.3\n";
code += "module load netcdf-mpi/4.3.3.1\n";
code += "\n";
code += "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_HOME/lib/\n";
code += "export OMP_NUM_THREADS=" + str(sys.argv[3]) + "\n";
code += "\n";
code += "poe ${HOME}/workspace/lbm/bin/lbm ${HOME}/workspace/lbm/configurations/hydra_" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".xml\n";
	
jobscript = open("jobscript.sh", "w")
jobscript.write(code)
jobscript.close()

