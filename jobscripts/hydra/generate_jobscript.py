#!/usr/bin/python
import sys

code =  "# @ shell=/bin/bash\n";
code += "#\n";
code += "# @ error            = " + str(sys.argv[5]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".e.txt\n";
code += "# @ output           = " + str(sys.argv[5]) + "/workspace/lbm/results/" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".o.txt\n";
code += "# @ job_type         = parallel\n";
code += "# @ requirements     = (Feature==\"gpu\")\n";
code += "# @ node_usage       = not_shared\n";
code += "#\n";
code += "# Number of nodes:\n";
code += "# @ node             = " + str(sys.argv[1]) + "\n";
code += "# Number of tasks/ranks/processes per node:\n";
code += "# @ tasks_per_node   = " + str(sys.argv[2]) + "\n";
code += "# Number of threads per task/rank/process:\n";
code += "# @ resources        = ConsumableCpus(" + str(sys.argv[3]) + ") ConsumableMemory(" + str(sys.argv[4]) + "gb)\n";
code += "# @ wall_clock_limit = 24:00:00\n";
code += "#\n";
code += "# @ network.MPI      = sn_all,not_shared,us\n";
code += "# @ notification     = complete\n";
code += "# @ notify_user      = $(user)@rzg.mpg.de\n";
code += "# @ queue\n";
code += "\n";
code += "module load cuda/6.5\n";
code += "module load mpi.ibm/1.4.0\n";
code += "module load netcdf-mpi/4.3.3.1\n";
code += "\n";
code += "poe ${HOME}/workspace/lbm/bin/lbm ${HOME}/workspace/lbm/configurations/hydra_" + str(int(sys.argv[1]) * int(sys.argv[2])) + ".xml\n";
	
jobscript = open("jobscript.sh", "w")
jobscript.write(code)
jobscript.close()

