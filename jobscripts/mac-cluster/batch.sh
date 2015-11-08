#!/bin/sh

# The MAC cluster uses SLURM as job scheduler: sbatch. As parameter it either
# only expects a job script (jobscript.sh) which contains all relevant
# information concerning the set up of the parallel environment and the
# execution of the parallel program or the relevant parameters can be passed by
# arguments.
# 
# We chose the first option. That's the reason why the job scripts are first
# generated by python before they are passed to sbatch.

let maxNumOfNodes=4
let numOfGPUsPerNode=2
let numOfRanksPerNode=${numOfGPUsPerNode}
let numOfThreadsPerProcess=1
let numOfThreadsPerNode=numOfThreadsPerProcess*numOfRanksPerNode

# The place option forces distribution of chunks among nodes even insufficient
# data is offered 
# free:    As many chunks as possible are scheduled on one single node (depends
#          on job queue (e.g. 24 for S queue)) before a second node is used
# scatter: Scatters chunks (numer after = sign) among different nodes, so
#          different chunks would be scheduled on the same node
# pack:    At maximum, one node is allocated
place="scatter"

echo "maxNumOfNodes:          "${maxNumOfNodes}
echo "numOfGPUsPerNode:       "${numOfGPUsPerNode}
echo "numOfRanksPerNode:      "${numOfRanksPerNode}
echo "numOfThreadsPerProcess: "${numOfThreadsPerProcess}
echo "numOfThreadsPerNode:    "${numOfThreadsPerNode}

for((numOfNodes=4; numOfNodes<=maxNumOfNodes; numOfNodes++))
do
	let numOfGPUs=numOfGPUsPerNode*numOfNodes
	let numOfRanks=numOfRanksPerNode*numOfNodes
	
	echo "Scheduling job on "${numOfNodes}" nodes, spawning "${numOfRanks}" MPI ranks to utilize "${numOfGPUs}" GPUs"

	python generate_configuration.py ${numOfNodes} ${numOfRanksPerNode} ${numOfThreadsPerProcess} ${HOME} ${SCRATCH}
	python generate_jobscript.py ${numOfNodes} ${numOfRanksPerNode} ${numOfThreadsPerProcess} ${HOME}
	sbatch jobscript.sh
done

