#!/bin/sh

# The MAC cluster uses SLURM as job scheduler: sbatch. As parameter it either
# only expects a job script (jobscript.sh) which contains all relevant
# information concerning the set up of the parallel environment and the
# execution of the parallel program or the relevant parameters can be passed by
# arguments.
# 
# We chose the first option. That's the reason why the XML configuration file
# and the job script are first generated by python before they are passed to
# sbatch.

let domainSizeX=652
let domainSizeY=652
let domainSizeZ=652
let domainLengthX=1
let domainLengthY=1
let domainLengthZ=1
let numOfSubdomainsX=1
let numOfSubdomainsY=4
let numOfSubdomainsZ=2
let loops=195600
let doBenchmark=1
let doLogging=0
let doValidation=0
let doVisualization=0
let visualizationRate=6520

let numOfGPUsPerNode=2
let numOfRanksPerNode=${numOfGPUsPerNode}
let numOfThreadsPerProcess=1
let numOfThreadsPerNode=numOfThreadsPerProcess*numOfRanksPerNode
let numOfNodes=(${numOfSubdomainsX}*${numOfSubdomainsY}*${numOfSubdomainsZ})/${numOfGPUsPerNode}
let numOfGPUs=${numOfGPUsPerNode}*${numOfNodes}
let numOfRanks=${numOfRanksPerNode}*${numOfNodes}
let numOfThreads=${numOfThreadsPerNode}*${numOfNodes}

# The place option forces distribution of chunks among nodes even insufficient
# data is offered 
# free:    As many chunks as possible are scheduled on one single node (depends
#          on job queue) before a second node is used
# scatter: Scatters chunks (numer after = sign) among different nodes, so
#          different chunks would be scheduled on the same node
# pack:    At maximum, one node is allocated
place="scatter"

echo "domain size:                   ["${domainSizeX}", "${domainSizeY}", "${domainSizeZ}"]"
echo "domain length:                 ["${domainLengthX}", "${domainLengthY}", "${domainLengthZ}"]"
echo "number of subdomains:          ["${numOfSubdomainsX}", "${numOfSubdomainsY}", "${numOfSubdomainsZ}"]"
echo "loops:                         "${loops}
echo "do benchmark:                  "${doBenchmark}
echo "do logging:                    "${doLogging}
echo "do validation:                 "${doValidation}
echo "do visualization:              "${doVisualization}
echo "visualization rate:            "${visualizationRate}
echo "GPUs per node:                 "${numOfGPUsPerNode}
echo "ranks per node:                "${numOfRanksPerNode}
echo "number of threads per process: "${numOfThreadsPerProcess}
echo "number of threads per node:    "${numOfThreadsPerNode}
echo "number of nodes:               "${numOfNodes}
echo "total number of GPUs:          "${numOfGPUs}
echo "total number of ranks:         "${numOfRanks}
echo "total number of threads:       "${numOfThreads}

python generate_configuration.py ${domainSizeX} ${domainSizeY} ${domainSizeZ} ${domainLengthX} ${domainLengthY} ${domainLengthZ} ${numOfSubdomainsX} ${numOfSubdomainsY} ${numOfSubdomainsZ} ${loops} ${doBenchmark} ${doLogging} ${doValidation} ${doVisualization} ${visualizationRate} ${numOfRanks} ${HOME} ${SCRATCH}
python generate_jobscript.py ${numOfNodes} ${numOfRanksPerNode} ${numOfThreadsPerProcess} ${HOME}
sbatch jobscript.sh

