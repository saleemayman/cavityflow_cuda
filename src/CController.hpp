/*
 * Copyright
 * 2010 Martin Schreiber
 * 2013 Arash Bakhtiari
 * 2016 Christoph Riesinger, Ayman Saleem
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CCONTROLLER_HPP
#define CCONTROLLER_HPP

#include <vector>

#include <cuda_runtime.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "libvis/CLbmVisualization.hpp"
#include "CConfiguration.hpp"
#include "CDomain.hpp"
#include "CLbmSolverCPU.hpp"
#include "CLbmSolverGPU.cuh"

/*
 * Class CConroller is responsible for controlling and managing of simulation and visualization
 * of a subdomain from the whole grid.
 */
template<typename T>
class CController
{
private:
    int id;
    CDomain<T> domain;
    std::vector<Flag> boundaryConditions;
    std::vector<CComm<T> > communication;
    CConfiguration<T>* configuration;
    CLbmSolverCPU<T>* solverCPU;
    CLbmSolverGPU<T>* solverGPU;
#ifdef USE_MPI
    std::vector<T*>* sendBuffers;
    std::vector<T*>* recvBuffers;
    MPI_Request* sendRequests;
    MPI_Request* recvRequests;
    std::vector<cudaStream_t>* streams;
#endif
    CLbmVisualization<T>* visualization;
    cudaStream_t defaultStream;
    int simulationStepCounter;

    CDomain<T> getDecomposedSubdomainCPU();
    CDomain<T> getDecomposedSubdomainGPU();
    void computeNextStep();
    void stepAlpha();
    void stepBeta();
    /*
#ifdef USE_MPI
    void syncAlpha();
    void syncBeta();
#endif
    */

public:
    CController(
            int id,
            CDomain<T> domain,
            std::vector<Flag> boundaryConditions,
            std::vector<CComm<T> > communication,
            CConfiguration<T>* configuration);
    ~CController();

    void setDrivenCavitySzenario();
    void run();

    int getId();
    CDomain<T>* getDomain();
    CLbmSolver<T>* getSolver();
};

#endif
