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

#include "libvis/ILbmVisualization.hpp"
#include "CComm.hpp"
#include "CConfiguration.hpp"
#include "CDomain.hpp"
#include "CLbmSolverCPU.hpp"
#include "CLbmSolverGPU.cuh"

#define MPI_TAG_ALPHA_SYNC  0
#define MPI_TAG_BETA_SYNC   1

/*
 * Class CConroller is responsible for controlling and managing of simulation and visualization
 * of a subdomain from the whole grid.
 */
template<typename T>
class CController
{
private:
    int id;
    CConfiguration<T>* configuration;
    CDomain<T> domain;
    // CLbmSolverCPU<T> *solverCPU;
    CLbmSolverGPU<T> *solverGPU;
    ILbmVisualization<T>* cLbmVisualization;
    std::vector<Flag> boundaryConditions;
    std::vector<CComm<T> > communication;
    int simulationStepCounter;

    T vector_checksum;

    /*
    CCL::CPlatforms* cPlatforms;
    CCL::CPlatform* cPlatform;
    CCL::CContext* cContext;
    CCL::CDevices* cDevices;
    CCL::CDevice* cDevice;
    CCL::CCommandQueue* cCommandQueue;
    */

    void computeNextStep();
    void syncAlpha();
    void syncBeta();

public:
    CController(
    		int id,
    		CDomain<T> domain,
    		std::vector<Flag> boundaryConditions,
    		std::vector<CComm<T> > communication,
    		CConfiguration<T>* configuration);
    ~CController();

    void setDrivenCavitySzenario();
    int run();

    int getId();
    CDomain<T>* getDomain();
    CLbmSolver<T>* getSolver();
};

#endif
