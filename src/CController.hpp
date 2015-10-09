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

#include "libcuda/CCL.hpp"
#include "libvis/ILbmVisualization.hpp"
#include "CComm.hpp"
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
    int _UID; ///< Unique ID of each controller
    CDomain<T> _domain; ///< Domain data
    ILbmVisualization<T>* cLbmVisualization; ///< Visualization class
    CLbmSolverGPU<T> *cLbmPtr;
    std::vector<Flag> boundaryConditions; ///< Boundary conditions. First index specifies the dimension and second the upper or the lower boundary.
    std::vector<CComm<T>*> _comm_container; ///< A std::Vector containing all the communcation objects for the subdomain
    T vector_checksum;
    CCL::CPlatforms* cPlatforms;
    CCL::CPlatform* cPlatform;
    CCL::CContext* cContext;
    CCL::CDevices* cDevices;
    CCL::CDevice* cDevice;
    CCL::CCommandQueue* cCommandQueue;

    int simulationStepCounter;

    void outputDD(int dd_i);
    int initLBMSolver();

public:
    CController(int UID, CDomain<T> domain, std::vector<Flag> boundaryConditions);
    ~CController();

    void syncAlpha();
    void syncBeta();
    void computeNextStep();
    int run();
    void addCommunication(CComm<T>* comm);
    /*
     * TODO
     * This piece of code should be obsolete since the ghost layer sizes are now
     * set implicitly by CLbmSolver depending on the domain size.
     */
    // void addCommToSolver();
    void setGeometry();
    CLbmSolver<T>* getSolver() const;
    void setSolver(CLbmSolverGPU<T>* lbmPtr);
    CDomain<T> getDomain() const;
    int getUid() const;
};

#endif
