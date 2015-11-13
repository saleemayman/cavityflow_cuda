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

#ifndef CLBMSOLVERCPU_HPP
#define CLBMSOLVERCPU_HPP

#include "CLbmSolver.hpp"

// #include "gpukernels/copy_buffer_rect.cuh"
#include "CComm.hpp"
#include "CLbmSolverGPU.cuh"
#include "cpukernels/CLbmInitCPU.hpp"

template<typename T>
class CLbmSolverCPU : public CLbmSolver<T>
{
private:
    CLbmSolverGPU<T>* solverGPU;
//    CLbmInitCPU<T>* initKernelCPU;

    std::vector<T*> densityDistributions;
    std::vector<Flag*> flags;
    std::vector<T*> velocities;
    std::vector<T*> densities;
    
    /*
     * Maybe these members are not required if data copy operations during sync
     * operations and "inner cells operations"/"boundary cells operations" can
     * be achieved "in place", e.g. via memcpy().
     */
    /*
    T** getDensityDistributionsIntraHalo, setDensityDistributionsIntraHalo;
    T** getDensityDistributionsInterHalo, setDensityDistributionsInterHalo;
    */
    std::vector<CComm<T>*> commContainer;

public:
    CLbmSolverCPU();
    CLbmSolverCPU(
            int id,
            CVector<3, T> &globalLength,
            CDomain<T> &domain,
            std::vector<Flag> boundaryConditions,
            CLbmSolverGPU<T>* solverGPU,
            T timestepSize,
            CVector<3, T>& velocity,
            CVector<3, T>& acceleration,
            T viscocity,
            T maxVelocityDimLess,
            T maxAccelerationDimLess,
            bool storeDensities,
            bool storeVelocities,
            bool doLogging);
    ~CLbmSolverCPU();

    void simulationStepAlpha();
    void simulationStepAlphaRect(CVector<3, int> origin, CVector<3, int> size);
    void simulationStepBeta();
    void simulationStepBetaRect(CVector<3, int> origin, CVector<3, int> size);
    void getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* dst);
    void getDensityDistributions(T* dst);
    void setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* src);
    void setDensityDistributions(T* src);
    void getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* src);
    void getFlags(Flag* src);
    void setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* dst);
    void setFlags(Flag* dst);
    void getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* src);
    void getVelocities(T* src);
    void setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* dst);
    void setVelocities(T* dst);
    void getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* src);
    void getDensities(T* src);
    void setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* dst);
    void setDensities(T* dst);
};

#endif
