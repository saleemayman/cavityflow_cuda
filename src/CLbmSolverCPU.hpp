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

#include "CComm.hpp"
#include "CLbmSolverGPU.cuh"
#include "cpukernels/CLbmInitCPU.hpp"
#include "cpukernels/CLbmAlphaCPU.hpp"
#include "cpukernels/CLbmBetaCPU.hpp"

template<typename T>
class CLbmSolverCPU : public CLbmSolver<T>
{
private:
    using CLbmSolver<T>::id;
    using CLbmSolver<T>::globalLength;
    using CLbmSolver<T>::domain;
    using CLbmSolver<T>::velocity;
    using CLbmSolver<T>::velocityDimLess;
    using CLbmSolver<T>::acceleration;
    using CLbmSolver<T>::accelerationDimLess;
    using CLbmSolver<T>::storeDensities;
    using CLbmSolver<T>::storeVelocities;
    using CLbmSolver<T>::doLogging;
    using CLbmSolver<T>::tauInv;

    CLbmSolverGPU<T>* solverGPU;
    CLbmInitCPU<T>* initLbmCPU;
    CLbmAlphaCPU<T>* alphaLbmCPU;
    CLbmBetaCPU<T>* betaLbmCPU;

    CVector<3, int> hollowCPULeftLimit;
    CVector<3, int> hollowCPURightLimit;

    std::vector<T> densityDistributions;
    std::vector<Flag> flags;
    std::vector<T> velocities;
    std::vector<T> densities;

    int domainCellsCPUWithHalo;    

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

    using CLbmSolver<T>::simulationStepAlpha;
    using CLbmSolver<T>::simulationStepBeta;
    using CLbmSolver<T>::getDensityDistributions;
    using CLbmSolver<T>::setDensityDistributions;

    void simulationStepAlpha();
    void simulationStepAlpha(CVector<3, int> origin, CVector<3, int> size);
    void simulationStepBeta();
    void simulationStepBeta(CVector<3, int> origin, CVector<3, int> size);
    void getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* dst);
    void getDensityDistributions(T* dst);
    void setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* src);
    void setDensityDistributions(T* src);
    void getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* dst);
    void getFlags(Flag* dst);
    void setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* src);
    void setFlags(Flag* src);
    void getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* dst);
    void getVelocities(T* dst);
    void setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* src);
    void setVelocities(T* src);
    void getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* dst);
    void getDensities(T* dst);
    void setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* src);
    void setDensities(T* srsrcc);

    CVector<3, int> getHollowCPULeftLimits();
    CVector<3, int> getHollowCPURightLimits();
};

#endif
