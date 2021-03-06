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

#ifndef CLBMSOLVERGPU_CUH
#define CLBMSOLVERGPU_CUH

#include "CLbmSolver.hpp"

template<typename T>
class CLbmSolverGPU : public CLbmSolver<T>
{
private:
    using CLbmSolver<T>::id;
    using CLbmSolver<T>::domain;
    using CLbmSolver<T>::configuration;
    using CLbmSolver<T>::velocityDimLess;
    using CLbmSolver<T>::accelerationDimLess;
    using CLbmSolver<T>::storeDensities;
    using CLbmSolver<T>::storeVelocities;
    using CLbmSolver<T>::tauInv;

    T* densityDistributions;
    Flag* flags;
    T* velocities;
    T* densities;

    /*
     * Four slots for the parallel setup (threads per block) for the four different GPU kernels.
     * threadsPerBlock[0]: number of threads per block for kernel lbm_init()
     * threadsPerBlock[1]: number of threads per block for kernel lbm_alpha()
     * threadsPerBlock[2]: number of threads per block for kernel lbm_beta()
     */
    std::vector<dim3> threadsPerBlock;

public:
    CLbmSolverGPU();
    CLbmSolverGPU(
            int id,
            CDomain<T> &domain,
            std::vector<Flag> boundaryConditions,
            CConfiguration<T>* configuration);
    ~CLbmSolverGPU();

    using CLbmSolver<T>::simulationStepAlpha;
    using CLbmSolver<T>::simulationStepBeta;
    using CLbmSolver<T>::getDensityDistributions;
    using CLbmSolver<T>::setDensityDistributions;

    void simulationStepAlpha(cudaStream_t* stream);
    void simulationStepAlpha();
    void simulationStepAlpha(CVector<3, int> origin, CVector<3, int> size, cudaStream_t* stream);
    void simulationStepAlpha(CVector<3, int> origin, CVector<3, int> size);
    void simulationStepBeta(cudaStream_t* stream);
    void simulationStepBeta();
    void simulationStepBeta(CVector<3, int> origin, CVector<3, int> size, cudaStream_t* stream);
    void simulationStepBeta(CVector<3, int> origin, CVector<3, int> size);
    void getDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, T* hDensityDistributions, cudaStream_t* stream);
    void getDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, T* hDensityDistributions);
    void getDensityDistributions(T* hDensityDistributions, cudaStream_t* stream);
    void getDensityDistributions(T* hDensityDistributions);
    void setDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, Direction direction, T* hDensityDistributions, cudaStream_t* stream);
    void setDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, Direction direction, T* hDensityDistributions);
    void setDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, T* hDensityDistributions, cudaStream_t* stream);
    void setDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, T* hDensityDistributions);
    void setDensityDistributions(T* hDensityDistributions, cudaStream_t* stream);
    void setDensityDistributions(T* hDensityDistributions);
    void getFlags(CVector<3, int>& origin, CVector<3, int>& size, Flag* hFlags);
    void getFlags(Flag* hFlags);
    void setFlags(CVector<3, int>& origin, CVector<3, int>& size, Flag* hFlags);
    void setFlags(Flag* hFlags);
    void getVelocities(CVector<3, int>& origin, CVector<3, int>& size, T* hVelocities);
    void getVelocities(T* hVelocities);
    void setVelocities(CVector<3, int>& origin, CVector<3, int>& size, T* hVelocities);
    void setVelocities(T* hVelocities);
    void getDensities(CVector<3, int>& origin, CVector<3, int>& size, T* hDensities);
    void getDensities(T* hDensities);
    void setDensities(CVector<3, int>& origin, CVector<3, int>& size, T* hDensities);
    void setDensities(T* hDensities);
};

#endif
