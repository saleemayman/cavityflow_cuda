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

#include "CLbmSolverCPU.hpp"

template <class T>
CLbmSolverCPU<T>::CLbmSolverCPU(
        int id,
        CVector<3, T> &globalLength,
        CDomain<T> &domain,
        std::vector<Flag> boundaryConditions,
        CLbmSolverGPU<T>* solverGPU,
        T timestepSize,
        CVector<3, T>& gravitation,
        CVector<3, T>& drivenCavityVelocity,
        T viscocity,
        T maxGravitationDimLess,
        bool storeDensities,
        bool storeVelocities,
        bool doLogging) :
        CLbmSolver<T>(id, globalLength,
        		domain, boundaryConditions,
                timestepSize, gravitation, drivenCavityVelocity,
                viscocity, maxGravitationDimLess,
                storeDensities, storeVelocities, doLogging),
        solverGPU(solverGPU)
{
    /*
     * Allocate memory for density distributions, density, flags and velocities.
     */


    /*
     * instantiate the cpu init-kernel class
     */

}

template <class T>
CLbmSolverCPU<T>::~CLbmSolverCPU()
{
//    delete initKernelCPU;
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlpha()
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlphaRect(CVector<3, int> origin, CVector<3, int> size)
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBeta()
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBetaRect(CVector<3, int> origin, CVector<3, int> size)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensityDistributions(T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensityDistributions(T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getFlags(Flag* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setFlags(Flag* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensities(T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensities(T* dst)
{
}

template class CLbmSolverCPU<float>;
template class CLbmSolverCPU<double>;
