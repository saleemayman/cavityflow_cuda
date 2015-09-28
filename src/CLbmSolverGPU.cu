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

#include "CLbmSolverGPU.cuh"

template <class T>
CLbmSolverGPU<T>::CLbmSolverGPU()
{
}

template <class T>
CLbmSolverGPU<T>::~CLbmSolverGPU()
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlpha()
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlphaRect()
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBeta()
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBetaRect()
{
}

template <class T>
void CLbmSolverGPU<T>::reset()
{
}

template <class T>
void CLbmSolverGPU<T>::reload()
{
}

template <class T>
void CLbmSolverGPU<T>::getDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* src)
{
}

template <class T>
void CLbmSolverGPU<T>::getFlags(CVector<3,int> &origin, CVector<3,int> &size, int* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setFlags(CVector<3,int> &origin, CVector<3,int> &size, int* src)
{
}

template <class T>
void CLbmSolverGPU<T>::getVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* src)
{
}

template <class T>
void CLbmSolverGPU<T>::getDensities(CVector<3,int> &origin, CVector<3,int> &size, T* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setDensities(CVector<3,int> &origin, CVector<3,int> &size, T* src)
{
}

template class CLbmSolverGPU<float>;
template class CLbmSolverGPU<double>;
