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
		CDomain<T> &domain,
		std::array<Flag,6> boundaryConditions,
		T timestepSize,
		CVector<3,T> &gravitation,
		CVector<4,T> &drivenCavityVelocity,
		T viscocity,
		T massExchangeFactor,
		T maxSimGravitationLength,
		T tau,
		bool storeDensities,
		bool storeVelocities) :
		CLbmSolver<T>(id, domain,
				boundaryConditions, timestepSize,
				gravitation, drivenCavityVelocity, viscocity,
				massExchangeFactor, maxSimGravitationLength, tau,
				storeDensities, storeVelocities)
{
}

template <class T>
CLbmSolverCPU<T>::~CLbmSolverCPU()
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlpha()
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlphaRect(CVector<3,int> origin, CVector<3,int> size)
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBeta()
{
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBetaRect(CVector<3,int> origin, CVector<3,int> size)
{
}

template <class T>
void CLbmSolverCPU<T>::reset()
{
}

template <class T>
void CLbmSolverCPU<T>::getDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getFlags(CVector<3,int> &origin, CVector<3,int> &size, Flag* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setFlags(CVector<3,int> &origin, CVector<3,int> &size, Flag* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensities(CVector<3,int> &origin, CVector<3,int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensities(CVector<3,int> &origin, CVector<3,int> &size, T* dst)
{
}

template class CLbmSolverCPU<float>;
template class CLbmSolverCPU<double>;
