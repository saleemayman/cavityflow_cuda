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

#include <vector>

// #include "gpukernels/copy_buffer_rect.cuh"
#include "CComm.hpp"
#include "CLbmSolverGPU.cuh"

template<typename T>
class CLbmSolverCPU : public CLbmSolver<T>
{
private:
	CLbmSolverGPU<T> solverGPU;

	std::vector<T*> densityDistributions;
	std::vector<Flag*> flags;
	std::vector<T*> velocities;
	std::vector<T*> densities;
	
	T** getDensityDistributionsIntraHalo, setDensityDistributionsIntraHalo;
	T** getDensityDistributionsInterHalo, setDensityDistributionsInterHalo;
	std::vector<CComm<T>*> commContainer;

public:
	CLbmSolverCPU();
	CLbmSolverCPU(
			int id,
			CDomain<T> &domain,
			std::array<Flag,6> boundaryConditions,
			T timestepSize,
			CVector<3,T> &gravitation,
			CVector<4,T> &drivenCavityVelocity,
			T viscocity,
			T massExchangeFactor = MASS_EXCHANGE_FACTOR,
			T maxSimGravitationLength = MAX_SIM_GRAVITATION_LENGTH,
			T tau = TAU,
			bool storeDensities = false,
			bool storeVelocities = false);
	~CLbmSolverCPU();

	void simulationStepAlpha();
	void simulationStepAlphaRect(CVector<3,int> origin, CVector<3,int> size);
	void simulationStepBeta();
	void simulationStepBetaRect(CVector<3,int> origin, CVector<3,int> size);
	void reset();
	void reload();
	void getDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* hDensityDistributions);
	void setDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* hDensityDistributions);
	void getFlags(CVector<3,int> &origin, CVector<3,int> &size, int* hFlags);
	void setFlags(CVector<3,int> &origin, CVector<3,int> &size, int* hFlags);
	void getVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* hVelocities);
	void setVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* hVelocities);
	void getDensities(CVector<3,int> &origin, CVector<3,int> &size, T* hDensities);
	void setDensities(CVector<3,int> &origin, CVector<3,int> &size, T* hDensities);
};

#endif
