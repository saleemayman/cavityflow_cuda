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

#include <vector>

#include <vector_types.h>

template<typename T>
class CLbmSolverGPU : public CLbmSolver<T>
{
private:
	using CLbmSolver<T>::domain;
	using CLbmSolver<T>::storeDensities;
	using CLbmSolver<T>::storeVelocities;

	T* densityDistributions;
	Flag* flags;
	T* velocities;
	T* densities;
	
	/*
	 * Six slots for the halo layers of the six faces of a cuboid.
	 * [g/s]etDensityDistributionsHalo[0]: halo layer for right face
	 * [g/s]etDensityDistributionsHalo[1]: halo layer for left face
	 * [g/s]etDensityDistributionsHalo[2]: halo layer for top face
	 * [g/s]etDensityDistributionsHalo[3]: halo layer for bottom face
	 * [g/s]etDensityDistributionsHalo[4]: halo layer for front face
	 * [g/s]etDensityDistributionsHalo[5]: halo layer for back face
	 */
	std::array<T*, 6> getDensityDistributionsHalo;
	std::array<T*, 6> setDensityDistributionsHalo;

	/*
	 * Four slots for the parallel setup (threads per block) for the four different GPU kernels.
	 * threadsPerBlock[0]: number of threads per block for kernel lbm_init()
	 * threadsPerBlock[1]: number of threads per block for kernel lbm_alpha()
	 * threadsPerBlock[2]: number of threads per block for kernel lbm_beta()
	 */
	std::array<dim3,3> threadsPerBlock;

public:
	CLbmSolverGPU();
	CLbmSolverGPU(
			int id,
			CDomain<T> &domain,
			std::array<Flag,6> boundaryConditions,
			T timestepSize,
			CVector<3, T> &gravitation,
			CVector<4, T> &drivenCavityVelocity,
			T viscocity,
			T massExchangeFactor = MASS_EXCHANGE_FACTOR,
			T maxSimGravitationLength = MAX_SIM_GRAVITATION_LENGTH,
			T tau = TAU,
			bool storeDensities = false,
			bool storeVelocities = false);
	~CLbmSolverGPU();

	void simulationStepAlpha();
	void simulationStepAlphaRect(CVector<3, int> origin, CVector<3, int> size);
	void simulationStepBeta();
	void simulationStepBetaRect(CVector<3, int> origin, CVector<3, int> size);
	void reset();
	void getDesityDistributions(int i, T* hDensityDistributions);
	void getDesityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions);
	void setDesityDistributions(int i, T* hDensityDistributions);
	void setDesityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions);
	void getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* hFlags);
	void setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* hFlags);
	void getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* hVelocities);
	void setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* hVelocities);
	void getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* hDensities);
	void setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* hDensities);
};

#endif
