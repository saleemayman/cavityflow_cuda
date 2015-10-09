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

#ifndef CLBMSOLVER_HPP
#define CLBMSOLVER_HPP

#include <vector>

#include "libmath/CVector.hpp"
#include "CDomain.hpp"
#include "common.h"

#define MASS_EXCHANGE_FACTOR       (1.0)
#define MAX_SIM_GRAVITATION_LENGTH (0.0001)
#define NUM_LATTICE_VECTORS        19
#define TAU                        (0.953575)

#define STORE_VELOCITIES           (true)
#define STORE_DENSITIES            (true)

template<typename T>
class CLbmSolver
{
protected:
	int id;
	CDomain<T> domain;
	
	/*
	 * Vector for the boundary conditions of the six faces of a cuboid.
	 * boundaryConditions[0]: boundary condition for right face (smallest x-coordinate)
	 * boundaryConditions[1]: boundary condition for left face (largest x-coordinate)
	 * boundaryConditions[2]: boundary condition for top face (smallest y-coordinate)
	 * boundaryConditions[3]: boundary condition for bottom face (largest y-coordinate)
	 * boundaryConditions[4]: boundary condition for front face (smallest z-coordinate)
	 * boundaryConditions[5]: boundary condition for back face (largest z-coordinate)
	 */
	std::vector<Flag> boundaryConditions;
	T timestepSize;
	
	CVector<3, T> gravitation;
	CVector<4, T> drivenCavityVelocity;
	T viscocity;
	T massExchangeFactor;
	T maxSimGravitationLength;
	T tau;

	bool storeDensities;
	bool storeVelocities;

	T reynolds, tauInv, tauInvTrt;

public:
	CLbmSolver();
	CLbmSolver(
			int id,
			CDomain<T> &domain,
			std::vector<Flag> boundaryConditions,
			T timestepSize,
			CVector<3, T> &gravitation,
			CVector<4, T> &drivenCavityVelocity,
			T viscocity,
			T massExchangeFactor = MASS_EXCHANGE_FACTOR,
			T maxSimGravitationLength = MAX_SIM_GRAVITATION_LENGTH,
			T tau = TAU,
			bool storeDensities = false,
			bool storeVelocities = false);
	virtual ~CLbmSolver() {}
	
	virtual void simulationStepAlpha() {}
	virtual void simulationStepAlphaRect(CVector<3, int> origin, CVector<3, int> size) {}
	virtual void simulationStepBeta() {}
	virtual void simulationStepBetaRect(CVector<3, int> origin, CVector<3, int> size) {}
	virtual void reset() {}
	CDomain<T>* getDomain();
	virtual void getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* dst) {}
	virtual void getDensityDistributions(T* dst) {}
	virtual void setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* src) {}
	virtual void setDensityDistributions(T* src) {}
	virtual void getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* dst) {}
	virtual void getFlags(Flag* dst) {}
	virtual void setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* src) {}
	virtual void setFlags(Flag* src) {}
	virtual void getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* dst) {}
	virtual void getVelocities(T* dst) {}
	virtual void setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* src) {}
	virtual void setVelocities(T* src) {}
	virtual void getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* dst) {}
	virtual void getDensities(T* dst) {}
	virtual void setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* src) {}
	virtual void setDensities(T* src) {}
};

#endif
