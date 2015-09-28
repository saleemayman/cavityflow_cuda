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

#include <list>
#include <vector>

#include "CComm.hpp"
#include "CLbmSolverGPU.cuh"

template<typename T>
class CLbmSolverCPU : public CLbmSolver<T>
{
private:
	CLbmSolverGPU<T> solverGPU;

	std::vector<T*> densityDistributions;
	std::vector<int*> flags;
	std::vector<T*> velocities;
	std::vector<T*> densities;
	
	std::vector<CComm<T>*> commContainer;

public:
	CLbmSolverCPU();
	~CLbmSolverCPU();

	void simulationStepAlpha();
	void simulationStepAlphaRect();
	void simulationStepBeta();
	void simulationStepBetaRect();
	void reset();
	void reload();
	void getDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* dst);
	void setDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* src);
	void getFlags(CVector<3,int> &origin, CVector<3,int> &size, int* dst);
	void setFlags(CVector<3,int> &origin, CVector<3,int> &size, int* src);
	void getVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* dst);
	void setVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* src);
	void getDensities(CVector<3,int> &origin, CVector<3,int> &size, T* dst);
	void setDensities(CVector<3,int> &origin, CVector<3,int> &size, T* src);
};

#endif
