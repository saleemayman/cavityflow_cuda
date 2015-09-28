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

#include "libmath/CVector.hpp"
#include "CDomain.hpp"

template<typename T>
class CLbmSolver
{
protected:
	int id;
	CDomain<T> domain;
	
	int boundaryConditions[3][2];
	T timestepSize;
	
	CVector<3,T> gravitation;
	CVector<4,T> drivenCavityVelocity;
	T viscocity;
	T massExchangeVactor;
	T renoldsNumber;
	T tau, tauInv, tauInvTrt;
	T maxSimGravitationLength;

public:
	CLbmSolver();
	virtual ~CLbmSolver();
	
	virtual void simulationStepAlpha();
	virtual void simulationStepAlphaRect();
	virtual void simulationStepBeta();
	virtual void simulationStepBetaRect();
	virtual void reset();
	virtual void reload();
	CDomain<T*> getDomain();
	virtual void getDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* dst);
	virtual void setDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* src);
	virtual void getFlags(CVector<3,int> &origin, CVector<3,int> &size, int* dst);
	virtual void setFlags(CVector<3,int> &origin, CVector<3,int> &size, int* src);
	virtual void getVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* dst);
	virtual void setVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* src);
	virtual void getDensities(CVector<3,int> &origin, CVector<3,int> &size, T* dst);
	virtual void setDensities(CVector<3,int> &origin, CVector<3,int> &size, T* src);
};

#endif
