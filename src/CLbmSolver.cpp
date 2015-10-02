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

#include "CLbmSolver.hpp"

template <class T>
CLbmSolver<T>::CLbmSolver(
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
		id(id), domain(domain),
		boundaryConditions(boundaryConditions), timestepSize(timestepSize),
		gravitation(gravitation), drivenCavityVelocity(drivenCavityVelocity), viscocity(viscocity),
		massExchangeFactor(massExchangeFactor), maxSimGravitationLength(maxSimGravitationLength), tau(tau),
		storeDensities(storeDensities), storeVelocities(storeVelocities)
{
	this->drivenCavityVelocity *= timestepSize;
	reynolds = domain.getLength()[0] * drivenCavityVelocity[0] / viscocity;
	tauInv = (T)1 / tau;
	tauInvTrt = (T)1 / ((T)0.5 + (T)3 / ((T)16 * tau - (T)8));

#if DEBUG
	std::cout << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
	std::cout << "id:                         " << this->id << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "domain size:                " << this->domain.getSize() << std::endl;
	std::cout << "domain length:              " << this->domain.getLength() << std::endl;
	std::cout << "domain origin:              " << this->domain.getOrigin() << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "timestep size:              " << this->timestepSize << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "gravitation:                " << this->gravitation << std::endl;
	std::cout << "driven cavity velocity:     " << this->drivenCavityVelocity << std::endl;
	std::cout << "viscocity:                  " << this->viscocity << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "mass exchange factor:       " << this->massExchangeFactor << std::endl;
	std::cout << "max sim gravitation length: " << this->maxSimGravitationLength << std::endl;
	std::cout << "tau:                        " << this->tau << std::endl;
	std::cout << "inv tau:                    " << this->tauInv << std::endl;
	std::cout << "inv trt tau:                " << this->tauInvTrt << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "store densities:            " << this->massExchangeFactor << std::endl;
	std::cout << "store velocities:           " << this->maxSimGravitationLength << std::endl;
	std::cout << "---------------------------------------" << std::endl;
#endif
}

template <class T>
CLbmSolver<T>::~CLbmSolver()
{
}

template <class T>
CDomain<T>* CLbmSolver<T>::getDomain()
{
	return &domain;
}
