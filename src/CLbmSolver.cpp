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

#include <limits>

#include "libmath/CMath.hpp"

template <class T>
CLbmSolver<T>::CLbmSolver(
        int id,
        CDomain<T> &domain,
        std::vector<Flag> boundaryConditions,
        T timestepSize,
        CVector<3, T> &gravitation,
        CVector<3, T> &drivenCavityVelocity,
        T viscosity,
        T maxGravitationDimLess,
        bool storeDensities,
        bool storeVelocities,
        bool doLogging) :
        id(id), domain(domain),
        boundaryConditions(boundaryConditions),
        timestepSize(timestepSize), gravitation(gravitation), drivenCavityVelocity(drivenCavityVelocity),
        viscosity(viscosity), maxGravitationDimLess(maxGravitationDimLess),
        storeDensities(storeDensities), storeVelocities(storeVelocities), doLogging(doLogging)
{
	bool viscosityGiven = this->viscosity > (T)0;
    T cellLength = domain.getLength()[0] / (T)domain.getSize()[0];

    while (true)
    {
    	// (4.11)
		gravitationDimLess = this->gravitation * ((this->timestepSize * this->timestepSize) / cellLength);
		drivenCavityVelocityDimLess = this->drivenCavityVelocity * (this->timestepSize / cellLength);

		/*
		 * If no viscosity was passed to this program, a default viscosity
		 * basing on a predefined reynolds number is set. This artificial
		 * viscosity can be further adapted if the timestep size has to adapted.
		 * The predefined reynolds number is always kept and stays constant.
		 */
		if (!viscosityGiven)
		{
			std::cerr << "No viscosity has been passed!" << std::endl;
			std::cerr << "Artificial viscosity is set!" << std::endl;

			this->viscosity = this->domain.getLength()[0] * this->drivenCavityVelocity[0] / REYNOLDS_DEFAULT;

			std::cerr << "viscosity:         " << this->viscosity << std::endl;
			std::cerr << "---------------------------------------" << std::endl;
		}

		// (4.9)
		viscosityDimLess = this->viscosity * (this->timestepSize / (cellLength * cellLength));

		/*
		 * If the dimension less gravity is larger than the specified maximum value,
		 * the simulation becomes unstable. In such a case, the timestep size is
		 * adapted accordingly and a valid dimension less gravity is set in the next iteration of this loop.
		 */
		if (gravitationDimLess.length() > this->maxGravitationDimLess + std::numeric_limits<T>::epsilon())
		{
			std::cerr << "Gravitation (dimension less) is too large so the simulation could get unstable!" << std::endl;
			std::cerr << "Timestep size is adapted!" << std::endl;
	        std::cerr << "old timestep size: " << this->timestepSize << std::endl;

	        // (4.12)
	        this->timestepSize = CMath<T>::sqrt((this->maxGravitationDimLess * cellLength) / this->gravitation.length());

	        std::cerr << "new timestep size: " << this->timestepSize << std::endl;
			std::cerr << "---------------------------------------" << std::endl;

			continue;
		}

		tau = (T)0.5 * ((T)6 * viscosityDimLess + (T)1);

		/*
		 * If tau is not within a certain range, the simulation becomes
		 * unstable. In such a case, the timestep size is adapted accordingly
		 * and a valid tau is set in the next iteration of this loop.
		 */
		if (tau - std::numeric_limits<T>::epsilon() < (T)TAU_LOWER_LIMIT || tau + std::numeric_limits<T>::epsilon() > (T)TAU_UPPER_LIMIT)
		{

			std::cerr << "Tau " << tau << " not within the range [" << TAU_LOWER_LIMIT <<"; " << TAU_UPPER_LIMIT << "]." << std::endl;
			std::cerr << "Timestep size is adapted!" << std::endl;
	        std::cerr << "old timestep size: " << this->timestepSize << std::endl;

			this->timestepSize = (cellLength * cellLength) * (((T)2 * TAU_DEFAULT - (T)1) / ((T)6 * this->viscosity));

	        std::cerr << "new timestep size: " << this->timestepSize << std::endl;
			std::cerr << "---------------------------------------" << std::endl;

			continue;
		}

		break;
    }

    tauInv = (T)1 / tau;
    tauInvTrt = (T)1 / ((T)0.5 + (T)3 / ((T)16 * tau - (T)8));

    int reynolds = this->domain.getLength()[0] * this->drivenCavityVelocity[0] / this->viscosity;

    if (this->doLogging)
    {
        std::cout << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
        std::cout << "id:                                      " << this->id << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "domain size:                             " << this->domain.getSize() << std::endl;
        std::cout << "domain length:                           " << this->domain.getLength() << std::endl;
        std::cout << "domain origin:                           " << this->domain.getOrigin() << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "timestep size:                           " << this->timestepSize << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "gravitation:                             " << this->gravitation << std::endl;
        std::cout << "gravitation (dimension less):            " << gravitationDimLess << std::endl;
        std::cout << "driven cavity velocity:                  " << this->drivenCavityVelocity << std::endl;
        std::cout << "driven cavity velocity (dimension less): " << drivenCavityVelocityDimLess << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "viscosity:                               " << this->viscosity << std::endl;
        std::cout << "viscosity (dimension less):              " << viscosityDimLess << std::endl;
        std::cout << "tau:                                     " << this->tau << std::endl;
        std::cout << "reynolds number (dimension less):        " << reynolds << std::endl;
        std::cout << "max gravitation length (dimension less): " << this->maxGravitationDimLess << std::endl;
        std::cout << "inv tau:                                 " << this->tauInv << std::endl;
        std::cout << "inv trt tau:                             " << this->tauInvTrt << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "store densities:                         " << this->storeDensities << std::endl;
        std::cout << "store velocities:                        " << this->storeVelocities << std::endl;
        std::cout << "---------------------------------------" << std::endl;
    }
}

template <class T>
CDomain<T>* CLbmSolver<T>::getDomain()
{
    return &domain;
}

template class CLbmSolver<float>;
template class CLbmSolver<double>;
