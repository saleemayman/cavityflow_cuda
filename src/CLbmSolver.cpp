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
        CVector<3, T> &globalLength,
        CDomain<T> &domain,
        std::vector<Flag> boundaryConditions,
        T timestepSize,
        CVector<3, T> &velocity,
        CVector<3, T> &acceleration,
        T viscosity,
        T maxVelocityDimLess,
        T maxAccelerationDimLess,
        bool storeDensities,
        bool storeVelocities,
        bool doLogging) :
        id(id), globalLength(globalLength),
        domain(domain), boundaryConditions(boundaryConditions),
        timestepSize(timestepSize), velocity(velocity), acceleration(acceleration), viscosity(viscosity),
        maxVelocityDimLess(maxVelocityDimLess), maxAccelerationDimLess(maxAccelerationDimLess),
        storeDensities(storeDensities), storeVelocities(storeVelocities), doLogging(doLogging)
{
	bool limitedByVelocity = false;
	bool limitedByAcceleration = false;
	T cellLength = domain.getLength()[0] / (T)domain.getSize()[0];
	T oldTimestepSize;

    if (this->doLogging)
    {
        std::cout << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
        std::cout << "id:                                      " << this->id << std::endl;
        std::cout << "---------------------------------------" << std::endl;
    }

	/*
	 * If no viscosity was passed to this program, a default viscosity
	 * basing on a predefined reynolds number is set. This artificial
	 * viscosity can be further adapted if the timestep size has to adapted.
	 * The predefined reynolds number is always kept and stays constant.
	 */
	if (this->viscosity <= (T)0)
	{
		if (this->doLogging)
		{
			std::cout << "No viscosity has been passed!" << std::endl;
			std::cout << "Artificial viscosity is set!" << std::endl;
		}

		// definition reynolds number
		this->viscosity = this->globalLength.max() * this->velocity[0] / REYNOLDS_DEFAULT;

		if (this->doLogging)
		{
			std::cout << "viscosity:         " << this->viscosity << std::endl;
			std::cout << "---------------------------------------" << std::endl;
		}
	}

    while (true)
    {
		if (this->doLogging)
		{
			std::cout << "New iteration of finding timestep size started." << std::endl;
			std::cout << "---------------------------------------" << std::endl;
        }

    	// (4.11)
		velocityDimLess = this->velocity * (this->timestepSize / cellLength);
		accelerationDimLess = this->acceleration * ((this->timestepSize * this->timestepSize) / cellLength);
		// (4.9)
		viscosityDimLess = this->viscosity * (this->timestepSize / (cellLength * cellLength));

		/*
		 * If the dimension less velocity is larger than the specified maximum
		 * value, the simulation becomes unstable. In such a case, the timestep
		 * size is adapted accordingly and a valid dimension less velocity is
		 * set in the next iteration of this loop.
		 */
		if (velocityDimLess.length() > this->maxVelocityDimLess + std::numeric_limits<T>::epsilon())
		{
			limitedByVelocity= true;
			oldTimestepSize = this->timestepSize;
	        this->timestepSize = (this->maxVelocityDimLess * cellLength) / this->velocity.length();

			if (this->doLogging)
			{
				std::cout << "Velocity (dimension less) is too large so the simulation could get unstable!" << std::endl;
				std::cout << "Timestep size is adapted!" << std::endl;
				std::cout << "old timestep size: " << oldTimestepSize << std::endl;
		        std::cout << "new timestep size: " << this->timestepSize << std::endl;
				std::cout << "---------------------------------------" << std::endl;
	        }

			continue;
		}

		/*
		 * If the dimension less acceleration is larger than the specified
		 * maximum value, the simulation becomes unstable. In such a case, the
		 * timestep size is adapted accordingly and a valid dimension less
		 * acceleration is set in the next iteration of this loop.
		 */
		if (accelerationDimLess.length() > this->maxAccelerationDimLess + std::numeric_limits<T>::epsilon())
		{
			limitedByAcceleration = true;
			oldTimestepSize = this->timestepSize;
	        // (4.12)
	        this->timestepSize = CMath<T>::sqrt((this->maxAccelerationDimLess * cellLength) / this->acceleration.length());

			if (this->doLogging)
			{
				std::cout << "Acceleration (dimension less) is too large so the simulation could get unstable!" << std::endl;
				std::cout << "Timestep size is adapted!" << std::endl;
				std::cout << "old timestep size: " << oldTimestepSize << std::endl;
		        std::cout << "new timestep size: " << this->timestepSize << std::endl;
				std::cout << "---------------------------------------" << std::endl;
	        }

			continue;
		}

		// (4.7)
		tau = (T)0.5 * ((T)6 * viscosityDimLess + (T)1);

		/*
		 * If tau is not within a certain range, the simulation becomes
		 * unstable. In such a case, the timestep size is adapted accordingly
		 * and a valid tau is set in the next iteration of this loop.
		 */
		if (tau - std::numeric_limits<T>::epsilon() < (T)TAU_LOWER_LIMIT || tau + std::numeric_limits<T>::epsilon() > (T)TAU_UPPER_LIMIT)
		{
			oldTimestepSize = this->timestepSize;
			// 4.9 & 4.10
	        this->timestepSize = (cellLength * cellLength) * (((T)2 * TAU_DEFAULT - (T)1) / ((T)6 * this->viscosity));

			if (this->doLogging)
			{
				std::cout << "Tau " << tau << " not within the range [" << TAU_LOWER_LIMIT <<"; " << TAU_UPPER_LIMIT << "]." << std::endl;
				std::cout << "Timestep size is adapted!" << std::endl;
				std::cout << "old timestep size: " << oldTimestepSize << std::endl;
				std::cout << "new timestep size: " << this->timestepSize << std::endl;
				std::cout << "---------------------------------------" << std::endl;
	        }

			if ((limitedByVelocity || limitedByAcceleration) && this->timestepSize > oldTimestepSize) {
		        std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
		        std::cerr << "No valid timestep size could be determined which satisfies" << std::endl;
		        std::cerr << "- viscosity:                         " << this->viscosity << std::endl;
		        std::cerr << "- max velocity (dimension less):     " << this->maxVelocityDimLess << std::endl;
		        std::cerr << "- max acceleration (dimension less): " << this->maxAccelerationDimLess << std::endl;
		        std::cerr << "- tau:                               " << tau << std::endl;
		        std::cerr << "so the simulation stays stable!" << std::endl;

		        exit (EXIT_FAILURE);
			} else {
				continue;
			}
		}

		break;
    }

    tauInv = (T)1 / tau;
    // tauInv = (T)1 / ((T)0.5 + (T)3 / ((T)16 * tau - (T)8));

    int reynolds = this->domain.getLength()[0] * this->velocity[0] / this->viscosity;

    if (this->doLogging)
    {
        std::cout << "global length (without halo):      " << this->globalLength << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "domain size (without halo):        " << this->domain.getSize() << std::endl;
        std::cout << "domain length (without halo):      " << this->domain.getLength() << std::endl;
        std::cout << "domain origin (without halo):      " << this->domain.getOrigin() << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "timestep size:                     " << this->timestepSize << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "velocity:                          " << this->velocity << std::endl;
        std::cout << "velocity (dimension less):         " << velocityDimLess << std::endl;
        std::cout << "acceleration:                      " << this->acceleration << std::endl;
        std::cout << "acceleration (dimension less):     " << accelerationDimLess << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "viscosity:                         " << this->viscosity << std::endl;
        std::cout << "viscosity (dimension less):        " << viscosityDimLess << std::endl;
        std::cout << "tau:                               " << this->tau << std::endl;
        std::cout << "reynolds number (dimension less):  " << reynolds << std::endl;
        std::cout << "max velocity (dimension less):     " << this->maxVelocityDimLess << std::endl;
        std::cout << "max acceleration (dimension less): " << this->maxAccelerationDimLess << std::endl;
        std::cout << "inv tau:                           " << this->tauInv << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "store densities:                   " << this->storeDensities << std::endl;
        std::cout << "store velocities:                  " << this->storeVelocities << std::endl;
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
