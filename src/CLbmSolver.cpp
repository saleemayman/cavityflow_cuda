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

#include "libmath/CMath.hpp"

template <class T>
CLbmSolver<T>::CLbmSolver(
        int id,
        CDomain<T> &domain,
        std::vector<Flag> boundaryConditions,
        T timestepSize,
        CVector<3, T> &gravitation,
        CVector<4, T> &drivenCavityVelocity,
        T viscocity,
        T massExchangeFactor,
        T maxGravitationDimLess,
        bool storeDensities,
        bool storeVelocities,
        bool doLogging) :
        id(id), domain(domain),
        boundaryConditions(boundaryConditions),
        timestepSize(timestepSize), gravitation(gravitation), drivenCavityVelocity(drivenCavityVelocity),
        viscocity(viscocity), massExchangeFactor(massExchangeFactor), maxGravitationDimLess(maxGravitationDimLess),
        storeDensities(storeDensities), storeVelocities(storeVelocities), doLogging(doLogging)
{
    if (this->doLogging)
    {
        std::cout << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
        std::cout << "id:                                      " << this->id << std::endl;
        std::cout << "---------------------------------------" << std::endl;
    }

    T cellLength = domain.getLength()[0] / (T)domain.getSize()[0];

    /*
     * If an invalid timestep size was passed to this program, a default
     * timestep size is set.
     */
    if (this->timestepSize < (T)0)
    {
        this->timestepSize = (cellLength * cellLength) * ((T)2 * this->tau - (T)1) / ((T)6 * this->viscocity * CMath<T>::sqrt(this->massExchangeFactor));

        if (this->doLogging)
        {
            std::cout << "No valid timestep size has been passed. A default timestep size set!" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }
    }

    gravitationDimLess = this->gravitation * ((this->timestepSize * this->timestepSize) / cellLength);

    /*
     * Limit the gravitation parameter for the simulation to avoid large
     * velocities and thus an unstable simulation.
     */
    if (gravitationDimLess.length() >= this->maxGravitationDimLess)
    {
        this->timestepSize = CMath<T>::sqrt((this->maxGravitationDimLess * cellLength) / this->gravitation.length());
        gravitationDimLess = this->gravitation * ((this->timestepSize * this->timestepSize) / cellLength);

        if (this->doLogging)
        {
            std::cout << "Gravitation has been limited to avoid large velocities and thus instability!" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }
    }

    tau = (T)0.5 * (this->timestepSize * this->viscocity * CMath<T>::sqrt(this->massExchangeFactor) * (T)6) / (cellLength * cellLength) + (T)0.5;

    if (tau < (T)TAU_LOWER_LIMIT || tau > (T)TAU_UPPER_LIMIT)
    {
        std::cerr <<    "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
        std::cerr << "Tau " << tau << " not within the boundary [" << TAU_LOWER_LIMIT <<"; " << TAU_UPPER_LIMIT <<"]." << std::endl;
        std::cerr << "Simulation becomes unstable!" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "---------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
    }

    drivenCavityVelocityDimLess = this->drivenCavityVelocity * this->timestepSize;

    tauInv = (T)1 / tau;
    tauInvTrt = (T)1 / ((T)0.5 + (T)3 / ((T)16 * tau - (T)8));

    if (this->doLogging)
    {
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
        std::cout << "viscocity:                               " << this->viscocity << std::endl;
        std::cout << "tau:                                     " << this->tau << std::endl;
        std::cout << "mass exchange factor:                    " << this->massExchangeFactor << std::endl;
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
