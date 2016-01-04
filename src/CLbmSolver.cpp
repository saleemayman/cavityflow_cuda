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

#include <fstream>
#include <limits>
#include <sstream>

#include "libmath/CMath.hpp"

template <class T>
CLbmSolver<T>::CLbmSolver(
        int id,
        CDomain<T> &domain,
        std::vector<Flag> boundaryConditions,
        CConfiguration<T>* configuration) :
        id(id),
        domain(domain),
        boundaryConditions(boundaryConditions),
        configuration(configuration),
        timestepSize(configuration->timestep),
        viscosity(configuration->viscosity),
        storeDensities(configuration->doValidation || configuration->doVisualization),
        storeVelocities(configuration->doValidation || configuration->doVisualization)
{
    bool limitedByVelocity = false;
    bool limitedByAcceleration = false;
    T cellLength = domain.getLength()[0] / (T)domain.getSize()[0];
    T oldTimestepSize;

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
        	loggingFile << "id:                                      " << this->id << std::endl;
        	loggingFile << "---------------------------------------" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    /*
     * If no viscosity was passed to this program, a default viscosity
     * basing on a predefined reynolds number is set. This artificial
     * viscosity can be further adapted if the timestep size has to adapted.
     * The predefined reynolds number is always kept and stays constant.
     */
    if (viscosity <= (T)0)
    {
        if (configuration->doLogging)
        {
            std::stringstream loggingFileName;
            loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
            std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
            if (loggingFile.is_open())
            {
            	loggingFile << "No viscosity has been passed!" << std::endl;
            	loggingFile << "Artificial viscosity is set!" << std::endl;
            	loggingFile.close();
            } else {
                std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                std::cerr << "There is no open file to write logs." << std::endl;
                std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                std::cerr << "---------------------------------------" << std::endl;

                exit (EXIT_FAILURE);
            }
        }

        // definition reynolds number
        viscosity = configuration->domainLength.max() * configuration->velocity[0] / REYNOLDS_DEFAULT;

        if (configuration->doLogging)
        {
            std::stringstream loggingFileName;
            loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
            std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
            if (loggingFile.is_open())
            {
            	loggingFile << "viscosity:         " << viscosity << std::endl;
            	loggingFile << "---------------------------------------" << std::endl;
            	loggingFile.close();
            } else {
                std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                std::cerr << "There is no open file to write logs." << std::endl;
                std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                std::cerr << "---------------------------------------" << std::endl;

                exit (EXIT_FAILURE);
            }
        }
    }

    while (true)
    {
        if (configuration->doLogging)
        {
            std::stringstream loggingFileName;
            loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
            std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
            if (loggingFile.is_open())
            {
            	loggingFile << "New iteration of finding timestep size started." << std::endl;
            	loggingFile << "---------------------------------------" << std::endl;
            	loggingFile.close();
            } else {
                std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                std::cerr << "There is no open file to write logs." << std::endl;
                std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                std::cerr << "---------------------------------------" << std::endl;

                exit (EXIT_FAILURE);
            }
        }

        // (4.11)
        velocityDimLess = configuration->velocity * (timestepSize / cellLength);
        accelerationDimLess = configuration->acceleration * ((timestepSize * timestepSize) / cellLength);
        // (4.9)
        viscosityDimLess = viscosity * (timestepSize / (cellLength * cellLength));

        /*
         * If the dimension less velocity is larger than the specified maximum
         * value, the simulation becomes unstable. In such a case, the timestep
         * size is adapted accordingly and a valid dimension less velocity is
         * set in the next iteration of this loop.
         */
        if (velocityDimLess.length() > configuration->maxVelocityDimLess + std::numeric_limits<T>::epsilon())
        {
            limitedByVelocity= true;
            oldTimestepSize = timestepSize;
            timestepSize = (configuration->maxVelocityDimLess * cellLength) / configuration->velocity.length();

            if (configuration->doLogging)
            {
                std::stringstream loggingFileName;
                loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
                std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
                if (loggingFile.is_open())
                {
                	loggingFile << "Velocity (dimension less) is too large so the simulation could get unstable!" << std::endl;
                	loggingFile << "Timestep size is adapted!" << std::endl;
                	loggingFile << "old timestep size: " << oldTimestepSize << std::endl;
                	loggingFile << "new timestep size: " << this->timestepSize << std::endl;
                	loggingFile << "---------------------------------------" << std::endl;
                	loggingFile.close();
                } else {
                    std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                    std::cerr << "There is no open file to write logs." << std::endl;
                    std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                    std::cerr << "---------------------------------------" << std::endl;

                    exit (EXIT_FAILURE);
                }
            }

            continue;
        }

        /*
         * If the dimension less acceleration is larger than the specified
         * maximum value, the simulation becomes unstable. In such a case, the
         * timestep size is adapted accordingly and a valid dimension less
         * acceleration is set in the next iteration of this loop.
         */
        if (accelerationDimLess.length() > configuration->maxAccelerationDimLess + std::numeric_limits<T>::epsilon())
        {
            limitedByAcceleration = true;
            oldTimestepSize = timestepSize;
            // (4.12)
            timestepSize = CMath<T>::sqrt((configuration->maxAccelerationDimLess * cellLength) / configuration->acceleration.length());

            if (configuration->doLogging)
            {
                std::stringstream loggingFileName;
                loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
                std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
                if (loggingFile.is_open())
                {
                	loggingFile << "Acceleration (dimension less) is too large so the simulation could get unstable!" << std::endl;
                	loggingFile << "Timestep size is adapted!" << std::endl;
                	loggingFile << "old timestep size: " << oldTimestepSize << std::endl;
                	loggingFile << "new timestep size: " << timestepSize << std::endl;
                	loggingFile << "---------------------------------------" << std::endl;
                	loggingFile.close();
                } else {
                    std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                    std::cerr << "There is no open file to write logs." << std::endl;
                    std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                    std::cerr << "---------------------------------------" << std::endl;

                    exit (EXIT_FAILURE);
                }
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
            oldTimestepSize = timestepSize;
            // 4.9 & 4.10
            timestepSize = (cellLength * cellLength) * (((T)2 * TAU_DEFAULT - (T)1) / ((T)6 * viscosity));

            if (configuration->doLogging)
            {
                std::stringstream loggingFileName;
                loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
                std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
                if (loggingFile.is_open())
                {
                	loggingFile << "Tau " << tau << " not within the range [" << TAU_LOWER_LIMIT <<"; " << TAU_UPPER_LIMIT << "]." << std::endl;
                	loggingFile << "Timestep size is adapted!" << std::endl;
                	loggingFile << "old timestep size: " << oldTimestepSize << std::endl;
                	loggingFile << "new timestep size: " << timestepSize << std::endl;
                	loggingFile << "---------------------------------------" << std::endl;
                	loggingFile.close();
                } else {
                    std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                    std::cerr << "There is no open file to write logs." << std::endl;
                    std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                    std::cerr << "---------------------------------------" << std::endl;

                    exit (EXIT_FAILURE);
                }
            }

            if ((limitedByVelocity || limitedByAcceleration) && timestepSize > oldTimestepSize) {
                std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
                std::cerr << "No valid timestep size could be determined which satisfies" << std::endl;
                std::cerr << "- viscosity:                         " << viscosity << std::endl;
                std::cerr << "- max velocity (dimension less):     " << configuration->maxVelocityDimLess << std::endl;
                std::cerr << "- max acceleration (dimension less): " << configuration->maxAccelerationDimLess << std::endl;
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

    int reynolds = configuration->domainLength.max() * configuration->velocity[0] / viscosity;

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "global length (without halo):      " << configuration->domainLength << std::endl;
        	loggingFile << "---------------------------------------" << std::endl;
        	loggingFile << "domain size (without halo):        " << this->domain.getSize() << std::endl;
        	loggingFile << "domain length (without halo):      " << this->domain.getLength() << std::endl;
            loggingFile << "domain origin (without halo):      " << this->domain.getOrigin() << std::endl;
            loggingFile << "---------------------------------------" << std::endl;
            loggingFile << "timestep size:                     " << timestepSize << std::endl;
            loggingFile << "---------------------------------------" << std::endl;
            loggingFile << "velocity:                          " << configuration->velocity << std::endl;
            loggingFile << "velocity (dimension less):         " << velocityDimLess << std::endl;
            loggingFile << "acceleration:                      " << configuration->acceleration << std::endl;
            loggingFile << "acceleration (dimension less):     " << accelerationDimLess << std::endl;
            loggingFile << "---------------------------------------" << std::endl;
            loggingFile << "viscosity:                         " << viscosity << std::endl;
            loggingFile << "viscosity (dimension less):        " << viscosityDimLess << std::endl;
            loggingFile << "tau:                               " << tau << std::endl;
            loggingFile << "reynolds number (dimension less):  " << reynolds << std::endl;
            loggingFile << "max velocity (dimension less):     " << configuration->maxVelocityDimLess << std::endl;
            loggingFile << "max acceleration (dimension less): " << configuration->maxAccelerationDimLess << std::endl;
            loggingFile << "inv tau:                           " << tauInv << std::endl;
            loggingFile << "---------------------------------------" << std::endl;
            loggingFile << "store densities:                   " << storeDensities << std::endl;
            loggingFile << "store velocities:                  " << storeVelocities << std::endl;
            loggingFile << "---------------------------------------" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolver<T>::CLbmSolver() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }
}

template <class T>
CDomain<T>* CLbmSolver<T>::getDomain()
{
    return &domain;
}

template class CLbmSolver<float>;
template class CLbmSolver<double>;
