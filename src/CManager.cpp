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

#include "CManager.hpp"

#include <fstream>
#include <sstream>

#include "libmath/CVector.hpp"

template <class T>
CManager<T>::CManager(int rank, CConfiguration<T>* configuration) :
        domain(rank, configuration->domainSize, CVector<3, int>(0), configuration->domainLength)
{
    /*
     * Determine parameters for subdomain managed by this class.
     */
    int id = rank;
    int subdomainX, subdomainY, subdomainZ;

    subdomainX = id % configuration->numOfSubdomains[0];
    id /= configuration->numOfSubdomains[0];
    subdomainY = id % configuration->numOfSubdomains[1];
    id /= configuration->numOfSubdomains[1];
    subdomainZ = id;

    /*
     * Create the subdomain instances for the whole domain
     */
    CVector<3, int> subdomainSize(
            configuration->domainSize[0] / configuration->numOfSubdomains[0],
            configuration->domainSize[1] / configuration->numOfSubdomains[1],
            configuration->domainSize[2] / configuration->numOfSubdomains[2]);
    CVector<3, int> subdomainOrigin(
            subdomainX * subdomainSize[0],
            subdomainY * subdomainSize[1],
            subdomainZ * subdomainSize[2]);
    CVector<3, T> subdomainLength(
            configuration->domainLength[0] / (T)configuration->numOfSubdomains[0],
            configuration->domainLength[1] / (T)configuration->numOfSubdomains[1],
            configuration->domainLength[2] / (T)configuration->numOfSubdomains[2]);
    CDomain<T> *subdomain = new CDomain<T>(rank, subdomainSize, subdomainOrigin, subdomainLength);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << rank << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "----- CManager<T>::CManager() -----" << std::endl;
        	loggingFile << "domain size:            " << domain.getSize() << std::endl;
        	loggingFile << "domain length:          " << domain.getLength() << std::endl;
        	loggingFile << "domain origin:          " << domain.getOrigin() << std::endl;
        	loggingFile << "-----------------------------------" << std::endl;
        	loggingFile << "subdomain coordinates:  [" << subdomainX << ", " << subdomainY << ", " << subdomainZ << "]" << std::endl;
        	loggingFile << "subdomain size:         " << subdomainSize << std::endl;
        	loggingFile << "subdomain length:       " << subdomainLength << std::endl;
        	loggingFile << "subdomain origin:       " << subdomainOrigin << std::endl;
        	loggingFile << "-----------------------------------" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CManager<T>::CManager() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-----------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    /*
     * Setting the boundary conditions for the current Controller
     */
    std::vector<Flag> boundaryConditions(6, GHOST_LAYER);

    if (subdomainX == 0)
        boundaryConditions[0] = OBSTACLE;
    if (subdomainX == (configuration->numOfSubdomains[0] - 1))
        boundaryConditions[1] = OBSTACLE;
    if (subdomainY == 0)
        boundaryConditions[2] = OBSTACLE;
    if (subdomainY == (configuration->numOfSubdomains[1] - 1))
        boundaryConditions[3] = OBSTACLE;
    if (subdomainZ == 0)
        boundaryConditions[4] = OBSTACLE;
    if (subdomainZ == (configuration->numOfSubdomains[2] - 1))
        boundaryConditions[5] = OBSTACLE;

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << rank << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "boundaryConditions[0 = LEFT]:   " << boundaryConditions[0] << std::endl;
        	loggingFile << "boundaryConditions[1 = RIGHT]:  " << boundaryConditions[1] << std::endl;
        	loggingFile << "boundaryConditions[2 = BOTTOM]: " << boundaryConditions[2] << std::endl;
        	loggingFile << "boundaryConditions[3 = TOP]:    " << boundaryConditions[3] << std::endl;
        	loggingFile << "boundaryConditions[4 = BACK]:   " << boundaryConditions[4] << std::endl;
        	loggingFile << "boundaryConditions[5 = FRONT]:  " << boundaryConditions[5] << std::endl;
        	loggingFile << "-----------------------------------" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CManager<T>::CManager() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-----------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    /*
     * Initializing the Controller's communication classes based on the already
     * computed boundary conditions
     */
    std::vector<CComm<T> > communication;

    if (boundaryConditions[0] == GHOST_LAYER) {
        int commDestination = rank - 1;
        CVector<3, int> sendSize(1, subdomain->getSizeWithHalo()[1], subdomain->getSizeWithHalo()[2]);
        CVector<3, int> recvSize(1, subdomain->getSizeWithHalo()[1], subdomain->getSizeWithHalo()[2]);
        CVector<3, int> sendOrigin(1, 0, 0);
        CVector<3, int> recvOrigin(0, 0, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, LEFT);

        communication.push_back(comm);
    }
    if (boundaryConditions[1] == GHOST_LAYER) {
        int commDestination = rank + 1;
        CVector<3, int> sendSize(1, subdomain->getSizeWithHalo()[1], subdomain->getSizeWithHalo()[2]);
        CVector<3, int> recvSize(1, subdomain->getSizeWithHalo()[1], subdomain->getSizeWithHalo()[2]);
        CVector<3, int> sendOrigin(subdomain->getSizeWithHalo()[0] - 2, 0, 0);
        CVector<3, int> recvOrigin(subdomain->getSizeWithHalo()[0] - 1, 0, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, RIGHT);

        communication.push_back(comm);
    }
    if (boundaryConditions[2] == GHOST_LAYER) {
        int commDestination = rank - configuration->numOfSubdomains[0];
        CVector<3, int> sendSize(subdomain->getSizeWithHalo()[0], 1, subdomain->getSizeWithHalo()[2]);
        CVector<3, int> recvSize(subdomain->getSizeWithHalo()[0], 1, subdomain->getSizeWithHalo()[2]);
        CVector<3, int> sendOrigin(0, 1, 0);
        CVector<3, int> recvOrigin(0, 0, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, BOTTOM);

        communication.push_back(comm);
    }
    if (boundaryConditions[3] == GHOST_LAYER) {
        int commDestination = rank + configuration->numOfSubdomains[0];
        CVector<3, int> sendSize(subdomain->getSizeWithHalo()[0], 1, subdomain->getSizeWithHalo()[2]);
        CVector<3, int> recvSize(subdomain->getSizeWithHalo()[0], 1, subdomain->getSizeWithHalo()[2]);
        CVector<3, int> sendOrigin(0, subdomain->getSizeWithHalo()[1] - 2, 0);
        CVector<3, int> recvOrigin(0, subdomain->getSizeWithHalo()[1] - 1, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, TOP);

        communication.push_back(comm);
    }
    if (boundaryConditions[4] == GHOST_LAYER) {
        int commDestination = rank - configuration->numOfSubdomains[0] * configuration->numOfSubdomains[1];
        CVector<3, int> sendSize(subdomain->getSizeWithHalo()[0], subdomain->getSizeWithHalo()[1], 1);
        CVector<3, int> recvSize(subdomain->getSizeWithHalo()[0], subdomain->getSizeWithHalo()[1], 1);
        CVector<3, int> sendOrigin(0, 0, 1);
        CVector<3, int> recvOrigin(0, 0, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, BACK);

        communication.push_back(comm);
    }
    if (boundaryConditions[5] == GHOST_LAYER) {
        int commDestination = rank + configuration->numOfSubdomains[0] * configuration->numOfSubdomains[1];
        CVector<3, int> sendSize(subdomain->getSizeWithHalo()[0], subdomain->getSizeWithHalo()[1], 1);
        CVector<3, int> recvSize(subdomain->getSizeWithHalo()[0], subdomain->getSizeWithHalo()[1], 1);
        CVector<3, int> sendOrigin(0, 0, subdomain->getSizeWithHalo()[2] - 2);
        CVector<3, int> recvOrigin(0, 0, subdomain->getSizeWithHalo()[2] - 1);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, FRONT);

        communication.push_back(comm);
    }

    controller = new CController<T>(rank, *subdomain, boundaryConditions, communication, configuration);

    if (subdomainY == configuration->numOfSubdomains[1] - 1) {
        controller->setDrivenCavitySzenario();
    }
}

template <class T>
CManager<T>::~CManager()
{
    delete controller;
}

template <class T>
void CManager<T>::run()
{
    controller->run();
}

template <class T>
CDomain<T>* CManager<T>::getDomain()
{
    return &domain;
}

template <class T>
CDomain<T>* CManager<T>::getSubDomain()
{
    return controller->getDomain();
}

template <class T>
CController<T>* CManager<T>::getController()
{
    return controller;
}

template class CManager<double>;
template class CManager<float>;
