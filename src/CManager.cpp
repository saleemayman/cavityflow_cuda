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
        std::cout << "----- CManager<T>::CManager() -----" << std::endl;
        std::cout << "id:                     " << rank << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "domain size:            " << domain.getSize() << std::endl;
        std::cout << "domain length:          " << domain.getLength() << std::endl;
        std::cout << "domain origin:          " << domain.getOrigin() << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "subdomain coordinates:  [" << subdomainX << ", " << subdomainY << ", " << subdomainZ << "]" << std::endl;
        std::cout << "subdomain size:         " << subdomainSize << std::endl;
        std::cout << "subdomain length:       " << subdomainLength << std::endl;
        std::cout << "subdomain origin:       " << subdomainOrigin << std::endl;
        std::cout << "-----------------------------------" << std::endl;
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
        std::cout << "boundaryConditions[0 = LEFT]:   " << boundaryConditions[0] << std::endl;
        std::cout << "boundaryConditions[1 = RIGHT]:  " << boundaryConditions[1] << std::endl;
        std::cout << "boundaryConditions[2 = BOTTOM]: " << boundaryConditions[2] << std::endl;
        std::cout << "boundaryConditions[3 = TOP]:    " << boundaryConditions[3] << std::endl;
        std::cout << "boundaryConditions[4 = BACK]:   " << boundaryConditions[4] << std::endl;
        std::cout << "boundaryConditions[5 = FRONT]:  " << boundaryConditions[5] << std::endl;
        std::cout << "-----------------------------------" << std::endl;
    }

    /*
     * Initializing the Controller's communication classes based on the already
     * computed boundary conditions
     */
    std::vector<CComm<T> > communication;

    if (boundaryConditions[0] == GHOST_LAYER) {
        int commDestination = rank - 1;
        CVector<3, int> sendSize(1, subdomainSize[1], subdomainSize[2]);
        CVector<3, int> recvSize(1, subdomainSize[1], subdomainSize[2]);
        CVector<3, int> sendOrigin(1, 0, 0);
        CVector<3, int> recvOrigin(0, 0, 0);
        CVector<3, int> commDirection(-1, 0, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, commDirection);

        communication.push_back(comm);
    }
    if (boundaryConditions[1] == GHOST_LAYER) {
        int commDestination = rank + 1;
        CVector<3, int> sendSize(1, subdomainSize[1], subdomainSize[2]);
        CVector<3, int> recvSize(1, subdomainSize[1], subdomainSize[2]);
        CVector<3, int> sendOrigin(subdomainSize[0] - 2, 0, 0);
        CVector<3, int> recvOrigin(subdomainSize[0] - 1, 0, 0);
        CVector<3, int> commDirection(1, 0, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, commDirection);

        communication.push_back(comm);
    }
    if (boundaryConditions[2] == GHOST_LAYER) {
        int commDestination = rank - configuration->numOfSubdomains[0];
        CVector<3, int> sendSize(subdomainSize[0], 1, subdomainSize[2]);
        CVector<3, int> recvSize(subdomainSize[0], 1, subdomainSize[2]);
        CVector<3, int> sendOrigin(0, 1, 0);
        CVector<3, int> recvOrigin(0, 0, 0);
        CVector<3, int> commDirection(0, -1, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, commDirection);

        communication.push_back(comm);
    }
    if (boundaryConditions[3] == GHOST_LAYER) {
        int commDestination = rank + configuration->numOfSubdomains[0];
        CVector<3, int> sendSize(subdomainSize[0], 1, subdomainSize[2]);
        CVector<3, int> recvSize(subdomainSize[0], 1, subdomainSize[2]);
        CVector<3, int> sendOrigin(0, subdomainSize[1] - 2, 0);
        CVector<3, int> recvOrigin(0, subdomainSize[1] - 1, 0);
        CVector<3, int> commDirection(0, 1, 0);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, commDirection);

        communication.push_back(comm);
    }
    if (boundaryConditions[4] == GHOST_LAYER) {
        int commDestination = rank - configuration->numOfSubdomains[0] * configuration->numOfSubdomains[1];
        CVector<3, int> sendSize(subdomainSize[0], subdomainSize[1], 1);
        CVector<3, int> recvSize(subdomainSize[0], subdomainSize[1], 1);
        CVector<3, int> sendOrigin(0, 0, 1);
        CVector<3, int> recvOrigin(0, 0, 0);
        CVector<3, int> commDirection(0, 0, -1);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, commDirection);

        communication.push_back(comm);
    }
    if (boundaryConditions[5] == GHOST_LAYER) {
        int commDestination = rank + configuration->numOfSubdomains[0] * configuration->numOfSubdomains[1];
        CVector<3, int> sendSize(subdomainSize[0], subdomainSize[1], 1);
        CVector<3, int> recvSize(subdomainSize[0], subdomainSize[1], 1);
        CVector<3, int> sendOrigin(0, 0, subdomainSize[2] - 2);
        CVector<3, int> recvOrigin(0, 0, subdomainSize[2] - 1);
        CVector<3, int> commDirection(0, 0, 1);

        CComm<T> comm(commDestination, sendSize, recvSize, sendOrigin, recvOrigin, commDirection);

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
