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

#include "CController.hpp"

#include <fstream>
#include <sstream>
#include <typeinfo>
#include <sys/time.h>

#include <mpi.h>

#include "libvis/CLbmVisualizationNetCDF.hpp"
#include "libvis/CLbmVisualizationVTK.hpp"

template <class T>
CController<T>::CController(
        int id,
        CDomain<T> domain,
        std::vector<Flag> boundaryConditions,
        std::vector<CComm<T> > communication,
        CConfiguration<T>* configuration) :
        id(id),
        configuration(configuration),
        domain(domain),
        boundaryConditions(boundaryConditions),
        communication(communication),
        simulationStepCounter(0)
{
	CDomain<T> domainGPU = decomposeSubdomain();

	solverGPU = new CLbmSolverGPU<T>(
            this->id,
            this->configuration->threadsPerBlock,
            // &domainGPU,
            this->domain,
            this->boundaryConditions,
            this->configuration->timestep,
            this->configuration->gravitation,
            this->configuration->cavityVelocity,
            this->configuration->viscosity,
            this->configuration->massExchangeFactor,
            this->configuration->maxGravitationDimLess,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doLogging);
    solverCPU = new CLbmSolverCPU<T>(
            this->id,
            this->domain,
            this->boundaryConditions,
            solverGPU,
            this->configuration->timestep,
            this->configuration->gravitation,
            this->configuration->cavityVelocity,
            this->configuration->viscosity,
            this->configuration->massExchangeFactor,
            this->configuration->maxGravitationDimLess,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doLogging);

    if (this->configuration->doVisualization)
#ifdef PAR_NETCDF
        visualization = new CLbmVisualizationNetCDF<T>(
        		id,
        		this->configuration->visualizationRate,
        		getSolver(),
        		this->configuration->numOfSubdomains,
        		this->configuration->visualizationOutputDir);
#else
        visualization = new CLbmVisualizationVTK<T>(
        		id,
        		this->configuration->visualizationRate,
        		getSolver(),
        		this->configuration->visualizationOutputDir);
#endif
}

template <class T>
CController<T>::~CController()
{
    if (configuration->doVisualization)
        delete visualization;

    delete solverCPU;
    delete solverGPU;
}

template<class T>
CDomain<T> CController<T>::decomposeSubdomain()
{
	CVector<3, int> sizeGPU(
			configuration->CPUSubdomainRatio[0] * domain.getSize()[0],
			configuration->CPUSubdomainRatio[1] * domain.getSize()[1],
			configuration->CPUSubdomainRatio[2] * domain.getSize()[2]);
	CVector<3, int> originGPU(domain.getOrigin() + ((domain.getSize() - sizeGPU) / 2));
	CVector<3, T> lengthGPU(
			domain.getLength()[0] * (T)sizeGPU[0] / (T)domain.getSize()[0],
			domain.getLength()[1] * (T)sizeGPU[1] / (T)domain.getSize()[1],
			domain.getLength()[2] * (T)sizeGPU[2] / (T)domain.getSize()[2]);

	CDomain<T> domainGPU(id, sizeGPU, originGPU, lengthGPU);

    if (configuration->doLogging)
    {
        std::cout << "----- CController<T>::decomposeSubdomain() -----" << std::endl;
        std::cout << "id:                   " << id << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "GPU subdomain size:   " << domainGPU.getSize() << std::endl;
        std::cout << "GPU subdomain length: " << domainGPU.getLength() << std::endl;
        std::cout << "GPU subdomain origin: " << domainGPU.getOrigin() << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << "CPU subdomain size:   " << domain.getSize() << std::endl;
        std::cout << "CPU subdomain length: " << domain.getLength() << std::endl;
        std::cout << "CPU subdomain origin: " << domain.getOrigin() << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
    }

	return domainGPU;
}


template <class T>
void CController<T>::syncAlpha()
{
    CVector<3, int> sendSize, recvSize;
    CVector<3, int> sendOrigin, recvOrigin;
    int dstId;
    MPI_Request request[2];
    MPI_Status status[2];
    int sendBufferSize, recvBufferSize;
    T* sendBuffer;
    T* recvBuffer;

    if (configuration->doLogging)
    {
        std::cout << "----- CController<T>::syncAlpha() -----" << std::endl;
        std::cout << "id:                  " << id << std::endl;
        std::cout << "---------------------------------------" << std::endl;
    }

    for (typename std::vector<CComm<T> >::iterator it = communication.begin(); it != communication.end(); it++)
    {
        sendSize = it->getSendSize();
        recvSize = it->getRecvSize();
        sendOrigin = it->getSendOrigin();
        recvOrigin = it->getRecvOrigin();
        dstId = it->getDstId();
        sendBufferSize = NUM_LATTICE_VECTORS * sendSize.elements();
        recvBufferSize = NUM_LATTICE_VECTORS * recvSize.elements();
        sendBuffer = new T[sendBufferSize];
        recvBuffer = new T[recvBufferSize];

        if (configuration->doLogging)
        {
            std::cout << "destination rank:    " << dstId << std::endl;
            std::cout << "send origin:         " << sendOrigin << std::endl;
            std::cout << "receive origin:      " << recvOrigin << std::endl;
            std::cout << "send buffer size:    " << ((T)sendBufferSize / (T)(1<<20)) << " MBytes" << std::endl;
            std::cout << "receive buffer size: " << ((T)recvBufferSize / (T)(1<<20)) << " MBytes" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }

        solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffer);

        MPI_Isend(
                sendBuffer,
                sendBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                MPI_TAG_ALPHA_SYNC,
                MPI_COMM_WORLD,
                &request[0]);
        MPI_Irecv(
                recvBuffer,
                recvBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                MPI_TAG_ALPHA_SYNC,
                MPI_COMM_WORLD,
                &request[1]);
        MPI_Waitall(2, request, status);

        if (configuration->doLogging)
        {
            std::cout << "Alpha synchronization successful." << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }

        solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffer);

        delete[] recvBuffer;
        delete[] sendBuffer;
    }
}

template <class T>
void CController<T>::syncBeta()
{
    CVector<3, int> sendSize, recvSize;
    CVector<3, int> sendOrigin, recvOrigin;
    int dstId;
    MPI_Request request[2];
    MPI_Status status[2];
    int sendBufferSize, recvBufferSize;
    T* sendBuffer;
    T* recvBuffer;

    if (configuration->doLogging)
    {
        std::cout << "----- CController<T>::syncBeta() -----" << std::endl;
        std::cout << "id:                  " << id << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }

    for (typename std::vector<CComm<T> >::iterator it = communication.begin(); it != communication.end(); it++)
    {
        sendSize = it->getRecvSize();
        recvSize = it->getSendSize();
        sendOrigin = it->getRecvOrigin();
        recvOrigin = it->getSendOrigin();
        dstId = it->getDstId();
        sendBufferSize = NUM_LATTICE_VECTORS * sendSize.elements();
        recvBufferSize = NUM_LATTICE_VECTORS * recvSize.elements();
        sendBuffer = new T[sendBufferSize];
        recvBuffer = new T[recvBufferSize];

        if (configuration->doLogging)
        {
            std::cout << "destination rank:    " << dstId << std::endl;
            std::cout << "send origin:         " << sendOrigin << std::endl;
            std::cout << "receive origin:      " << recvOrigin << std::endl;
            std::cout << "send buffer size:    " << ((T)sendBufferSize / (T)(1<<20)) << " MBytes" << std::endl;
            std::cout << "receive buffer size: " << ((T)recvBufferSize / (T)(1<<20)) << " MBytes" << std::endl;
            std::cout << "--------------------------------------" << std::endl;
        }

        solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffer);

        MPI_Isend(
                sendBuffer,
                sendBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                MPI_TAG_BETA_SYNC,
                MPI_COMM_WORLD,
                &request[0]);
        MPI_Irecv(
                recvBuffer,
                recvBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                MPI_TAG_BETA_SYNC,
                MPI_COMM_WORLD,
                &request[1]);
        MPI_Waitall(2, request, status);

        if (configuration->doLogging)
        {
            std::cout << "Beta synchronization successful." << std::endl;
            std::cout << "--------------------------------------" << std::endl;
        }

        solverGPU->setDensityDistributions(recvOrigin, recvSize, it->getDirection(), recvBuffer);

        delete[] recvBuffer;
        delete[] sendBuffer;
    }
}

template <class T>
void CController<T>::computeNextStep()
{
    if (simulationStepCounter & 1) {
        solverGPU->simulationStepBeta();
        syncBeta();
    } else {
        solverGPU->simulationStepAlpha();
        syncAlpha();
    }
    simulationStepCounter++;
}

template <class T>
void CController<T>::setDrivenCavitySzenario()
{
    if (configuration->doLogging)
    {
        std::cout << "----- CController<T>::setDrivenCavitySzenario() -----" << std::endl;
        std::cout << "Cells where velocity injection takes place due to driven cavity scenario are marked accordingly." << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "id:     " << id << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    CVector<3, int> origin(1, domain.getSize()[1] - 2, 1);
    CVector<3, int> size(domain.getSize()[0] - 2, 1, domain.getSize()[2] - 2);

    if (configuration->doLogging)
    {
        std::cout << "origin: " << origin << std::endl;
        std::cout << "size:   " << size << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    Flag* src = new Flag[size.elements()];

    for (int i = 0; i < size.elements(); i++)
    {
        src[i] = VELOCITY_INJECTION;
    }

    solverGPU->setFlags(origin, size, src);

    delete[] src;
}

/*
 * This function starts the simulation for the particular subdomain corresponded to
 * this class.
 */
template <class T>
void CController<T>::run()
{
    int usedDataSize;
    timeval start, end;
    T elapsed;

    usedDataSize = 0;

    if (configuration->doLogging)
    {
        std::cout << "----- CController<T>::run() -----" << std::endl;
        std::cout << "id:     " << id << std::endl;
        std::cout << "---------------------------------" << std::endl;
    }

    if (configuration->doBenchmark)
    {
        /*
         * Density distributions are read and written
         */
        usedDataSize += NUM_LATTICE_VECTORS * 2 * sizeof(T);
        /*
         * Flag is read
         */
        usedDataSize += 1 * sizeof(Flag);


        if (configuration->doVisualization || configuration->doValidation)
            /*
             * Velocity (3) and density (1) are written
             */
            usedDataSize += 4 * sizeof(T);
    }

    if (configuration->doBenchmark)
        gettimeofday(&start, NULL);

    for (int i = 0; i < configuration->loops || configuration->loops < 0; i++) {
        if (configuration->doLogging)
            std::cout << "Do iteration " << i << ":" << std::endl;

        computeNextStep();

        if (configuration->doVisualization)
            visualization->render(i);
    }

    if (configuration->doBenchmark)
    {
        gettimeofday(&end, NULL);
        elapsed = (T)(end.tv_sec - start.tv_sec) + (T)(end.tv_usec - start.tv_usec) * (T)0.000001;

        T iterationsPerSecond = (T)(configuration->loops) / elapsed;
        T glups = iterationsPerSecond * (T)configuration->domainSize.elements() * (T)0.000000001;
        T gbandwidth = glups * (T)usedDataSize;

        std::stringstream benchmarkFileName;
        benchmarkFileName << configuration->benchmarkOutputDir << "/benchmark_" << id << ".txt";
        std::ofstream benchmarkFile(benchmarkFileName.str().c_str(), std::ios::out);
        if (benchmarkFile.is_open())
        {
            benchmarkFile << "loops:           " << configuration->loops << std::endl;
            benchmarkFile << "time:            " << elapsed << "s" << std::endl;
            benchmarkFile << "iterations:      " << iterationsPerSecond << "s^-1" << std::endl;
            benchmarkFile << "lattice updates: " << glups << "GLUPS" << std::endl;
            benchmarkFile << "bandwidth:       " << gbandwidth << "GB/s" << std::endl;
            benchmarkFile.close();
        } else {
            std::cerr << "----- CController<T>::run() -----" << std::endl;
            std::cerr << "There is no open file to write benchmark results." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    if (configuration->doLogging)
    {
        std::cout << "---------------------------------" << std::endl;
    }
}

template <class T>
int CController<T>::getId() {
    return id;
}

template <class T>
CDomain<T>* CController<T>::getDomain() {
    return &domain;
}

template <class T>
CLbmSolver<T>* CController<T>::getSolver() {
    // return solverGPU;
    return solverGPU;
}

template class CController<double>;
template class CController<float>;
