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

#include <cassert>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <sys/time.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "libvis/CLbmVisualizationNetCDF.hpp"

template <class T>
CController<T>::CController(
        int id,
        CDomain<T> domain,
        std::vector<Flag> boundaryConditions,
        std::vector<CComm<T> > communication,
        CConfiguration<T>* configuration) :
        id(id),
        domain(domain),
        boundaryConditions(boundaryConditions),
        communication(communication),
        configuration(configuration),
        simulationStepCounter(0)
{
    CVector<3, int> originGPU = domain.getOrigin();
    CVector<3, int> sizeCPU = domain.getSize();
    sizeCPU[1] *= configuration->CPUSubdomainRatio;
    CVector<3, int> originCPU = domain.getOrigin();
    originGPU[1] += sizeCPU[1];
    CVector<3, int> sizeGPU = domain.getSize();
    sizeGPU[1] -= sizeCPU[1];
    CVector<3, T> lengthCPU = domain.getLength();
    lengthCPU[1] = domain.getLength()[1] * (T)sizeCPU[1] / (T)domain.getSize()[1];
    CVector<3, T> lengthGPU = domain.getLength();
    lengthGPU[1] -= lengthCPU[1];

    CDomain<T> domainCPU(id, sizeCPU, originCPU, lengthCPU);
    CDomain<T> domainGPU(id, sizeGPU, originGPU, lengthGPU);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "----- CController<T>::CController() -----" << std::endl;
        	loggingFile << "CPU subdomain size:   " << domainCPU.getSize() << std::endl;
        	loggingFile << "CPU subdomain length: " << domainCPU.getLength() << std::endl;
        	loggingFile << "CPU subdomain origin: " << domainCPU.getOrigin() << std::endl;
        	loggingFile << "-----------------------------------------" << std::endl;
        	loggingFile << "GPU subdomain size:   " << domainGPU.getSize() << std::endl;
        	loggingFile << "GPU subdomain length: " << domainGPU.getLength() << std::endl;
        	loggingFile << "GPU subdomain origin: " << domainGPU.getOrigin() << std::endl;
        	loggingFile << "-----------------------------------------" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CController<T>::CController() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-----------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    solverGPU = new CLbmSolverGPU<T>(
            this->id,
            domainGPU,
            this->boundaryConditions,
            this->configuration);
    solverCPU = new CLbmSolverCPU<T>(
            this->id,
            domainCPU,
            solverGPU,
            this->boundaryConditions,
            this->configuration);

#ifdef USE_MPI
    sendBuffers = new std::vector<T*>(this->communication.size());
    recvBuffers = new std::vector<T*>(this->communication.size());
    sendRequests = new MPI_Request[this->communication.size()];
    recvRequests = new MPI_Request[this->communication.size()];
    streams = new std::vector<cudaStream_t>(this->communication.size());
#endif
    GPU_ERROR_CHECK(cudaStreamCreate(&defaultStream))

#ifdef USE_MPI
    for (unsigned int i = 0; i < this->communication.size(); i++)
    {
        sendBuffers->at(i) = new T[NUM_LATTICE_VECTORS * this->communication[i].getSendSize().elements()];
        recvBuffers->at(i) = new T[NUM_LATTICE_VECTORS * this->communication[i].getRecvSize().elements()];
        GPU_ERROR_CHECK(cudaStreamCreate(&streams->at(i)))
    }
#endif

    if (this->configuration->doVisualization)
        /*
        visualization = new CLbmVisualizationVTK<T>(
                id,
                this->configuration->visualizationRate,
                getSolver(),
                this->configuration->visualizationOutputDir);
        */
        visualization = new CLbmVisualizationNetCDF<T>(
                id,
                this->configuration->visualizationRate,
                getSolver(),
#if defined(USE_MPI) && defined(PAR_NETCDF)
                this->configuration->numOfSubdomains,
#endif
                this->configuration->visualizationOutputDir);
}

template <class T>
CController<T>::~CController()
{
	if (configuration->doVisualization)
        delete visualization;

/*
#ifdef USE_MPI
	for (unsigned int i = communication.size() - 1; i >= 0; i--)
    {
    	GPU_ERROR_CHECK(cudaStreamDestroy(streams->at(i)))
        delete[] recvBuffers->at(i);
        delete[] sendBuffers->at(i);
    }
#endif
*/

	GPU_ERROR_CHECK(cudaStreamDestroy(defaultStream))
#ifdef USE_MPI
    delete streams;
	delete[] recvRequests;
    delete[] sendRequests;
    delete recvBuffers;
    delete sendBuffers;
#endif

    delete solverCPU;
    delete solverGPU;
}

template <class T>
void CController<T>::stepAlpha()
{
    CVector<3, int> innerOrigin(0);
    CVector<3, int> innerSize(domain.getSizeWithHalo());

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "----- CController<T>::stepAlpha() -----" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CController<T>::stepAlpha() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

#ifdef USE_MPI
    CVector<3, int> boundaryOrigin;
    CVector<3, int> boundarySize;

    for (unsigned int i = 0; i < communication.size(); i++)
    {
        MPI_Irecv(
                recvBuffers->at(i),
                NUM_LATTICE_VECTORS * communication[i].getRecvSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[i].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &recvRequests[i]);
    }

    for (unsigned int i = 0; i < communication.size(); i++)
    {
        switch(communication[i].getDirection())
        {
        case LEFT:
            boundaryOrigin = innerOrigin;
            boundarySize.set(2, innerSize[1], innerSize[2]);
            innerOrigin[0] += 2;
            innerSize[0] -= 2;
            break;
        case RIGHT:
            boundaryOrigin = innerOrigin;
            boundaryOrigin[0] = domain.getSizeWithHalo()[0] - 2;
            boundarySize.set(2, innerSize[1], innerSize[2]);
            innerSize[0] -= 2;
            break;
        case BOTTOM:
            boundaryOrigin = innerOrigin;
            boundarySize.set(innerSize[0], 2, innerSize[2]);
            innerOrigin[1] += 2;
            innerSize[1] -= 2;
            break;
        case TOP:
            boundaryOrigin = innerOrigin;
            boundaryOrigin[1] = domain.getSizeWithHalo()[1] - 2;
            boundarySize.set(innerSize[0], 2, innerSize[2]);
            innerSize[1] -= 2;
            break;
        case BACK:
            boundaryOrigin = innerOrigin;
            boundarySize.set(innerSize[0], innerSize[1], 2);
            innerOrigin[2] += 2;
            innerSize[2] -= 2;
            break;
        case FRONT:
            boundaryOrigin = innerOrigin;
            boundaryOrigin[2] = domain.getSizeWithHalo()[2] - 2;
            boundarySize.set(innerSize[0], innerSize[1], 2);
            innerSize[2] -= 2;
            break;
        }

        if (configuration->doLogging)
        {
            std::stringstream loggingFileName;
            loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
            std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
            if (loggingFile.is_open())
            {
            	loggingFile << "destination rank:       " << communication[i].getDstId() << std::endl;
            	loggingFile << "boundary origin:        " << boundaryOrigin << std::endl;
            	loggingFile << "boundary size:          " << boundarySize << std::endl;
            	loggingFile << "---------------------------------------" << std::endl;
            	loggingFile.close();
            } else {
                std::cerr << "----- CController<T>::stepAlpha() -----" << std::endl;
                std::cerr << "There is no open file to write logs." << std::endl;
                std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
                std::cerr << "---------------------------------------" << std::endl;

                exit (EXIT_FAILURE);
            }
        }

        // solverGPU->simulationStepAlpha(boundaryOrigin, boundarySize, &streams->at(i));
		#pragma omp parallel
		#pragma omp single
		{
			solverCPU->simulationStepAlpha(boundaryOrigin, boundarySize);
		}
    }
    GPU_ERROR_CHECK(cudaDeviceSynchronize())
#endif

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
        	loggingFile << "inner origin:           " << innerOrigin << std::endl;
        	loggingFile << "inner size:             " << innerSize << std::endl;
        	loggingFile << "---------------------------------------" << std::endl;
        	loggingFile.close();
        } else {
            std::cerr << "----- CController<T>::stepAlpha() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    // solverGPU->simulationStepAlpha(innerOrigin, innerSize, &defaultStream);
	#pragma omp parallel
	#pragma omp single
	{
		solverCPU->simulationStepAlpha(innerOrigin, innerSize);
	}

#ifdef USE_MPI
    int sendIdx, recvIdx;
    CVector<3, int> sendOrigin, recvOrigin;
    CVector<3, int> sendSize, recvSize;

    /*
    for (unsigned int i = 0; i < communication.size(); i++)
    {
        sendOrigin = communication[i].getSendOrigin();
        sendSize = communication[i].getSendSize();
        recvOrigin = communication[i].getRecvOrigin();
        recvSize = communication[i].getRecvSize();

        if (configuration->doLogging)
        {
			std::stringstream loggingFileName;
			loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
			std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
			if (loggingFile.is_open())
			{
				loggingFile << "destination rank:       " << communication[i].getDstId() << std::endl;
				loggingFile << "send origin:            " << sendOrigin << std::endl;
				loggingFile << "send size:              " << sendSize << std::endl;
				loggingFile << "receive origin:         " << recvOrigin << std::endl;
				loggingFile << "receive size:           " << recvSize << std::endl;
				loggingFile << "amount of send data:    " << ((T)(sendSize.elements() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "amount of receive data: " << ((T)(recvSize.elements() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "---------------------------------------" << std::endl;
				loggingFile.close();
			} else {
				std::cerr << "----- CController<T>::stepAlpha() -----" << std::endl;
				std::cerr << "There is no open file to write logs." << std::endl;
				std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
				std::cerr << "---------------------------------------" << std::endl;

				exit (EXIT_FAILURE);
			}
        }

        solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(i), &streams->at(i));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(i)))

        MPI_Isend(
                sendBuffers->at(i),
                NUM_LATTICE_VECTORS * communication[i].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[i].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[i]);

        MPI_Wait(&recvRequests[i], MPI_STATUS_IGNORE);

        solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(i), &streams->at(i));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(i)))
    }
    */

    sendIdx = recvIdx = 0;

    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == LEFT)
    {
        sendOrigin = communication[sendIdx].getSendOrigin();
        sendSize = communication[sendIdx].getSendSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == RIGHT)
    {
        sendOrigin = communication[sendIdx].getSendOrigin();
        sendSize = communication[sendIdx].getSendSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == LEFT)
    {
        recvOrigin = communication[recvIdx].getRecvOrigin();
        recvSize = communication[recvIdx].getRecvSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == RIGHT)
    {
        recvOrigin = communication[recvIdx].getRecvOrigin();
        recvSize = communication[recvIdx].getRecvSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == BOTTOM)
    {
        sendOrigin = communication[sendIdx].getSendOrigin();
        sendSize = communication[sendIdx].getSendSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == TOP)
    {
        sendOrigin = communication[sendIdx].getSendOrigin();
        sendSize = communication[sendIdx].getSendSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == BOTTOM)
    {
        recvOrigin = communication[recvIdx].getRecvOrigin();
        recvSize = communication[recvIdx].getRecvSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == TOP)
    {
        recvOrigin = communication[recvIdx].getRecvOrigin();
        recvSize = communication[recvIdx].getRecvSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == BACK)
    {
        sendOrigin = communication[sendIdx].getSendOrigin();
        sendSize = communication[sendIdx].getSendSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == FRONT)
    {
        sendOrigin = communication[sendIdx].getSendOrigin();
        sendSize = communication[sendIdx].getSendSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == BACK)
    {
        recvOrigin = communication[recvIdx].getRecvOrigin();
        recvSize = communication[recvIdx].getRecvSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == FRONT)
    {
        recvOrigin = communication[recvIdx].getRecvOrigin();
        recvSize = communication[recvIdx].getRecvSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }

    MPI_Waitall(communication.size(), sendRequests, MPI_STATUS_IGNORE);
#endif
    GPU_ERROR_CHECK(cudaStreamSynchronize(defaultStream))

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "Alpha step successful." << std::endl;
			loggingFile << "---------------------------------------" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::stepAlpha() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "---------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
    }
}

/*
template <class T>
void CController<T>::stepAlpha()
{
    solverGPU->simulationStepAlpha();
#ifdef USE_MPI
    syncAlpha();
#endif
}
*/

template <class T>
void CController<T>::stepBeta()
{
    CVector<3, int> innerOrigin(0);
    CVector<3, int> innerSize(domain.getSizeWithHalo());

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "----- CController<T>::stepBeta() -----" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::stepBeta() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "--------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
    }

#ifdef USE_MPI
    CVector<3, int> boundaryOrigin;
    CVector<3, int> boundarySize;

    for (unsigned int i = 0; i < communication.size(); i++)
    {
        MPI_Irecv(
                recvBuffers->at(i),
                NUM_LATTICE_VECTORS * communication[i].getRecvSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[i].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &recvRequests[i]);
    }

    for (unsigned int i = 0; i < communication.size(); i++)
    {
        switch(communication[i].getDirection())
        {
        case LEFT:
            boundaryOrigin = innerOrigin;
            boundarySize.set(2, innerSize[1], innerSize[2]);
            innerOrigin[0] += 2;
            innerSize[0] -= 2;
            break;
        case RIGHT:
            boundaryOrigin = innerOrigin;
            boundaryOrigin[0] = domain.getSizeWithHalo()[0] - 2;
            boundarySize.set(2, innerSize[1], innerSize[2]);
            innerSize[0] -= 2;
            break;
        case BOTTOM:
            boundaryOrigin = innerOrigin;
            boundarySize.set(innerSize[0], 2, innerSize[2]);
            innerOrigin[1] += 2;
            innerSize[1] -= 2;
            break;
        case TOP:
            boundaryOrigin = innerOrigin;
            boundaryOrigin[1] = domain.getSizeWithHalo()[1] - 2;
            boundarySize.set(innerSize[0], 2, innerSize[2]);
            innerSize[1] -= 2;
            break;
        case BACK:
            boundaryOrigin = innerOrigin;
            boundarySize.set(innerSize[0], innerSize[1], 2);
            innerOrigin[2] += 2;
            innerSize[2] -= 2;
            break;
        case FRONT:
            boundaryOrigin = innerOrigin;
            boundaryOrigin[2] = domain.getSizeWithHalo()[2] - 2;
            boundarySize.set(innerSize[0], innerSize[1], 2);
            innerSize[2] -= 2;
            break;
        }

        if (configuration->doLogging)
        {
    		std::stringstream loggingFileName;
    		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
    		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
    		if (loggingFile.is_open())
    		{
    			loggingFile << "destination rank:       " << communication[i].getDstId() << std::endl;
    			loggingFile << "boundary origin:        " << boundaryOrigin << std::endl;
    			loggingFile << "boundary size:          " << boundarySize << std::endl;
    			loggingFile << "--------------------------------------" << std::endl;
    			loggingFile.close();
    		} else {
    			std::cerr << "----- CController<T>::stepBeta() -----" << std::endl;
    			std::cerr << "There is no open file to write logs." << std::endl;
    			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
    			std::cerr << "--------------------------------------" << std::endl;

    			exit (EXIT_FAILURE);
    		}
        }

        // solverGPU->simulationStepBeta(boundaryOrigin, boundarySize, &streams->at(i));
		#pragma omp parallel
		#pragma omp single
		{
			solverCPU->simulationStepBeta(boundaryOrigin, boundarySize);
		}
    }
    GPU_ERROR_CHECK(cudaDeviceSynchronize())
#endif

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "inner origin:           " << innerOrigin << std::endl;
			loggingFile << "inner size:             " << innerSize << std::endl;
			loggingFile << "--------------------------------------" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::stepBeta() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "--------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
    }

    // solverGPU->simulationStepBeta(innerOrigin, innerSize, &defaultStream);
	#pragma omp parallel
	#pragma omp single
	{
		solverCPU->simulationStepBeta(innerOrigin, innerSize);
	}

#ifdef USE_MPI
    int sendIdx, recvIdx;
    CVector<3, int> sendOrigin, recvOrigin;
    CVector<3, int> sendSize, recvSize;

    /*
    for (unsigned int i = 0; i < communication.size(); i++)
    {
        sendOrigin = communication[i].getRecvOrigin();
        sendSize = communication[i].getRecvSize();
        recvOrigin = communication[i].getSendOrigin();
        recvSize = communication[i].getSendSize();

        if (configuration->doLogging)
        {
			std::stringstream loggingFileName;
			loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
			std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
			if (loggingFile.is_open())
			{
				loggingFile << "destination rank:       " << communication[i].getDstId() << std::endl;
				loggingFile << "send origin:            " << sendOrigin << std::endl;
				loggingFile << "send size:              " << sendSize << std::endl;
				loggingFile << "receive origin:         " << recvOrigin << std::endl;
				loggingFile << "receive size:           " << recvSize << std::endl;
				loggingFile << "amount of send data:    " << ((T)(sendSize.elements() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "amount of receive data: " << ((T)(recvSize.elements() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "--------------------------------------" << std::endl;
				loggingFile.close();
			} else {
				std::cerr << "----- CController<T>::stepBeta() -----" << std::endl;
				std::cerr << "There is no open file to write logs." << std::endl;
				std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
				std::cerr << "--------------------------------------" << std::endl;

				exit (EXIT_FAILURE);
			}
        }

        solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(i), &streams->at(i));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(i)))

        MPI_Isend(
                sendBuffers->at(i),
                NUM_LATTICE_VECTORS * communication[i].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[i].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[i]);

        MPI_Wait(&recvRequests[i], MPI_STATUS_IGNORE);

        solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[i].getDirection(), recvBuffers->at(i), &streams->at(i));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(i)))
    }
    */

    sendIdx = recvIdx = 0;

    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == LEFT)
    {
        sendOrigin = communication[sendIdx].getRecvOrigin();
        sendSize = communication[sendIdx].getRecvSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == RIGHT)
    {
        sendOrigin = communication[sendIdx].getRecvOrigin();
        sendSize = communication[sendIdx].getRecvSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == LEFT)
    {
        recvOrigin = communication[recvIdx].getSendOrigin();
        recvSize = communication[recvIdx].getSendSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == RIGHT)
    {
        recvOrigin = communication[recvIdx].getSendOrigin();
        recvSize = communication[recvIdx].getSendSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == BOTTOM)
    {
        sendOrigin = communication[sendIdx].getRecvOrigin();
        sendSize = communication[sendIdx].getRecvSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == TOP)
    {
        sendOrigin = communication[sendIdx].getRecvOrigin();
        sendSize = communication[sendIdx].getRecvSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == BOTTOM)
    {
        recvOrigin = communication[recvIdx].getSendOrigin();
        recvSize = communication[recvIdx].getSendSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == TOP)
    {
        recvOrigin = communication[recvIdx].getSendOrigin();
        recvSize = communication[recvIdx].getSendSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == BACK)
    {
        sendOrigin = communication[sendIdx].getRecvOrigin();
        sendSize = communication[sendIdx].getRecvSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (sendIdx < communication.size() && communication[sendIdx].getDirection() == FRONT)
    {
        sendOrigin = communication[sendIdx].getRecvOrigin();
        sendSize = communication[sendIdx].getRecvSize();

        // solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx), &streams->at(sendIdx));
        solverCPU->getDensityDistributions(sendOrigin, sendSize, sendBuffers->at(sendIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(sendIdx)))

        MPI_Isend(
                sendBuffers->at(sendIdx),
                NUM_LATTICE_VECTORS * communication[sendIdx].getSendSize().elements(),
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                communication[sendIdx].getDstId(),
                simulationStepCounter,
                MPI_COMM_WORLD,
                &sendRequests[sendIdx]);

        sendIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == BACK)
    {
        recvOrigin = communication[recvIdx].getSendOrigin();
        recvSize = communication[recvIdx].getSendSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }
    if (recvIdx < communication.size() && communication[recvIdx].getDirection() == FRONT)
    {
        recvOrigin = communication[recvIdx].getSendOrigin();
        recvSize = communication[recvIdx].getSendSize();

        MPI_Wait(&recvRequests[recvIdx], MPI_STATUS_IGNORE);

        // solverGPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx), &streams->at(recvIdx));
        solverCPU->setDensityDistributions(recvOrigin, recvSize, communication[recvIdx].getDirection(), recvBuffers->at(recvIdx));
        GPU_ERROR_CHECK(cudaStreamSynchronize(streams->at(recvIdx)))

        recvIdx++;
    }

    MPI_Waitall(communication.size(), sendRequests, MPI_STATUS_IGNORE);
#endif
    GPU_ERROR_CHECK(cudaStreamSynchronize(defaultStream))

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "Beta step successful." << std::endl;
			loggingFile << "--------------------------------------" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::stepBeta() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "--------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
    }
}

/*
template <class T>
void CController<T>::stepBeta()
{
    solverGPU->simulationStepBeta();
#ifdef USE_MPI
    syncBeta();
#endif
}
*/

/*
#ifdef USE_MPI
template <class T>
void CController<T>::syncAlpha()
{
    CVector<3, int> sendSize, recvSize;
    CVector<3, int> sendOrigin, recvOrigin;
    int dstId;
    MPI_Request request[2];
    int sendBufferSize, recvBufferSize;
    T* sendBuffer;
    T* recvBuffer;

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "----- CController<T>::syncAlpha() -----" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::syncAlpha() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "---------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
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
			std::stringstream loggingFileName;
			loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
			std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
			if (loggingFile.is_open())
			{
				loggingFile << "destination rank:           " << dstId << std::endl;
				loggingFile << "send origin (with halo):    " << sendOrigin << std::endl;
				loggingFile << "receive origin (with halo): " << recvOrigin << std::endl;
				loggingFile << "send buffer size:           " << ((T)(sendBufferSize * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "receive buffer size:        " << ((T)(recvBufferSize * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "---------------------------------------" << std::endl;
				loggingFile.close();
			} else {
				std::cerr << "----- CController<T>::syncAlpha() -----" << std::endl;
				std::cerr << "There is no open file to write logs." << std::endl;
				std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
				std::cerr << "---------------------------------------" << std::endl;

				exit (EXIT_FAILURE);
			}
        }

        solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffer);

        MPI_Isend(
                sendBuffer,
                sendBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                simulationStepCounter,
                MPI_COMM_WORLD,
                &request[0]);
        MPI_Irecv(
                recvBuffer,
                recvBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                simulationStepCounter,
                MPI_COMM_WORLD,
                &request[1]);
        MPI_Waitall(2, request, MPI_STATUS_IGNORE);

        if (configuration->doLogging)
        {
			std::stringstream loggingFileName;
			loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
			std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
			if (loggingFile.is_open())
			{
				loggingFile << "Alpha synchronization successful." << std::endl;
				loggingFile << "---------------------------------------" << std::endl;
				loggingFile.close();
			} else {
				std::cerr << "----- CController<T>::syncAlpha() -----" << std::endl;
				std::cerr << "There is no open file to write logs." << std::endl;
				std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
				std::cerr << "---------------------------------------" << std::endl;

				exit (EXIT_FAILURE);
			}
        }

        solverGPU->setDensityDistributions(recvOrigin, recvSize, recvBuffer);

        delete[] recvBuffer;
        delete[] sendBuffer;
    }
}
#endif
*/

/*
#ifdef USE_MPI
template <class T>
void CController<T>::syncBeta()
{
    CVector<3, int> sendSize, recvSize;
    CVector<3, int> sendOrigin, recvOrigin;
    int dstId;
    MPI_Request request[2];
    int sendBufferSize, recvBufferSize;
    T* sendBuffer;
    T* recvBuffer;

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "----- CController<T>::syncBeta() -----" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::syncBeta() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "--------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
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
			std::stringstream loggingFileName;
			loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
			std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
			if (loggingFile.is_open())
			{
				loggingFile << "destination rank:           " << dstId << std::endl;
				loggingFile << "send origin (with halo):    " << sendOrigin << std::endl;
				loggingFile << "receive origin (with halo): " << recvOrigin << std::endl;
				loggingFile << "send buffer size:           " << ((T)(sendBufferSize * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "receive buffer size:        " << ((T)(sendBufferSize * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
				loggingFile << "--------------------------------------" << std::endl;
				loggingFile.close();
			} else {
				std::cerr << "----- CController<T>::syncBeta() -----" << std::endl;
				std::cerr << "There is no open file to write logs." << std::endl;
				std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
				std::cerr << "--------------------------------------" << std::endl;

				exit (EXIT_FAILURE);
			}
        }

        solverGPU->getDensityDistributions(sendOrigin, sendSize, sendBuffer);

        MPI_Isend(
                sendBuffer,
                sendBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                simulationStepCounter,
                MPI_COMM_WORLD,
                &request[0]);
        MPI_Irecv(
                recvBuffer,
                recvBufferSize,
                ((typeid(T) == typeid(float)) ? MPI_FLOAT : MPI_DOUBLE),
                dstId,
                simulationStepCounter,
                MPI_COMM_WORLD,
                &request[1]);
        MPI_Waitall(2, request, MPI_STATUS_IGNORE);

        if (configuration->doLogging)
        {
			std::stringstream loggingFileName;
			loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
			std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
			if (loggingFile.is_open())
			{
				loggingFile << "Beta synchronization successful." << std::endl;
				loggingFile << "--------------------------------------" << std::endl;
				loggingFile.close();
			} else {
				std::cerr << "----- CController<T>::syncBeta() -----" << std::endl;
				std::cerr << "There is no open file to write logs." << std::endl;
				std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
				std::cerr << "--------------------------------------" << std::endl;

				exit (EXIT_FAILURE);
			}
        }

        solverGPU->setDensityDistributions(recvOrigin, recvSize, it->getDirection(), recvBuffer);

        delete[] recvBuffer;
        delete[] sendBuffer;
    }
}
#endif
*/

template <class T>
void CController<T>::computeNextStep()
{
    if (simulationStepCounter & 1) {
        stepBeta();
    } else {
        stepAlpha();
    }
    simulationStepCounter++;
}

template <class T>
void CController<T>::setDrivenCavitySzenario()
{
    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "----- CController<T>::setDrivenCavitySzenario() -----" << std::endl;
			loggingFile << "Cells where velocity injection takes place due to driven cavity scenario are marked accordingly." << std::endl;
			loggingFile << "-----------------------------------------------------" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::setDrivenCavitySzenario() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "-----------------------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
    }

    CVector<3, int> origin(1, solverCPU->getDomain()->getSizeWithHalo()[1] - 2, 1);
    CVector<3, int> size(solverCPU->getDomain()->getSize()[0], 1, solverCPU->getDomain()->getSize()[2]);

    if (configuration->doLogging)
    {
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "origin: " << origin << std::endl;
			loggingFile << "size:   " << size << std::endl;
			loggingFile << "-----------------------------------------------------" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::setDrivenCavitySzenario() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "-----------------------------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
    }

    Flag* src = new Flag[size.elements()];

    for (int i = 0; i < size.elements(); i++)
    {
        src[i] = VELOCITY_INJECTION;
    }

    // solverGPU->setFlags(origin, size, src);
    solverCPU->setFlags(origin, size, src);

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
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "----- CController<T>::run() -----" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::run() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "---------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
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

    if (configuration->doVisualization)
        visualization->render(0);

    if (configuration->doBenchmark)
        gettimeofday(&start, NULL);

    for (int i = 0; i < configuration->loops || configuration->loops < 0; i++) {
        if (configuration->doLogging)
        {
    		std::stringstream loggingFileName;
    		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
    		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
    		if (loggingFile.is_open())
    		{
    			loggingFile << "Do iteration " << i << ":" << std::endl;
    			loggingFile.close();
    		} else {
    			std::cerr << "----- CController<T>::run() -----" << std::endl;
    			std::cerr << "There is no open file to write logs." << std::endl;
    			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
    			std::cerr << "---------------------------------" << std::endl;

    			exit (EXIT_FAILURE);
    		}
        }

        computeNextStep();

        if (configuration->doVisualization)
            visualization->render(simulationStepCounter);
            
        if (configuration->doLogging)
        {
    		std::stringstream loggingFileName;
    		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
    		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
    		if (loggingFile.is_open())
    		{
    			loggingFile << "Iteration " << i << " successful." << std::endl;
    			loggingFile.close();
    		} else {
    			std::cerr << "----- CController<T>::run() -----" << std::endl;
    			std::cerr << "There is no open file to write logs." << std::endl;
    			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
    			std::cerr << "---------------------------------" << std::endl;

    			exit (EXIT_FAILURE);
    		}
        }
    }

    if (configuration->doBenchmark)
    {
        gettimeofday(&end, NULL);
        elapsed = (T)(end.tv_sec - start.tv_sec) + (T)(end.tv_usec - start.tv_usec) * (T)0.000001;

#ifdef USE_MPI
    	int numOfRanks;

        MPI_Comm_size(MPI_COMM_WORLD, &numOfRanks);
#endif

        T iterationsPerSecond = (T)(configuration->loops) / elapsed;
        T lups = iterationsPerSecond * (T)configuration->domainSize[0] * (T)configuration->domainSize[1] * (T)configuration->domainSize[2] * (T)0.000000001;
        T bandwidth = lups * (T)usedDataSize;

        std::stringstream benchmarkFileName;
        benchmarkFileName << configuration->benchmarkOutputDir << "/benchmark_";
#ifdef USE_MPI
        benchmarkFileName << numOfRanks << "_";
#endif
        benchmarkFileName << id << ".txt";
        std::ofstream benchmarkFile(benchmarkFileName.str().c_str(), std::ios::out);
        if (benchmarkFile.is_open())
        {
            benchmarkFile << "loops:           " << configuration->loops << std::endl;
            benchmarkFile << "time:            " << elapsed << "s" << std::endl;
            benchmarkFile << "iterations:      " << iterationsPerSecond << "s^-1" << std::endl;
            benchmarkFile << "lattice updates: " << lups << "GLUPS" << std::endl;
            benchmarkFile << "bandwidth:       " << bandwidth << "GBytes/s" << std::endl;
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
		std::stringstream loggingFileName;
		loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
		std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
		if (loggingFile.is_open())
		{
			loggingFile << "---------------------------------" << std::endl;
			loggingFile.close();
		} else {
			std::cerr << "----- CController<T>::run() -----" << std::endl;
			std::cerr << "There is no open file to write logs." << std::endl;
			std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
			std::cerr << "---------------------------------" << std::endl;

			exit (EXIT_FAILURE);
		}
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
    return solverCPU;
    // return solverGPU;
}

template class CController<double>;
template class CController<float>;
