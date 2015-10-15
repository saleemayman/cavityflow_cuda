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

#include <mpi.h>

#include "libtools/CStopwatch.hpp"
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
    solverGPU = new CLbmSolverGPU<T>(
            this->id,
            this->configuration->threadsPerBlock,
            this->domain,
            this->boundaryConditions,
            this->configuration->timestep,
            this->configuration->gravitation,
            this->configuration->cavityVelocity,
            this->configuration->viscosity,
            this->configuration->tau,
            this->configuration->massExchangeFactor,
            this->configuration->maxGravitationDimLess,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doLogging);
    /*
    solverCPU = new CLbmSolverCPU<T>(
            this->id,
            this->domain,
            this->boundaryConditions,
            solverGPU,
            this->configuration->timestep,
            this->configuration->gravitation,
            this->configuration->cavityVelocity,
            this->configuration->viscosity,
            this->configuration->tau,
            this->configuration->massExchangeFactor,
            this->configuration->maxGravitationDimLess,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doValidation || this->configuration->doVisualization,
            this->configuration->doLogging);
    */

	if (this->configuration->doVisualization)
        visualization = new CLbmVisualizationVTK<T>(id, getSolver(), this->configuration->visualizationOutputDir);
}

template <class T>
CController<T>::~CController()
{
    if (this->configuration->doVisualization)
        delete visualization;

    // delete solverCPU;
    delete solverGPU;
}

template <class T>
void CController<T>::syncAlpha()
{
    /*
     * TODO
     * Rework this function.
     */
#if DEBUG
    // std::cout << "--> Sync alpha" << std::endl;
#endif
    // TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking way.
    int i = 0;
    typename std::vector<CComm<T> >::iterator it = communication.begin();
    for (; it != communication.end(); it++, i++)
    {
        CVector<3, int> send_size = it->getSendSize();
        CVector<3, int> recv_size = it->getRecvSize();
        // CVector<3, int> send_origin = (*it)->getSendOrigin();
        // CVector<3, int> recv_origin = (*it)->getRecvOrigin();
        int dst_rank = it->getDstId();
#if DEBUG
        // std::cout << "ALPHA RANK: " << dst_rank << std::endl;
#endif

        // send buffer
        int send_buffer_size = send_size.elements() * NUM_LATTICE_VECTORS;
        int recv_buffer_size = recv_size.elements() * NUM_LATTICE_VECTORS;
        T* send_buffer = new T[send_buffer_size];
        T* recv_buffer = new T[recv_buffer_size];

        //printf("syncAlph--> rank: %i, send_size: %li, recv_size: %li \n", _UID, send_size.elements(), recv_size.elements());

        MPI_Request req[2];
        MPI_Status status[2];

        // Download data from device to host
        solverGPU->getDensityDistributions(Direction(i), send_buffer);
        //cLbmPtr->wait();
        int my_rank, num_procs;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /// Get current process id
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs); /// get number of processes

        if (typeid(T) == typeid(float)) {
            MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank,
                    MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank,
                    MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
        } else if (typeid(T) == typeid(double)) {
            MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank,
                    MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank,
                    MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
        } else {
            throw "Type id of MPI send/receive buffer is unknown!";
        }
        MPI_Waitall(2, req, status);

        solverGPU->setDensityDistributions(Direction(i), recv_buffer);
        //cLbmPtr->wait();

        delete[] send_buffer;
        delete[] recv_buffer;
    }
}

template <class T>
void CController<T>::syncBeta()
{
    /*
     * TODO
     * Rework this function.
     */
#if DEBUG
    // std::cout << "--> Sync beta" << std::endl;
#endif
    int i = 0;
    // TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking form.
    typename std::vector<CComm<T> >::iterator it = communication.begin();
    for (; it != communication.end(); it++, i++)
    {
        // the send and receive values in beta sync is the opposite values of
        // Comm instance related to current communication, since the ghost layer data
        // will be sent back to their origin
        CVector<3, int> send_size = it->getRecvSize();
        CVector<3, int> recv_size = it->getSendSize();
        // CVector<3, int> send_origin = (*it)->getRecvOrigin();
        // CVector<3, int> recv_origin = (*it)->getSendOrigin();
        // CVector<3, int> normal = (*it)->getCommDirection();
        int dst_rank = it->getDstId();
#if DEBUG
        // std::cout << "BETA RANK: " << dst_rank << std::endl;
#endif
        // send buffer
        int send_buffer_size = send_size.elements() * NUM_LATTICE_VECTORS;
        int recv_buffer_size = recv_size.elements() * NUM_LATTICE_VECTORS;
        T* send_buffer = new T[send_buffer_size];
        T* recv_buffer = new T[recv_buffer_size];

        //printf("syncBeta--> rank: %i, send_size: %li, recv_size: %li \n", _UID, send_size.elements(), recv_size.elements());

        MPI_Request req[2];
        MPI_Status status[2];

        // Download data from device to host
        solverGPU->getDensityDistributions(Direction(i), send_buffer);
        //cLbmPtr->wait();
        int my_rank, num_procs;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /// Get current process id
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs); /// get number of processes

        if (typeid(T) == typeid(float)) {
            MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank,
                    MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank,
                    MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
        } else if (typeid(T) == typeid(double)) {
            MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank,
                    MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank,
                    MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
        } else {
            throw "Type id of MPI send/receive buffer is unknown!";
        }
        MPI_Waitall(2, req, status);

        // TODO: OPTIMIZATION: you need to wait only for receiving to execute following command
        solverGPU->setDensityDistributions(Direction(i), recv_buffer);
        //cLbmPtr->wait();

        delete[] send_buffer;
        delete[] recv_buffer;
    }
}

template <class T>
void CController<T>::computeNextStep()
{
    if (simulationStepCounter & 1) {
        solverGPU->simulationStepBeta();
        // syncBeta();
    } else {
        solverGPU->simulationStepAlpha();
        // syncAlpha();
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

    // solverGPU->setFlags(origin, size, src);

    delete[] src;
}

/*
 * This function starts the simulation for the particular subdomain corresponded to
 * this class.
 */
template <class T>
int CController<T>::run()
{
    int usedDataSize;
    CStopwatch cStopwatch;

    usedDataSize = 0;

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
        cStopwatch.start();

    for (int i = 0; i < configuration->loops; i++) {
        computeNextStep();

        if (configuration->doVisualization)
            visualization->render(i);
    }

    if (configuration->doBenchmark)
    {
        cStopwatch.stop();

        T iterationsPerSecond = (T)configuration->loops / (T)cStopwatch.time;
        T glups = iterationsPerSecond * (T)configuration->domainSize.elements() * (T)0.000000001;
        T gbandwidth = glups * (T)usedDataSize;

        std::stringstream benchmarkFileName;
        benchmarkFileName << configuration->benchmarkOutputDir << "/benchmark_" << id << ".txt";
        std::ofstream benchmarkFile(benchmarkFileName.str().c_str(), std::ios::out);
        if (benchmarkFile.is_open())
        {
            benchmarkFile << "loops :          " << configuration->loops << std::endl;
            benchmarkFile << "time:            " << cStopwatch.time << "s" << std::endl;
            benchmarkFile << "iterations:      " << iterationsPerSecond << "s^-1" << std::endl;
            benchmarkFile << "lattice updates: " << glups << "GLUPS" << std::endl;
            benchmarkFile << "bandwidth:       " << gbandwidth << "GB/s" << std::endl;
            benchmarkFile.close();
        } else {
            /*
             * TODO
             */
        }
    }

    return EXIT_SUCCESS;
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
    return solverGPU;
}

template class CController<double>;
template class CController<float>;
