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

#include <iomanip>
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
    		MASS_EXCHANGE_FACTOR,
    		MAX_SIM_GRAVITATION_LENGTH,
    		TAU,
    		this->configuration->doValidation || this->configuration->doVisualization,
    		this->configuration->doValidation || this->configuration->doVisualization,
    		this->configuration->doLogging);
    solverCPU = new CLbmSolverCPU<T>(
    		this->id,
    		this->domain,
    		this->boundaryConditions,
    		this->configuration->timestep,
    		this->configuration->gravitation,
    		this->configuration->cavityVelocity,
    		this->configuration->viscosity,
    		solverGPU,
    		MASS_EXCHANGE_FACTOR,
    		MAX_SIM_GRAVITATION_LENGTH,
    		TAU,
    		this->configuration->doValidation || this->configuration->doVisualization,
    		this->configuration->doValidation || this->configuration->doVisualization,
    		this->configuration->doLogging);

    setDrivenCavitySzenario();

    if(this->configuration->doVisualization) {
    	cLbmVisualization = new CLbmVisualizationVTK<T>(this->id, this->configuration->visualizationOutputDir);
    }
}

template <class T>
CController<T>::~CController()
{
    if(this->configuration->doVisualization)
    	delete cLbmVisualization;

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
    	solverGPU->simulationStepAlpha();
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
int CController<T>::run()
{
	/*
	 * TODO
	 * Rework this function.
	 */
    CVector<3, int> domain_size = domain.getSize();
    int loops = configuration->loops;
    if (loops < 0)
        loops = 100;

    vector_checksum = 0;

    // approximate bandwidth
    double floats_per_cell = 0.0;

    // 19 density distribution which are read and written
    floats_per_cell += 19.0 * 2.0;

    // flag (obstacle, injection and fluid) is read
    floats_per_cell += 1.0;

    // velocity vector is also stored
    if (configuration->doVisualization
            /*|| configuration->debug_mode*/)
        floats_per_cell += 3;
    CStopwatch cStopwatch;

    // setting up the visualization
    std::string outputfilename = "OUTPUT";
    std::stringstream ss_file;
    ss_file << "./" << configuration->visualizationOutputDir << "/" << outputfilename;
    std::string outputfile = ss_file.str();
    if (configuration->doVisualization) {
        cLbmVisualization = new CLbmVisualizationVTK<T>(id, outputfile);
        cLbmVisualization->setup(solverGPU);
    }

    cStopwatch.start();
    for (int i = 0; i < loops; i++) {
        computeNextStep();
        //simulation
        if (configuration->doVisualization // && (i %  100 == 0)
            )
            cLbmVisualization->render(i);
    }
    /*
     * TODO
     * Check if this comment out affects the correctness of the code
     */
    // cLbmPtr->wait();
    cStopwatch.stop();
#if DEBUG
    /*
    if (domain_size.elements() <= 512) {
        solverGPU->debug_print();
    }
    */
#endif

#if BENCHMARK
    double ltime = cStopwatch.time;
    double gtime;
    //int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    MPI_Reduce(&ltime, &gtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (id == 0) {
        double gfps = (((double)loops) / gtime);
        double gmlups = ((double)gfps*(double)configuration->domain_size.elements())*(double)0.000001;
        double gbandwidth = (gmlups*floats_per_cell*(double)sizeof(T));
        std::stringstream benchmark_file_name;
        benchmark_file_name << "./" << BENCHMARK_OUTPUT_DIR << "/" <<
        "benchmark_" << configuration->subdomain_num.elements() //<< "_" << id
        << ".ini";
        const std::string& tmp = benchmark_file_name.str();
        const char* cstr = tmp.c_str();
        std::ofstream benchmark_file (cstr, std::ios::out | std::ios::app );
        if (benchmark_file.is_open())
        {
            //benchmark_file << "[RESULTS]" << std::endl;
            benchmark_file << "CUBE_X : " << configuration->domain_size[0] << std::endl;
            benchmark_file << "CUBE_Y : " << configuration->domain_size[1] << std::endl;
            benchmark_file << "CUBE_Z : " << configuration->domain_size[2] << std::endl;
            benchmark_file << "SECONDS : " << gtime << std::endl;
            benchmark_file << "FPS : " << gfps << std::endl;
            benchmark_file << "MLUPS : " << gmlups << std::endl;
            benchmark_file << "BANDWIDTH : " << gbandwidth// << " MB/s (RW, bidirectional)"
            << std::endl;
            benchmark_file << std::endl;

        }
        else std::cout << "Unable to open file";
    }
#endif // end of BENCHMARK
    std::cout << std::endl;
    std::cout << "Cube: " << domain_size << std::endl;
    std::cout << "Seconds: " << cStopwatch.time << std::endl;
    double fps = (((double) loops) / cStopwatch.time);
    std::cout << "FPS: " << fps << std::endl;
    double mlups =
            ((double) fps * (double) domain.getNumOfCells())
                    * (double) 0.000001;
    std::cout << "MLUPS: " << mlups << std::endl;
    std::cout << "Bandwidth: "
            << (mlups * floats_per_cell * (double) sizeof(T))
            << " MB/s (RW, bidirectional)" << std::endl;
    std::streamsize ss = std::cout.precision();
    std::cout.precision(8);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
#if 1
    // The velocity checksum is only stored in debug mode!
    /*
     * TODO
     * Alternative to getVelocityChecksum() with equivalent behavior has to be coded in CLbmSolver
     */
    // vector_checksum = cLbmPtr->getVelocityChecksum();
    std::cout << "Checksum: " << (vector_checksum*1000.0f) << std::endl;
#endif // end of DEBUG
    std::cout.precision(ss);
    std::cout << std::resetiosflags(std::ios::fixed);

    std::cout << "done." << std::endl;
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
