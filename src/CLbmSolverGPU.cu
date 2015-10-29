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

#include "CLbmSolverGPU.cuh"

#include <cassert>

#include "gpukernels/copy_buffer_rect.cuh"
#include "gpukernels/lbm_alpha.cuh"
#include "gpukernels/lbm_beta.cuh"
#include "gpukernels/lbm_init.cuh"

template <class T>
CLbmSolverGPU<T>::CLbmSolverGPU(
		int id,
		std::vector<dim3> threadsPerBlock,
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
        CLbmSolver<T>(id, domain,
                boundaryConditions,
                timestepSize, gravitation, drivenCavityVelocity,
                viscosity, maxGravitationDimLess,
                storeDensities, storeVelocities, doLogging),
        threadsPerBlock(threadsPerBlock)
{
    int numOfGPUsPerNode;

    GPU_ERROR_CHECK(cudaGetDeviceCount(&numOfGPUsPerNode))
    GPU_ERROR_CHECK(cudaSetDevice(this->id % numOfGPUsPerNode))

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::CLbmSolverGPU() -----" << std::endl;
        std::cout << "id:                                                 " << this->id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "number of GPUs per node:                            " << numOfGPUsPerNode << std::endl;
        std::cout << "number of selected GPU:                             " << (this->id % numOfGPUsPerNode) << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&densityDistributions, this->domain.getNumOfCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&flags, this->domain.getNumOfCells() * sizeof(Flag)))
    if(this->storeDensities)
        GPU_ERROR_CHECK(cudaMalloc(&densities, this->domain.getNumOfCells() * sizeof(T)))
    if(this->storeVelocities)
        GPU_ERROR_CHECK(cudaMalloc(&velocities, this->domain.getNumOfCells() * 3 * sizeof(T)))

    if (doLogging) {
        std::cout << "size of allocated memory for density distributions: " << ((T)(this->domain.getNumOfCells() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        std::cout << "size of allocated memory for flags:                 " << ((T)(this->domain.getNumOfCells() * sizeof(Flag)) / (T)(1<<20)) << " MBytes" << std::endl;
        if(this->storeDensities)
            std::cout << "size of allocated memory for velocities:            " << ((T)(this->domain.getNumOfCells() * 3 * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        if(this->storeVelocities)
            std::cout << "size of allocated memory for densities:             " << ((T)(this->domain.getNumOfCells() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
    }

	getDensityDistributionsHalo.resize(6);
	setDensityDistributionsHalo.resize(6);

    GPU_ERROR_CHECK(cudaMalloc(&(getDensityDistributionsHalo[0]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(setDensityDistributionsHalo[0]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(getDensityDistributionsHalo[1]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(setDensityDistributionsHalo[1]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(getDensityDistributionsHalo[2]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(setDensityDistributionsHalo[2]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(getDensityDistributionsHalo[3]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(setDensityDistributionsHalo[3]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(getDensityDistributionsHalo[4]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(setDensityDistributionsHalo[4]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(getDensityDistributionsHalo[5]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&(setDensityDistributionsHalo[5]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))

    if (doLogging) {
        std::cout << "size of allocated memory for LEFT/RIGHT halos:      " << ((T)(4 * this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        std::cout << "size of allocated memory for BOTTOM/TOP halos:      " << ((T)(4 * this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        std::cout << "size of allocated memory for BACK/FRONT halos:      " << ((T)(4 * this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    dim3 blocksPerGrid = getBlocksPerGrid(3, this->domain.getSize(), this->threadsPerBlock[0]);

    if (doLogging) {
        std::cout << "threads per block:                                  [" << this->threadsPerBlock[0].x << ", " << this->threadsPerBlock[0].y << ", " << this->threadsPerBlock[0].z << "]" << std::endl;
        std::cout << "blocks per grid:                                    [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    lbm_init<T><<<blocksPerGrid, this->threadsPerBlock[0]>>>(
        densityDistributions,
        flags,
        velocities,
        densities,
        boundaryConditions[0],
        boundaryConditions[1],
        boundaryConditions[2],
        boundaryConditions[3],
        boundaryConditions[4],
        boundaryConditions[5],
        drivenCavityVelocityDimLess.data[0],
        domain.getSize()[0],
        domain.getSize()[1],
        domain.getSize()[2]);
    GPU_ERROR_CHECK(cudaPeekAtLastError())
    
    if (doLogging) {
        std::cout << "Domain successfully initialized." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
}

template <class T>
CLbmSolverGPU<T>::~CLbmSolverGPU()
{
    GPU_ERROR_CHECK(cudaFree(setDensityDistributionsHalo[5]))
    GPU_ERROR_CHECK(cudaFree(getDensityDistributionsHalo[5]))
    GPU_ERROR_CHECK(cudaFree(setDensityDistributionsHalo[4]))
    GPU_ERROR_CHECK(cudaFree(getDensityDistributionsHalo[4]))
    GPU_ERROR_CHECK(cudaFree(setDensityDistributionsHalo[3]))
    GPU_ERROR_CHECK(cudaFree(getDensityDistributionsHalo[3]))
    GPU_ERROR_CHECK(cudaFree(setDensityDistributionsHalo[2]))
    GPU_ERROR_CHECK(cudaFree(getDensityDistributionsHalo[2]))
    GPU_ERROR_CHECK(cudaFree(setDensityDistributionsHalo[1]))
    GPU_ERROR_CHECK(cudaFree(getDensityDistributionsHalo[1]))
    GPU_ERROR_CHECK(cudaFree(setDensityDistributionsHalo[0]))
    GPU_ERROR_CHECK(cudaFree(getDensityDistributionsHalo[0]))

    if(storeVelocities)
        GPU_ERROR_CHECK(cudaFree(velocities))
    if(storeDensities)
        GPU_ERROR_CHECK(cudaFree(densities))
    GPU_ERROR_CHECK(cudaFree(flags))
    GPU_ERROR_CHECK(cudaFree(densityDistributions))
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlpha()
{
    dim3 blocksPerGrid = getBlocksPerGrid(3, this->domain.getSize(), this->threadsPerBlock[1]);

    if (doLogging)
    {
		std::cout << "----- CLbmSolverGPU<T>::simulationStepAlpha() -----" << std::endl;
		std::cout << "id:                " << id << std::endl;
		std::cout << "---------------------------------------------------" << std::endl;
		std::cout << "threads per block: [" << this->threadsPerBlock[1].x << ", " << this->threadsPerBlock[1].y << ", " << this->threadsPerBlock[1].z << "]" << std::endl;
		std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
		std::cout << "---------------------------------------------------" << std::endl;
    }

    lbm_kernel_alpha<T><<<blocksPerGrid, this->threadsPerBlock[1]>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            gravitationDimLess[0],
            gravitationDimLess[1],
            gravitationDimLess[2],
            drivenCavityVelocityDimLess[0],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            storeDensities,
            storeVelocities);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    if (doLogging)
    {
		std::cout << "Alpha kernel was successfully executed on the whole subdomain." << std::endl;
		std::cout << "---------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlphaRect(CVector<3, int> origin, CVector<3, int> size)
{
    /*
     * TODO
     */
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBeta()
{
    dim3 blocksPerGrid = getBlocksPerGrid(3, domain.getSize(), threadsPerBlock[2]);
    size_t sMemSize = 12 * sizeof(T) * getSize(threadsPerBlock[2]);

    if (doLogging)
    {
		std::cout << "----- CLbmSolverGPU<T>::simulationStepBeta() -----" << std::endl;
		std::cout << "id:                 " << id << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "threads per block:  [" << this->threadsPerBlock[2].x << ", " << this->threadsPerBlock[2].y << ", " << this->threadsPerBlock[2].z << "]" << std::endl;
		std::cout << "blocks per grid:    [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
		std::cout << "shared memory size: " << ((T)sMemSize / (T)(1<<10)) << " KB" << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
    }

    lbm_kernel_beta<T><<<blocksPerGrid, threadsPerBlock[2], sMemSize>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            gravitationDimLess[0],
            gravitationDimLess[1],
            gravitationDimLess[2],
            drivenCavityVelocityDimLess[0],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            getSize(threadsPerBlock[2]),
            isPowerOfTwo(domain.getNumOfCells()),
            isPowerOfTwo(getSize(threadsPerBlock[2])),
            storeDensities,
            storeVelocities);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    if (doLogging)
    {
		std::cout << "Beta kernel was successfully executed on the whole subdomain." << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBetaRect(CVector<3, int> origin, CVector<3, int> size)
{
    /*
     * TODO
     */
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    T* dDensityDistributions;
    cudaMemcpy3DParms params = {0};

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getDensityDistributions() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T)))

    for(int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
    	// domain location and size
    	params.srcPtr = make_cudaPitchedPtr(densityDistributions, domain.getSize()[0] * sizeof(T), domain.getSize()[0], domain.getSize()[1]);
    	// cuboid origin
    	params.srcPos = make_cudaPos(origin[0], origin[1], origin[2]);
    	// hDensityDistributions location and size
    	params.dstPtr = make_cudaPitchedPtr(hDensityDistributions, size[0] * sizeof(T), size[0], size[1]);
    	// hDensityDistributions origin
    	params.dstPos = make_cudaPos(0, 0, 0);
    	// cuboid size
    	params.extent = make_cudaExtent(size[0], size[1], size[2]);
    	params.kind = cudaMemcpyDeviceToHost;

    	// cudaMemcpy3D(&params);

        copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
                densityDistributions,
                latticeVector * domain.getNumOfCells(),
                origin[0],
                origin[1],
                origin[2],
                domain.getSize()[0],
                domain.getSize()[1],
                domain.getSize()[2],
                dDensityDistributions,
                latticeVector * size.elements(),
                0,
                0,
                0,
                size[0],
                size[1],
                size[2],
                size[0],
                size[1],
                size[2]);
        GPU_ERROR_CHECK(cudaPeekAtLastError())
    }

    GPU_ERROR_CHECK(cudaMemcpy(hDensityDistributions, dDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyDeviceToHost))

    GPU_ERROR_CHECK(cudaFree(dDensityDistributions))

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(T* hDensityDistributions)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    getDensityDistributions(origin, size, hDensityDistributions);
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, Direction direction, T* hDensityDistributions)
{
    assert(0 <= direction < 6);

    CVector<3, int> norm(0);

    switch(direction)
    {
    case LEFT:
        norm[0] = 1;
        break;
    case RIGHT:
        norm[0] = -1;
        break;
    case BOTTOM:
        norm[1] = 1;
        break;
    case TOP:
        norm[1] = -1;
        break;
    case BACK:
        norm[2] = 1;
        break;
    case FRONT:
        norm[2] = -1;
        break;
    }

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging) {
        std::cout << "----- CLbmSolverGPU<T>::setDensityDistributions() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "direction:         " << norm << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMemcpy(setDensityDistributionsHalo[direction], hDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        if(norm.dotProd(lbm_units[latticeVector]) > 0)
        {
            copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
                    setDensityDistributionsHalo[direction],
                    latticeVector * size.elements(),
                    0,
                    0,
                    0,
                    size[0],
                    size[1],
                    size[2],
                    densityDistributions,
                    latticeVector * domain.getNumOfCells(),
                    origin[0],
                    origin[1],
                    origin[2],
                    domain.getSize()[0],
                    domain.getSize()[1],
                    domain.getSize()[2],
                    size[0],
                    size[1],
                    size[2]);
        }
        GPU_ERROR_CHECK(cudaPeekAtLastError())
    }

    if (doLogging) {
        std::cout << "A copy operation from host to device for lattice vectors in direction " << direction << " was performed." << std::endl;
        std::cout << "No additional buffer memory was allocated. Instead, setDensityDistributionsHalo was used." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    T* dDensityDistributions;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T)))

    GPU_ERROR_CHECK(cudaMemcpy(dDensityDistributions, hDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
                dDensityDistributions,
                latticeVector * size.elements(),
                0,
                0,
                0,
                size[0],
                size[1],
                size[2],
                densityDistributions,
                latticeVector * domain.getNumOfCells(),
                origin[0],
                origin[1],
                origin[2],
                domain.getSize()[0],
                domain.getSize()[1],
                domain.getSize()[2],
                size[0],
                size[1],
                size[2]);
        GPU_ERROR_CHECK(cudaPeekAtLastError())
    }

    GPU_ERROR_CHECK(cudaFree(dDensityDistributions))

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(T* hDensityDistributions)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    setDensityDistributions(origin, size, hDensityDistributions);
}

template <class T>
void CLbmSolverGPU<T>::getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* hFlags)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    Flag* dFlags;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getFlags() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dFlags, size.elements() * sizeof(Flag)))

    copy_buffer_rect<Flag><<<blocksPerGrid, threadsPerBlock>>>(
            flags,
            0,
            origin[0],
            origin[1],
            origin[2],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            dFlags,
            0,
            0,
            0,
            0,
            size[0],
            size[1],
            size[2],
            size[0],
            size[1],
            size[2]);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    GPU_ERROR_CHECK(cudaMemcpy(hFlags, dFlags, size.elements() * sizeof(Flag), cudaMemcpyDeviceToHost))

    GPU_ERROR_CHECK(cudaFree(dFlags))

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getFlags(Flag* hFlags)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    getFlags(origin, size, hFlags);
}

template <class T>
void CLbmSolverGPU<T>::setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* hFlags)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    Flag* dFlags;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setFlags() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dFlags, size.elements() * sizeof(Flag)))

    GPU_ERROR_CHECK(cudaMemcpy(dFlags, hFlags, size.elements() * sizeof(Flag), cudaMemcpyHostToDevice))

    copy_buffer_rect<Flag><<<blocksPerGrid, threadsPerBlock>>>(
            dFlags,
            0,
            0,
            0,
            0,
            size[0],
            size[1],
            size[2],
            flags,
            0,
            origin[0],
            origin[1],
            origin[2],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            size[0],
            size[1],
            size[2]);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    GPU_ERROR_CHECK(cudaFree(dFlags))

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setFlags(Flag* hFlags)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    setFlags(origin, size, hFlags);
}

template <class T>
void CLbmSolverGPU<T>::getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* hVelocities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    T* dVelocities;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getVelocities() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dVelocities, 3 * size.elements() * sizeof(T)))

    for (int dim = 0; dim < 3; dim++)
    {
        copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
                velocities,
                dim * domain.getNumOfCells(),
                origin[0],
                origin[1],
                origin[2],
                domain.getSize()[0],
                domain.getSize()[1],
                domain.getSize()[2],
                dVelocities,
                dim * size.elements(),
                0,
                0,
                0,
                size[0],
                size[1],
                size[2],
                size[0],
                size[1],
                size[2]);
        GPU_ERROR_CHECK(cudaPeekAtLastError())
    }

    GPU_ERROR_CHECK(cudaMemcpy(hVelocities, dVelocities, 3 * size.elements() * sizeof(T), cudaMemcpyDeviceToHost))

    GPU_ERROR_CHECK(cudaFree(dVelocities))

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getVelocities(T* hVelocities)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    getVelocities(origin, size, hVelocities);
}

template <class T>
void CLbmSolverGPU<T>::setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* hVelocities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    T* dVelocities;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setVelocities() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dVelocities, 3 * size.elements() * sizeof(T)))

    GPU_ERROR_CHECK(cudaMemcpy(dVelocities, hVelocities, 3 * size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    for (int dim = 0; dim < 3; dim++)
    {
        copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
                dVelocities,
                dim * size.elements(),
                0,
                0,
                0,
                size[0],
                size[1],
                size[2],
                velocities,
                dim * domain.getNumOfCells(),
                origin[0],
                origin[1],
                origin[2],
                domain.getSize()[0],
                domain.getSize()[1],
                domain.getSize()[2],
                size[0],
                size[1],
                size[2]);
        GPU_ERROR_CHECK(cudaPeekAtLastError())
    }

    GPU_ERROR_CHECK(cudaFree(dVelocities))

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setVelocities(T* hVelocities)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    setVelocities(origin, size, hVelocities);
}

template <class T>
void CLbmSolverGPU<T>::getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* hDensities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    T* dDensities;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getDensities() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dDensities, size.elements() * sizeof(T)))

    copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
            densities,
            0,
            origin[0],
            origin[1],
            origin[2],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            dDensities,
            0,
            0,
            0,
            0,
            size[0],
            size[1],
            size[2],
            size[0],
            size[1],
            size[2]);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    GPU_ERROR_CHECK(cudaMemcpy(hDensities, dDensities, size.elements() * sizeof(T), cudaMemcpyDeviceToHost))

    GPU_ERROR_CHECK(cudaFree(dDensities))

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getDensities(T* hDensities)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    getDensities(origin, size, hDensities);
}

template <class T>
void CLbmSolverGPU<T>::setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* hDensities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0]);
    assert(origin[1] + size[1] <= domain.getSize()[1]);
    assert(origin[2] + size[2] <= domain.getSize()[2]);

    T* dDensities;

    dim3 threadsPerBlock(size[1], size[2], 1);
    if(size[1] > 1024)
    {
    	threadsPerBlock.x = 1024;
    	threadsPerBlock.y = 1;
    }
    if(size[1] * size[2] > 1024)
    	threadsPerBlock.y = 1024 / size[1];
    dim3 blocksPerGrid = getBlocksPerGrid(2, CVector<3, int>(size[1], size[2], 0), threadsPerBlock);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setDensities() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock.x << ", " << threadsPerBlock.y << ", " << threadsPerBlock.z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&dDensities, size.elements() * sizeof(T)))

    GPU_ERROR_CHECK(cudaMemcpy(dDensities, hDensities, size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    copy_buffer_rect<T><<<blocksPerGrid, threadsPerBlock>>>(
            dDensities,
            0,
            0,
            0,
            0,
            size[0],
            size[1],
            size[2],
            densities,
            0,
            origin[0],
            origin[1],
            origin[2],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            size[0],
            size[1],
            size[2]);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    GPU_ERROR_CHECK(cudaFree(dDensities))

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensities(T* hDensities)
{
    CVector<3, int> origin(0);
    CVector<3, int> size(domain.getSize());

    setDensities(origin, size, hDensities);
}

template class CLbmSolverGPU<float>;
template class CLbmSolverGPU<double>;
