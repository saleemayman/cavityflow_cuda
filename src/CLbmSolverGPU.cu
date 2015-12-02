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

#include "gpukernels/lbm_alpha.cuh"
#include "gpukernels/lbm_beta.cuh"
#include "gpukernels/lbm_init.cuh"

template <class T>
CLbmSolverGPU<T>::CLbmSolverGPU(
        int id,
        std::vector<dim3> threadsPerBlock,
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
        CLbmSolver<T>(id, globalLength,
                domain, boundaryConditions,
                timestepSize, velocity, acceleration,
                viscosity, maxVelocityDimLess, maxAccelerationDimLess,
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

    GPU_ERROR_CHECK(cudaMalloc(&densityDistributions, NUM_LATTICE_VECTORS * this->domain.getNumOfCellsWithHalo() * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&flags, this->domain.getNumOfCellsWithHalo() * sizeof(Flag)))
    if(this->storeDensities)
        GPU_ERROR_CHECK(cudaMalloc(&densities, this->domain.getNumOfCellsWithHalo() * sizeof(T)))
    if(this->storeVelocities)
        GPU_ERROR_CHECK(cudaMalloc(&velocities, 3 * this->domain.getNumOfCellsWithHalo() * sizeof(T)))

    if (doLogging) {
        std::cout << "size of allocated memory for density distributions: " << ((T)(NUM_LATTICE_VECTORS * this->domain.getNumOfCellsWithHalo() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        std::cout << "size of allocated memory for flags:                 " << ((T)(this->domain.getNumOfCellsWithHalo() * sizeof(Flag)) / (T)(1<<20)) << " MBytes" << std::endl;
        if(this->storeDensities)
            std::cout << "size of allocated memory for velocities:            " << ((T)(3 * this->domain.getNumOfCellsWithHalo() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        if(this->storeVelocities)
            std::cout << "size of allocated memory for densities:             " << ((T)(this->domain.getNumOfCellsWithHalo() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
    }

    dim3 blocksPerGrid = getBlocksPerGrid(3, this->domain.getSizeWithHalo(), this->threadsPerBlock[0]);

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
        velocityDimLess[0],
        domain.getSizeWithHalo()[0],
        domain.getSizeWithHalo()[1],
        domain.getSizeWithHalo()[2],
        storeDensities,
        storeVelocities);
    GPU_ERROR_CHECK(cudaPeekAtLastError())
    
    if (doLogging) {
        std::cout << "Domain successfully initialized." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
}

template <class T>
CLbmSolverGPU<T>::~CLbmSolverGPU()
{
    if(storeVelocities)
        GPU_ERROR_CHECK(cudaFree(velocities))
    if(storeDensities)
        GPU_ERROR_CHECK(cudaFree(densities))
    GPU_ERROR_CHECK(cudaFree(flags))
    GPU_ERROR_CHECK(cudaFree(densityDistributions))
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlpha(cudaStream_t* stream)
{
    dim3 blocksPerGrid = getBlocksPerGrid(3, domain.getSizeWithHalo(), threadsPerBlock[1]);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::simulationStepAlpha() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock[1].x << ", " << threadsPerBlock[1].y << ", " << threadsPerBlock[1].z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

    lbm_kernel_alpha<T><<<blocksPerGrid, threadsPerBlock[1], 0, ((stream == NULL) ? 0 : *stream)>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            accelerationDimLess[0],
            accelerationDimLess[1],
            accelerationDimLess[2],
            velocityDimLess[0],
            0,
            0,
            0,
            domain.getSizeWithHalo()[0],
            domain.getSizeWithHalo()[1],
            domain.getSizeWithHalo()[2],
            domain.getSizeWithHalo()[0],
            domain.getSizeWithHalo()[1],
            domain.getSizeWithHalo()[2],
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
void CLbmSolverGPU<T>::simulationStepAlpha()
{
    simulationStepAlpha(NULL);
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlpha(CVector<3, int> origin, CVector<3, int> size, cudaStream_t* stream)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    dim3 blocksPerGrid = getBlocksPerGrid(3, size, threadsPerBlock[1]);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::simulationStepAlpha() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
        std::cout << "threads per block: [" << threadsPerBlock[1].x << ", " << threadsPerBlock[1].y << ", " << threadsPerBlock[1].z << "]" << std::endl;
        std::cout << "blocks per grid:   [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

    lbm_kernel_alpha<T><<<blocksPerGrid, threadsPerBlock[1], 0, ((stream == NULL) ? 0 : *stream)>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            accelerationDimLess[0],
            accelerationDimLess[1],
            accelerationDimLess[2],
            velocityDimLess[0],
            origin[0],
            origin[1],
            origin[2],
            size[0],
            size[1],
            size[2],
            domain.getSizeWithHalo()[0],
            domain.getSizeWithHalo()[1],
            domain.getSizeWithHalo()[2],
            storeDensities,
            storeVelocities);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    if (doLogging)
    {
        std::cout << "Alpha kernel was successfully executed on the following subdomain:" << std::endl;
        std::cout << "origin:            " << origin << std::endl;
        std::cout << "size:              " << size << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlpha(CVector<3, int> origin, CVector<3, int> size)
{
    simulationStepAlpha(origin, size, NULL);
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBeta(cudaStream_t* stream)
{
    dim3 blocksPerGrid = getBlocksPerGrid(3, domain.getSizeWithHalo(), threadsPerBlock[2]);
    size_t sMemSize = 12 * sizeof(T) * getSize(threadsPerBlock[2]);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::simulationStepBeta() -----" << std::endl;
        std::cout << "id:                 " << id << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "threads per block:  [" << threadsPerBlock[2].x << ", " << threadsPerBlock[2].y << ", " << threadsPerBlock[2].z << "]" << std::endl;
        std::cout << "blocks per grid:    [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "shared memory size: " << ((T)sMemSize / (T)(1<<10)) << " KB" << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }

    lbm_kernel_beta<T><<<blocksPerGrid, threadsPerBlock[2], sMemSize, ((stream == NULL) ? 0 : *stream)>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            accelerationDimLess[0],
            accelerationDimLess[1],
            accelerationDimLess[2],
            velocityDimLess[0],
            0,
            0,
            0,
            domain.getSizeWithHalo()[0],
            domain.getSizeWithHalo()[1],
            domain.getSizeWithHalo()[2],
            domain.getSizeWithHalo()[0],
            domain.getSizeWithHalo()[1],
            domain.getSizeWithHalo()[2],
            getSize(threadsPerBlock[2]),
            isPowerOfTwo(domain.getNumOfCellsWithHalo()),
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
void CLbmSolverGPU<T>::simulationStepBeta()
{
    simulationStepBeta(NULL);
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBeta(CVector<3, int> origin, CVector<3, int> size, cudaStream_t* stream)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    dim3 blocksPerGrid = getBlocksPerGrid(3, size, threadsPerBlock[2]);
    size_t sMemSize = 12 * sizeof(T) * getSize(threadsPerBlock[2]);

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::simulationStepBeta() -----" << std::endl;
        std::cout << "id:                 " << id << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "threads per block:  [" << threadsPerBlock[2].x << ", " << threadsPerBlock[2].y << ", " << threadsPerBlock[2].z << "]" << std::endl;
        std::cout << "blocks per grid:    [" << blocksPerGrid.x << ", " << blocksPerGrid.y << ", " << blocksPerGrid.z << "]" << std::endl;
        std::cout << "shared memory size: " << ((T)sMemSize / (T)(1<<10)) << " KB" << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }

    lbm_kernel_beta<T><<<blocksPerGrid, threadsPerBlock[2], sMemSize, ((stream == NULL) ? 0 : *stream)>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            accelerationDimLess[0],
            accelerationDimLess[1],
            accelerationDimLess[2],
            velocityDimLess[0],
            origin[0],
            origin[1],
            origin[2],
            size[0],
            size[1],
            size[2],
            domain.getSizeWithHalo()[0],
            domain.getSizeWithHalo()[1],
            domain.getSizeWithHalo()[2],
            getSize(threadsPerBlock[2]),
            isPowerOfTwo(domain.getNumOfCellsWithHalo()),
            isPowerOfTwo(getSize(threadsPerBlock[2])),
            storeDensities,
            storeVelocities);
    GPU_ERROR_CHECK(cudaPeekAtLastError())

    if (doLogging)
    {
        std::cout << "Beta kernel was successfully executed on the following subdomain." << std::endl;
        std::cout << "origin:             " << origin << std::endl;
        std::cout << "size:               " << size << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBeta(CVector<3, int> origin, CVector<3, int> size)
{
    simulationStepBeta(origin, size, NULL);
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions, cudaStream_t* stream)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

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
    }

    for(int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        // domain location and size
        params.srcPtr = make_cudaPitchedPtr(&densityDistributions[latticeVector * domain.getNumOfCellsWithHalo()], domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
        // cuboid origin
        params.srcPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
        // hDensityDistributions location and size
        params.dstPtr = make_cudaPitchedPtr(&hDensityDistributions[latticeVector * size.elements()], size[0] * sizeof(T), size[0], size[1]);
        // hDensityDistributions origin
        params.dstPos = make_cudaPos(0, 0, 0);
        // cuboid size
        params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
        params.kind = cudaMemcpyDeviceToHost;

        GPU_ERROR_CHECK(cudaMemcpy3DAsync(&params, (stream == NULL) ? 0 : *stream))
    }

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, T* hDensityDistributions)
{
    getDensityDistributions(origin, size, hDensityDistributions, NULL);
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(T* hDensityDistributions, cudaStream_t* stream)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getDensityDistributions(origin, size, hDensityDistributions, stream);
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(T* hDensityDistributions)
{
    getDensityDistributions(hDensityDistributions, NULL);
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, Direction direction, T* hDensityDistributions, cudaStream_t* stream)
{
    assert(0 <= direction < 6);
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    CVector<3, int> norm(0);
    cudaMemcpy3DParms params = {0};

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
    }

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        if(norm.dotProd(lbm_units[latticeVector]) > 0)
        {
            // hDensityDistributions location and size
            params.srcPtr = make_cudaPitchedPtr(&hDensityDistributions[latticeVector * size.elements()], size[0] * sizeof(T), size[0], size[1]);
            // hDensityDistributions origin
            params.srcPos = make_cudaPos(0, 0, 0);
            // domain location and size
            params.dstPtr = make_cudaPitchedPtr(&densityDistributions[latticeVector * domain.getNumOfCellsWithHalo()], domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
            // cuboid origin
            params.dstPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
            // cuboid size
            params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
            params.kind = cudaMemcpyHostToDevice;

            GPU_ERROR_CHECK(cudaMemcpy3DAsync(&params, (stream == NULL) ? 0 : *stream))
        }
    }

    if (doLogging) {
        std::cout << "A copy operation from host to device for lattice vectors in direction " << direction << " was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, Direction direction, T* hDensityDistributions)
{
    setDensityDistributions(origin, size, direction, hDensityDistributions, NULL);
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions, cudaStream_t* stream)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        // hDensityDistributions location and size
        params.srcPtr = make_cudaPitchedPtr(&hDensityDistributions[latticeVector * size.elements()], size[0] * sizeof(T), size[0], size[1]);
        // hDensityDistributions origin
        params.srcPos = make_cudaPos(0, 0, 0);
        // domain location and size
        params.dstPtr = make_cudaPitchedPtr(&densityDistributions[latticeVector * domain.getNumOfCellsWithHalo()], domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
        // cuboid origin
        params.dstPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
        // cuboid size
        params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
        params.kind = cudaMemcpyHostToDevice;

        GPU_ERROR_CHECK(cudaMemcpy3DAsync(&params, (stream == NULL) ? 0 : *stream))
    }

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int>& origin, CVector<3, int>& size, T* hDensityDistributions)
{
    setDensityDistributions(origin, size, hDensityDistributions, NULL);
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(T* hDensityDistributions, cudaStream_t* stream)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setDensityDistributions(origin, size, hDensityDistributions, stream);
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(T* hDensityDistributions)
{
    setDensityDistributions(hDensityDistributions, NULL);
}

template <class T>
void CLbmSolverGPU<T>::getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* hFlags)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

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
    }

    // domain location and size
    params.srcPtr = make_cudaPitchedPtr(flags, domain.getSizeWithHalo()[0] * sizeof(Flag), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.srcPos = make_cudaPos(origin[0] * (sizeof(Flag) / sizeof(unsigned char)), origin[1], origin[2]);
    // hFlags location and size
    params.dstPtr = make_cudaPitchedPtr(hFlags, size[0] * sizeof(Flag), size[0], size[1]);
    // hFlags origin
    params.dstPos = make_cudaPos(0, 0, 0);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(Flag) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyDeviceToHost;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getFlags(Flag* hFlags)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getFlags(origin, size, hFlags);
}

template <class T>
void CLbmSolverGPU<T>::setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* hFlags)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setFlags() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

    // hFlags location and size
    params.srcPtr = make_cudaPitchedPtr(hFlags, size[0] * sizeof(Flag), size[0], size[1]);
    // hFlags origin
    params.srcPos = make_cudaPos(0, 0, 0);
    // domain location and size
    params.dstPtr = make_cudaPitchedPtr(flags, domain.getSizeWithHalo()[0] * sizeof(Flag), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.dstPos = make_cudaPos(origin[0] * (sizeof(Flag) / sizeof(unsigned char)), origin[1], origin[2]);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(Flag) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyHostToDevice;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setFlags(Flag* hFlags)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setFlags(origin, size, hFlags);
}

template <class T>
void CLbmSolverGPU<T>::getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* hVelocities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getVelocities() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    for (int dim = 0; dim < 3; dim++)
    {
        // domain location and size
        params.srcPtr = make_cudaPitchedPtr(&velocities[dim * domain.getNumOfCellsWithHalo()], domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
        // cuboid origin
        params.srcPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
        // hVelocities location and size
        params.dstPtr = make_cudaPitchedPtr(&hVelocities[dim * size.elements()], size[0] * sizeof(T), size[0], size[1]);
        // hVelocities origin
        params.dstPos = make_cudaPos(0, 0, 0);
        // cuboid size
        params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
        params.kind = cudaMemcpyDeviceToHost;

        GPU_ERROR_CHECK(cudaMemcpy3D(&params))
    }

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getVelocities(T* hVelocities)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getVelocities(origin, size, hVelocities);
}

template <class T>
void CLbmSolverGPU<T>::setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* hVelocities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setVelocities() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    for (int dim = 0; dim < 3; dim++)
    {
        // hVelocities location and size
        params.srcPtr = make_cudaPitchedPtr(&hVelocities[dim * size.elements()], size[0] * sizeof(T), size[0], size[1]);
        // hVelocities origin
        params.srcPos = make_cudaPos(0, 0, 0);
        // domain location and size
        params.dstPtr = make_cudaPitchedPtr(&velocities[dim * domain.getNumOfCellsWithHalo()], domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
        // cuboid origin
        params.dstPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
        // cuboid size
        params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
        params.kind = cudaMemcpyHostToDevice;

        GPU_ERROR_CHECK(cudaMemcpy3D(&params))
    }

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setVelocities(T* hVelocities)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setVelocities(origin, size, hVelocities);
}

template <class T>
void CLbmSolverGPU<T>::getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* hDensities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getDensities() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }

    // domain location and size
    params.srcPtr = make_cudaPitchedPtr(densities, domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.srcPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
    // hDensities location and size
    params.dstPtr = make_cudaPitchedPtr(hDensities, size[0] * sizeof(T), size[0], size[1]);
    // hDensities origin
    params.dstPos = make_cudaPos(0, 0, 0);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyDeviceToHost;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))

    if (doLogging)
    {
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getDensities(T* hDensities)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getDensities(origin, size, hDensities);
}

template <class T>
void CLbmSolverGPU<T>::setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* hDensities)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::setDensities() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }

    // hDensities location and size
    params.srcPtr = make_cudaPitchedPtr(hDensities, size[0] * sizeof(T), size[0], size[1]);
    // hDensities origin
    params.srcPos = make_cudaPos(0, 0, 0);
    // domain location and size
    params.dstPtr = make_cudaPitchedPtr(densities, domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.dstPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyHostToDevice;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))

    if (doLogging)
    {
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensities(T* hDensities)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setDensities(origin, size, hDensities);
}

template class CLbmSolverGPU<float>;
template class CLbmSolverGPU<double>;
