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
        CVector<4, T> &drivenCavityVelocity,
        T viscocity,
        T massExchangeFactor,
        T maxSimGravitationLength,
        T tau,
        bool storeDensities,
        bool storeVelocities,
        bool doLogging) :
        CLbmSolver<T>(id, domain,
                boundaryConditions, timestepSize,
                gravitation, drivenCavityVelocity, viscocity,
                massExchangeFactor, maxSimGravitationLength, tau,
                storeDensities, storeVelocities,
                doLogging),
        threadsPerBlock(threadsPerBlock)
{
    int numOfGPUsPerNode;

    GPU_ERROR_CHECK(cudaGetDeviceCount(&numOfGPUsPerNode))
    GPU_ERROR_CHECK(cudaSetDevice(this->id % numOfGPUsPerNode))

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::CLbmSolverGPU() -----" << std::endl;
        std::cout << "id:                      " << this->id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "number of GPUs per node: " << numOfGPUsPerNode << std::endl;
        std::cout << "number of selected GPU:  " << (this->id % numOfGPUsPerNode) << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc(&densityDistributions, this->domain.getNumOfCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc(&flags, this->domain.getNumOfCells() * sizeof(Flag)))
    if(this->storeDensities)
        GPU_ERROR_CHECK(cudaMalloc(&densities, this->domain.getNumOfCells() * sizeof(T)))
    if(this->storeVelocities)
        GPU_ERROR_CHECK(cudaMalloc(&velocities, this->domain.getNumOfCells() * 3 * sizeof(T)))

    if (doLogging) {
        std::cout << "size of allocated memory for density distributions: " << (this->domain.getNumOfCells() * NUM_LATTICE_VECTORS * sizeof(T)) << "Bytes" << std::endl;
        std::cout << "size of allocated memory for flags:                 " << (this->domain.getNumOfCells() * sizeof(Flag)) << "Bytes" << std::endl;
        if(this->storeDensities)
            std::cout << "size of allocated memory for velocities:            " << (this->domain.getNumOfCells() * 3 * sizeof(T)) << "Bytes" << std::endl;
        if(this->storeVelocities)
            std::cout << "size of allocated memory for densities:             " << (this->domain.getNumOfCells() * sizeof(T)) << "Bytes" << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[0]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[0]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[1]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[1]), this->domain.getNumOfXFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[2]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[2]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[3]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[3]), this->domain.getNumOfYFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[4]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[4]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[5]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))
    GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[5]), this->domain.getNumOfZFaceCells() * NUM_LATTICE_VECTORS * sizeof(T)))

    lbm_init<T><<<getBlocksPerGrid(3, this->domain.getSize(), this->threadsPerBlock[0]), this->threadsPerBlock[0]>>>(
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
        drivenCavityVelocity.data[0],
        domain.getSize()[0],
        domain.getSize()[1],
        domain.getSize()[2]);
    GPU_ERROR_CHECK(cudaPeekAtLastError())
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
    lbm_kernel_alpha<T><<<getBlocksPerGrid(3, domain.getSize(), threadsPerBlock[1]), threadsPerBlock[1]>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            gravitation[0],
            gravitation[1],
            gravitation[2],
            drivenCavityVelocity[0],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2]);
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
    lbm_kernel_beta<T><<<getBlocksPerGrid(3, domain.getSize(), threadsPerBlock[2]), threadsPerBlock[2]>>>(
            densityDistributions,
            flags,
            velocities,
            densities,
            tauInv,
            gravitation[0],
            gravitation[1],
            gravitation[2],
            drivenCavityVelocity[0],
            domain.getSize()[0],
            domain.getSize()[1],
            domain.getSize()[2],
            threadsPerBlock[2].x * threadsPerBlock[2].y * threadsPerBlock[2].z,
            isPowerOfTwo(domain.getNumOfCells()),
            isPowerOfTwo(threadsPerBlock[2].x * threadsPerBlock[2].y * threadsPerBlock[2].z));
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBetaRect(CVector<3, int> origin, CVector<3, int> size)
{
    /*
     * TODO
     */
}

template <class T>
void CLbmSolverGPU<T>::reset()
{
    /*
     * TODO
     */
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(Direction direction, T* hDensityDistributions)
{
    assert(0 <= direction < 6);

    CVector<3, int> origin(0);
    CVector<3, int> size(1);
    CVector<3, int> norm(0);

    switch(direction)
    {
    case LEFT:
        size[1] = domain.getSize()[1];
        size[2] = domain.getSize()[2];
        norm[0] = -1;
    case RIGHT:
        size[1] = domain.getSize()[1];
        size[2] = domain.getSize()[2];
        origin[0] = domain.getSize()[0] - 1;
        norm[0] = 1;
    case BOTTOM:
        size[0] = domain.getSize()[0];
        size[2] = domain.getSize()[2];
        norm[1] = -1;
    case TOP:
        size[0] = domain.getSize()[0];
        size[2] = domain.getSize()[2];
        origin[0] = domain.getSize()[1] - 1;
        norm[1] = 1;
    case BACK:
        size[0] = domain.getSize()[0];
        size[1] = domain.getSize()[1];
        norm[2] = -1;
    case FRONT:
        size[0] = domain.getSize()[0];
        size[1] = domain.getSize()[1];
        origin[0] = domain.getSize()[2] - 1;
        norm[2] = 1;
    }

    dim3 threadsPerBlock(size[1], size[2], 1);

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        if(norm.dotProd(lbm_units[latticeVector]) > 0)
        {
            copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
                    densityDistributions,
                    latticeVector * domain.getNumOfCells(),
                    origin[0],
                    origin[1],
                    origin[2],
                    domain.getSize()[0],
                    domain.getSize()[1],
                    domain.getSize()[2],
                    getDensityDistributionsHalo[direction],
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
        }
        GPU_ERROR_CHECK(cudaPeekAtLastError())
    }

    GPU_ERROR_CHECK(cudaMemcpy(hDensityDistributions, getDensityDistributionsHalo[direction], NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyDeviceToHost))

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation from device to host for lattice vectors in direction " << direction << " was performed." << std::endl;
        std::cout << "No additional buffer memory was allocated. Instead, getDensityDistributionsHalo was used." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
        std::cout << "direction:     " << norm << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    T* dDensityDistributions;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T)))

    for(int dim = 0; dim < NUM_LATTICE_VECTORS; dim++)
    {
        copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
                densityDistributions,
                dim * domain.getNumOfCells(),
                origin[0],
                origin[1],
                origin[2],
                domain.getSize()[0],
                domain.getSize()[1],
                domain.getSize()[2],
                dDensityDistributions,
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

    GPU_ERROR_CHECK(cudaMemcpy(hDensityDistributions, dDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyDeviceToHost))

    GPU_ERROR_CHECK(cudaFree(dDensityDistributions))

    if (doLogging)
    {
        std::cout << "----- CLbmSolverGPU<T>::getDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
void CLbmSolverGPU<T>::setDensityDistributions(Direction direction, T* hDensityDistributions)
{
    assert(0 <= direction < 6);

    CVector<3, int> origin(0);
    CVector<3, int> size(1);
    CVector<3, int> norm(0);

    switch(direction)
    {
    case LEFT:
        size[1] = domain.getSize()[1];
        size[2] = domain.getSize()[2];
        norm[0] = -1;
    case RIGHT:
        size[1] = domain.getSize()[1];
        size[2] = domain.getSize()[2];
        origin[0] = domain.getSize()[0] - 1;
        norm[0] = 1;
    case BOTTOM:
        size[0] = domain.getSize()[0];
        size[2] = domain.getSize()[2];
        norm[1] = -1;
    case TOP:
        size[0] = domain.getSize()[0];
        size[2] = domain.getSize()[2];
        origin[0] = domain.getSize()[1] - 1;
        norm[1] = 1;
    case BACK:
        size[0] = domain.getSize()[0];
        size[1] = domain.getSize()[1];
        norm[2] = -1;
    case FRONT:
        size[0] = domain.getSize()[0];
        size[1] = domain.getSize()[1];
        origin[0] = domain.getSize()[2] - 1;
        norm[2] = 1;
    }

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMemcpy(getDensityDistributionsHalo[direction], hDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        if(norm.dotProd(lbm_units[latticeVector]) > 0)
        {
            copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
                    getDensityDistributionsHalo[direction],
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
        std::cout << "----- CLbmSolverGPU<T>::setDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation from host to device for lattice vectors in direction " << direction << " was performed." << std::endl;
        std::cout << "No additional buffer memory was allocated. Instead, setDensityDistributionsHalo was used." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
        std::cout << "direction:     " << norm << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverGPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* hDensityDistributions)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    T* dDensityDistributions;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T)))

    GPU_ERROR_CHECK(cudaMemcpy(dDensityDistributions, hDensityDistributions, NUM_LATTICE_VECTORS * size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    for (int i = 0; i < NUM_LATTICE_VECTORS; i++)
    {
        copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
                dDensityDistributions,
                i * size.elements(),
                0,
                0,
                0,
                size[0],
                size[1],
                size[2],
                densityDistributions,
                i * domain.getNumOfCells(),
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
        std::cout << "----- CLbmSolverGPU<T>::setDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    Flag* dFlags;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dFlags, size.elements() * sizeof(Flag)))

    copy_buffer_rect<Flag><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
        std::cout << "----- CLbmSolverGPU<T>::getFlags() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    Flag* dFlags;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dFlags, size.elements() * sizeof(Flag)))

    GPU_ERROR_CHECK(cudaMemcpy(dFlags, hFlags, size.elements() * sizeof(Flag), cudaMemcpyHostToDevice))

    copy_buffer_rect<Flag><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
        std::cout << "----- CLbmSolverGPU<T>::setFlags() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    T* dVelocities;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dVelocities, 3 * size.elements() * sizeof(T)))

    for (int dim = 0; dim < 3; dim++)
    {
        copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
        std::cout << "----- CLbmSolverGPU<T>::getVelocities() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    T* dVelocities;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dVelocities, 3 * size.elements() * sizeof(T)))

    GPU_ERROR_CHECK(cudaMemcpy(dVelocities, hVelocities, 3 * size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    for (int dim = 0; dim < 3; dim++)
    {
        copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
        std::cout << "----- CLbmSolverGPU<T>::setVelocities() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    T* dDensities;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dDensities, size.elements() * sizeof(T)))

    copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
        std::cout << "----- CLbmSolverGPU<T>::getDensities() -----" << std::endl;
        std::cout << "A copy operation from device to host was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
    assert(origin[0] + size[0] <= domain.getSize()[0] - 1);
    assert(origin[1] + size[1] <= domain.getSize()[1] - 1);
    assert(origin[2] + size[2] <= domain.getSize()[2] - 1);

    T* dDensities;

    dim3 threadsPerBlock(size[1], size[2], 1);

    GPU_ERROR_CHECK(cudaMalloc(&dDensities, size.elements() * sizeof(T)))

    GPU_ERROR_CHECK(cudaMemcpy(dDensities, hDensities, size.elements() * sizeof(T), cudaMemcpyHostToDevice))

    copy_buffer_rect<T><<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
        std::cout << "----- CLbmSolverGPU<T>::setDensities() -----" << std::endl;
        std::cout << "A copy operation from host to device was performed." << std::endl;
        std::cout << "Additional buffer memory was allocated and freed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "id:            " << id << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "domain origin: " << domain.getOrigin() << std::endl;
        std::cout << "domain size:   " << domain.getSize() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "cuboid origin: " << origin << std::endl;
        std::cout << "cuboid size:   " << size << std::endl;
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
