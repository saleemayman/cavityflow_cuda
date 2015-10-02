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

#include <cuda.h>

#include "gpukernels/copy_buffer_rect.cuh"
#include "gpukernels/lbm_alpha.cuh"
#include "gpukernels/lbm_beta.cuh"
#include "gpukernels/lbm_init.cuh"

template <class T>
CLbmSolverGPU<T>::CLbmSolverGPU(
		int id,
		CDomain<T> &domain,
		std::array<Flag,6> boundaryConditions,
		T timestepSize,
		CVector<3,T> &gravitation,
		CVector<4,T> &drivenCavityVelocity,
		T viscocity,
		T massExchangeFactor,
		T maxSimGravitationLength,
		T tau,
		bool storeDensities,
		bool storeVelocities) :
		CLbmSolver<T>(id, domain,
				boundaryConditions, timestepSize,
				gravitation, drivenCavityVelocity, viscocity,
				massExchangeFactor, maxSimGravitationLength, tau,
				storeDensities, storeVelocities)
{
	GPU_ERROR_CHECK(cudaMalloc(&densityDistributions, this->domain.getNumOfCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc(&flags, this->domain.getNumOfCellsWithHalo() * sizeof(Flag)))
	if(this->storeDensities)
		GPU_ERROR_CHECK(cudaMalloc(&densities, this->domain.getNumOfCellsWithHalo() * 3 * sizeof(T)))
	if(this->storeVelocities)
		GPU_ERROR_CHECK(cudaMalloc(&velocities, this->domain.getNumOfCellsWithHalo() * sizeof(T)))

	GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[0]), this->domain.getNumOfXFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[0]), this->domain.getNumOfXFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[1]), this->domain.getNumOfXFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[1]), this->domain.getNumOfXFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[2]), this->domain.getNumOfYFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[2]), this->domain.getNumOfYFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[3]), this->domain.getNumOfYFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[3]), this->domain.getNumOfYFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[4]), this->domain.getNumOfZFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[4]), this->domain.getNumOfZFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(getDensityDistributionsHalo[5]), this->domain.getNumOfZFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))
	GPU_ERROR_CHECK(cudaMalloc<T>(&(setDensityDistributionsHalo[5]), this->domain.getNumOfZFaceCellsWithHalo() * NUM_LATTICE_VECTORS * sizeof(T)))

	lbm_init<<<getBlocksPerGrid(3, this->domain.getSize(), this->threadsPerBlock[0]), this->threadsPerBlock[0]>>>(
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
		domain.getSizeWithHalo()[0],
		domain.getSizeWithHalo()[1],
		domain.getSizeWithHalo()[2]);
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

	if(this->storeVelocities)
		GPU_ERROR_CHECK(cudaFree(velocities))
	if(this->storeDensities)
		GPU_ERROR_CHECK(cudaFree(densities))
	GPU_ERROR_CHECK(cudaFree(flags))
	GPU_ERROR_CHECK(cudaFree(densityDistributions))
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlpha()
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepAlphaRect(CVector<3,int> origin, CVector<3,int> size)
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBeta()
{
}

template <class T>
void CLbmSolverGPU<T>::simulationStepBetaRect(CVector<3,int> origin, CVector<3,int> size)
{
}

template <class T>
void CLbmSolverGPU<T>::reset()
{
}

template <class T>
void CLbmSolverGPU<T>::reload()
{
}

template <class T>
void CLbmSolverGPU<T>::getDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setDesityDistribution(CVector<3,int> &origin, CVector<3,int> &size, int i, T* src)
{
}

template <class T>
void CLbmSolverGPU<T>::getFlags(CVector<3,int> &origin, CVector<3,int> &size, int* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setFlags(CVector<3,int> &origin, CVector<3,int> &size, int* src)
{
}

template <class T>
void CLbmSolverGPU<T>::getVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* dst)
{
}

template <class T>
void CLbmSolverGPU<T>::setVelocities(CVector<3,int> &origin, CVector<3,int> &size, T* src)
{
}

template <class T>
void CLbmSolverGPU<T>::getDensities(CVector<3,int> &origin, CVector<3,int> &size, T* hDensities)
{
	assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
	assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
	assert(origin[0] + size[0] <= this->domain.getSize()[0] - 1);
	assert(origin[1] + size[1] <= this->domain.getSize()[1] - 1);
	assert(origin[2] + size[2] <= this->domain.getSize()[2] - 1);

	T* dDensities;

	dim3 threadsPerBlock(size[1], size[2], 1);

	GPU_ERROR_CHECK(cudaMalloc(&dDensities, size.elements() * sizeof(T)))

	copy_buffer_rect<<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
	        densities,
	        0,
	        origin[0],
	        origin[1],
	        origin[2],
	        this->domain.getSize()[0],
	        this->domain.getSize()[1],
	        this->domain.getSize()[2],
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
}

template <class T>
void CLbmSolverGPU<T>::setDensities(CVector<3,int> &origin, CVector<3,int> &size, T* hDensities)
{
	assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
	assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
	assert(origin[0] + size[0] <= this->domain.getSize()[0] - 1);
	assert(origin[1] + size[1] <= this->domain.getSize()[1] - 1);
	assert(origin[2] + size[2] <= this->domain.getSize()[2] - 1);

	T* dDensities;

	dim3 threadsPerBlock(size[1], size[2], 1);

	GPU_ERROR_CHECK(cudaMalloc(&dDensities, size.elements() * sizeof(T)))

	GPU_ERROR_CHECK(cudaMemcpy(dDensities, hDensities, size.elements() * sizeof(T), cudaMemcpyHostToDevice))

	copy_buffer_rect<<<getBlocksPerGrid(2, size, threadsPerBlock), threadsPerBlock>>>(
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
	        this->domain.getSize()[0],
	        this->domain.getSize()[1],
	        this->domain.getSize()[2],
	        size[0],
	        size[1],
	        size[2]);
	GPU_ERROR_CHECK(cudaPeekAtLastError())

	GPU_ERROR_CHECK(cudaFree(dDensities))
}

template class CLbmSolverGPU<float>;
template class CLbmSolverGPU<double>;
