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

#include <cassert>
#include <cstring>

#include "CLbmSolverCPU.hpp"

template <class T>
CLbmSolverCPU<T>::CLbmSolverCPU(
        int id,
        CVector<3, T> &globalLength,
        CDomain<T> &domain,
        std::vector<Flag> boundaryConditions,
        CLbmSolverGPU<T>* solverGPU,
        T timestepSize,
        CVector<3, T>& velocity,
        CVector<3, T>& acceleration,
        T viscocity,
        T maxVelocityDimLess,
        T maxAccelerationDimLess,
        bool storeDensities,
        bool storeVelocities,
        bool doLogging) :
        CLbmSolver<T>(id, globalLength,
                domain, boundaryConditions,
                timestepSize, velocity, acceleration,
                viscocity, maxVelocityDimLess, maxAccelerationDimLess,
                storeDensities, storeVelocities, doLogging),
        solverGPU(solverGPU)
{
    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::CLbmSolverCPU() -----" << std::endl;
        std::cout << "id:                                          " << this->id << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    domainSizeGPU = (solverGPU->getDomain())->getSize();
#if TOP_DOWN_DECOMPOSITION
	domainSizeCPUWithHalo.set(this->domain.getSize()[0] + 2, this->domain.getSize()[1] - domainSizeGPU[1] + 2, this->domain.getSize()[2] + 2);
#endif

#if !TOP_DOWN_DECOMPOSITION
    /*
     * Set the limits of the CPU domain w.r.t to the inner GPU domain.
     */
    hollowCPULeftLimit[0] = (this->domain.getSize()[0] - domainSizeGPU[0]) / 2 - 1;
    hollowCPURightLimit[0] = hollowCPULeftLimit[0] + domainSizeGPU[0] + 1;
    hollowCPULeftLimit[1] = (this->domain.getSize()[1] - domainSizeGPU[1]) / 2 - 1;
    hollowCPURightLimit[1] = hollowCPULeftLimit[1] + domainSizeGPU[1] + 1;
    hollowCPULeftLimit[2] = (this->domain.getSize()[2] - domainSizeGPU[2]) / 2 - 1;
    hollowCPURightLimit[2] = hollowCPULeftLimit[2] + domainSizeGPU[2] + 1;
#endif

    /*
     * Allocate memory for density distributions, density, flags and velocities.
     */
#if !TOP_DOWN_DECOMPOSITION
    domainCellsCPUWithHalo = this->domain.getNumOfCellsWithHalo() - ((domainSizeGPU[0] - 2) * (domainSizeGPU[1] - 2) * (domainSizeGPU[2] - 2));
#else
    domainCellsCPUWithHalo = domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1] * domainSizeCPUWithHalo[2];
#endif

    densityDistributions.resize(NUM_LATTICE_VECTORS * domainCellsCPUWithHalo); 
    flags.resize(domainCellsCPUWithHalo);

    if (storeDensities)
        densities.resize(domainCellsCPUWithHalo);
    if (storeVelocities)
        velocities.resize(3 * domainCellsCPUWithHalo);

    if (doLogging) {
        std::cout << "size of allocated memory for density distributions: " << ((T)(NUM_LATTICE_VECTORS * domainCellsCPUWithHalo * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        std::cout << "size of allocated memory for flags:                 " << ((T)(domainCellsCPUWithHalo * sizeof(Flag)) / (T)(1<<20)) << " MBytes" << std::endl;
        if(this->storeDensities)
            std::cout << "size of allocated memory for velocities:            " << ((T)(3 * domainCellsCPUWithHalo * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
        if(this->storeVelocities)
            std::cout << "size of allocated memory for densities:             " << ((T)(domainCellsCPUWithHalo * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
    }

    /*
     * instantiate the cpu init-kernel class and initialize the LBM simulation
     */
#if !TOP_DOWN_DECOMPOSITION
    initLbmCPU = new CLbmInitCPU<T>(
                        domainCellsCPUWithHalo,
                        this->domain.getSizeWithHalo(), 
                        hollowCPULeftLimit, 
                        hollowCPURightLimit, 
                        boundaryConditions);
#else
    initLbmCPU = new CLbmInitCPU<T>(
                        domainCellsCPUWithHalo,
                        domainSizeCPUWithHalo, 
                        boundaryConditions);
#endif
    initLbmCPU->initLbm(
                    densityDistributions, 
                    flags, 
                    velocities, 
                    densities, 
                    storeDensities, 
                    storeVelocities);

    if (doLogging) {
        std::cout << "CPU Domain successfully initialized!" <<  std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    /*
     * instantiate the alhpa and beta classes
     */
#if !TOP_DOWN_DECOMPOSITION
    alphaLbmCPU = new CLbmAlphaCPU<T>(
                            domainCellsCPUWithHalo,
                            this->domain.getSizeWithHalo(), 
                            //initLbmCPU->getCellIndexMap(),
                            initLbmCPU,
                            acceleration);
    betaLbmCPU = new CLbmBetaCPU<T>(
                            domainCellsCPUWithHalo,
                            this->domain.getSizeWithHalo(), 
                            //initLbmCPU->getCellIndexMap(),
                            initLbmCPU,
                            acceleration);
#else
    alphaLbmCPU = new CLbmAlphaCPU<T>(
                            domainCellsCPUWithHalo,
                            domainSizeCPUWithHalo, 
                            //initLbmCPU->getCellIndexMap(),
                            initLbmCPU,
                            acceleration);
    betaLbmCPU = new CLbmBetaCPU<T>(
                            domainCellsCPUWithHalo,
                            domainSizeCPUWithHalo, 
                            //initLbmCPU->getCellIndexMap(),
                            initLbmCPU,
                            acceleration);
#endif
}

template <class T>
CLbmSolverCPU<T>::~CLbmSolverCPU()
{
    delete initLbmCPU;
    delete alphaLbmCPU;
    delete betaLbmCPU;
}

template <class T>
void CLbmSolverCPU<T>::getVariable(CVector<3, int> &origin, CVector<3, int> &size, std::vector<T> &variableData, T* dst, int numDimensions)
{
    // get global linear index of starting cell
    int dstIndex;
    int srcIndex;
    int localIndex; 
    int globalIndex;
#if !TOP_DOWN_DECOMPOSITION
    int indexOffset = origin[0] + origin[1] * domain.getSizeWithHalo()[0] + origin[2] * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
#else
    int indexOffset = origin[0] + origin[1] * domainSizeCPUWithHalo[0] + origin[2] * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
#endif

    /* 
     * loop over the number of cell widths to copy in z- and y-dimensions, 
     * and then do a memcpy of all the size[0] cells in the x-dimension.
     */
    for(int dim = 0; dim < numDimensions; dim++)
    {
        for (int widthZ = 0; widthZ < size[2]; widthZ++)
        {
            for (int widthY = 0; widthY < size[1]; widthY++)
            {
#if !TOP_DOWN_DECOMPOSITION
                globalIndex = indexOffset + widthY * domain.getSizeWithHalo()[0] + widthZ * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
                localIndex = initLbmCPU->getLocalIndex(globalIndex);
#else
                globalIndex = indexOffset + widthY * domainSizeCPUWithHalo[0] + widthZ * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
                localIndex = globalIndex;
#endif                
                dstIndex = widthY * size[0] + widthZ * (size[0] * size[1]) + dim * size.elements();
                srcIndex = localIndex + dim * domainCellsCPUWithHalo;

                std::memcpy(dst + dstIndex, variableData.data() + srcIndex, size[0] * sizeof(T));
            }
        }
    } 
}

template <class T>
void CLbmSolverCPU<T>::getVariable(CVector<3, int> &origin, CVector<3, int> &size, std::vector<Flag> &variableData, Flag* dst, int numDimensions)
{
    // get global linear index of starting cell
    int dstIndex;
    int srcIndex;
    int localIndex; 
    int globalIndex;
#if !TOP_DOWN_DECOMPOSITION
    int indexOffset = origin[0] + origin[1] * domain.getSizeWithHalo()[0] + origin[2] * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
#else
    int indexOffset = origin[0] + origin[1] * domainSizeCPUWithHalo[0] + origin[2] * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
#endif
    /* 
     * loop over the number of cell widths to copy in z- and y-dimensions, 
     * and then do a memcpy of all the size[0] cells in the x-dimension.
     */
    for(int dim = 0; dim < numDimensions; dim++)
    {
        for (int widthZ = 0; widthZ < size[2]; widthZ++)
        {
            for (int widthY = 0; widthY < size[1]; widthY++)
            {
#if !TOP_DOWN_DECOMPOSITION
                globalIndex = indexOffset + widthY * domain.getSizeWithHalo()[0] + widthZ * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
                localIndex = initLbmCPU->getLocalIndex(globalIndex);
#else
                globalIndex = indexOffset + widthY * domainSizeCPUWithHalo[0] + widthZ * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
                localIndex = globalIndex;
#endif                
                dstIndex = widthY * size[0] + widthZ * (size[0] * size[1]) + dim * size.elements();
                srcIndex = localIndex + dim * domainCellsCPUWithHalo;

                std::memcpy(dst + dstIndex, variableData.data() + srcIndex, size[0] * sizeof(Flag));
            }
        }
    } 
}

template <class T>
void CLbmSolverCPU<T>::setVariable(CVector<3, int> &origin, CVector<3, int> &size, std::vector<T> &variableData, T* src, int numDimensions)
{
    // get global linear index of starting cell
    int dstIndex;
    int srcIndex;
    int localIndex; 
    int globalIndex;
#if !TOP_DOWN_DECOMPOSITION
    int indexOffset = origin[0] + origin[1] * domain.getSizeWithHalo()[0] + origin[2] * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
#else
    int indexOffset = origin[0] + origin[1] * domainSizeCPUWithHalo[0] + origin[2] * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
#endif
    /* 
     * loop over the number of cell widths to copy in z- and y-dimensions, 
     * and then do a memcpy of all the size[0] cells in the x-dimension.
     */
    for(int dim = 0; dim < numDimensions; dim++)
    {
        for (int widthZ = 0; widthZ < size[2]; widthZ++)
        {
            for (int widthY = 0; widthY < size[1]; widthY++)
            {
#if !TOP_DOWN_DECOMPOSITION
                globalIndex = indexOffset + widthY * domain.getSizeWithHalo()[0] + widthZ * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
                localIndex = initLbmCPU->getLocalIndex(globalIndex);
#else
                globalIndex = indexOffset + widthY * domainSizeCPUWithHalo[0] + widthZ * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
                localIndex = globalIndex;
#endif                
                srcIndex = widthY * size[0] + widthZ * (size[0] * size[1]) + dim * size.elements();
                dstIndex = localIndex + dim * domainCellsCPUWithHalo; 

                std::memcpy(variableData.data() + dstIndex, src + srcIndex, size[0] * sizeof(T));
            }
        }
    } 
}

template <class T>
void CLbmSolverCPU<T>::setVariable(CVector<3, int> &origin, CVector<3, int> &size, std::vector<Flag> &variableData, Flag* src, int numDimensions)
{
    // get global linear index of starting cell
    int dstIndex;
    int srcIndex;
    int localIndex; 
    int globalIndex;
#if !TOP_DOWN_DECOMPOSITION
    int indexOffset = origin[0] + origin[1] * domain.getSizeWithHalo()[0] + origin[2] * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
#else
    int indexOffset = origin[0] + origin[1] * domainSizeCPUWithHalo[0] + origin[2] * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
#endif

    /* 
     * loop over the number of cell widths to copy in z- and y-dimensions, 
     * and then do a memcpy of all the size[0] cells in the x-dimension.
     */
    for(int dim = 0; dim < numDimensions; dim++)
    {
        for (int widthZ = 0; widthZ < size[2]; widthZ++)
        {
            for (int widthY = 0; widthY < size[1]; widthY++)
            {
#if !TOP_DOWN_DECOMPOSITION
                globalIndex = indexOffset + widthY * domain.getSizeWithHalo()[0] + widthZ * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
                localIndex = initLbmCPU->getLocalIndex(globalIndex);
#else
                globalIndex = indexOffset + widthY * domainSizeCPUWithHalo[0] + widthZ * (domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1]);
                localIndex = globalIndex;
#endif                

                srcIndex = widthY * size[0] + widthZ * (size[0] * size[1]) + dim * size.elements();
                dstIndex = localIndex + dim * domainCellsCPUWithHalo; 

                std::memcpy(variableData.data() + dstIndex, src + srcIndex, size[0] * sizeof(Flag));
            }
        }
    } 
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlpha()
{
    // TODO: possible code duplication. Similar to simulationStepAlpha(origin, size). Check if possible to replace.
    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

    /*
     * launch the alpha-kernel for the CPU domain
     */
#if !TOP_DOWN_DECOMPOSITION
    alphaLbmCPU->alphaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    CVector<3, int>(0, 0, 0),
                    this->domain.getSizeWithHalo(), 
                    storeDensities,
                    storeVelocities);
#else
    alphaLbmCPU->alphaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    CVector<3, int>(0, 0, 0),
                    domainSizeCPUWithHalo, 
                    storeDensities,
                    storeVelocities);
#endif

    if (doLogging)
    {
        std::cout << "Alpha kernel was successfully executed on the whole CPU subdomain." << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlpha(CVector<3, int> origin, CVector<3, int> size)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

    /*
     * launch the alpha-kernel for the CPU domain
     */
    alphaLbmCPU->alphaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    origin,
                    size, 
                    storeDensities,
                    storeVelocities);

    if (doLogging)
    {
        std::cout << "Alpha kernel was successfully executed on the following CPU subdomain:" << std::endl;
        std::cout << "origin:            " << origin << std::endl;
        std::cout << "size:              " << size << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBeta()
{
    // TODO: possible code duplication. Similar to simulationStepBeta(origin, size). Check if possible to replace.
    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
        std::cout << "id:                 " << id << std::endl;
        std::cout << "---------------------------------------------------" << std::endl;
    }

#if !TOP_DOWN_DECOMPOSITION
    betaLbmCPU->betaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    CVector<3, int>(0, 0, 0),
                    this->domain.getSizeWithHalo(), 
                    isPowerOfTwo(this->domain.getNumOfCellsWithHalo()),
                    storeDensities,
                    storeVelocities);
#else
    betaLbmCPU->betaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    CVector<3, int>(0, 0, 0),
                    domainSizeCPUWithHalo, 
                    isPowerOfTwo(domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1] * domainSizeCPUWithHalo[2]),
                    storeDensities,
                    storeVelocities);
#endif

    if (doLogging)
    {
        std::cout << "Beta kernel was successfully executed on the whole CPU subdomain." << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBeta(CVector<3, int> origin, CVector<3, int> size)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
        std::cout << "id:                 " << id << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }

#if !TOP_DOWN_DECOMPOSITION
    betaLbmCPU->betaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    origin,
                    size, 
                    isPowerOfTwo(this->domain.getNumOfCellsWithHalo()),
                    storeDensities,
                    storeVelocities);
#else
    betaLbmCPU->betaKernelCPU(
                    densityDistributions,
                    flags,
                    velocities,
                    densities,
                    tauInv,
                    velocityDimLess[0],
                    origin,
                    size, 
                    isPowerOfTwo(domainSizeCPUWithHalo[0] * domainSizeCPUWithHalo[1] * domainSizeCPUWithHalo[2]),
                    storeDensities,
                    storeVelocities);
#endif

    if (doLogging)
    {
        std::cout << "Beta kernel was successfully executed on the following CPU subdomain." << std::endl;
        std::cout << "origin:             " << origin << std::endl;
        std::cout << "size:               " << size << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::getDensityDistributions() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
    
	getVariable(origin, size, densityDistributions, dst, NUM_LATTICE_VECTORS);

    if (doLogging)
    {
        std::cout << "A copy operation with in the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::getDensityDistributions(T* dst)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getDensityDistributions(origin, size, dst);
}

template <class T>
void CLbmSolverCPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
        std::cout << "A copy operation within the host will be performed." << std::endl;
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

	setVariable(origin, size, densityDistributions, src, NUM_LATTICE_VECTORS);

    if (doLogging)
    {
        std::cout << "A copy operation with in the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::setDensityDistributions(T* src)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setDensityDistributions(origin, size, src);
}

template <class T>
void CLbmSolverCPU<T>::getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::getFlags() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

	getVariable(origin, size, flags, dst, 1);

    if (doLogging)
    {
        std::cout << "A copy operation within the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::getFlags(Flag* dst)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getFlags(origin, size, dst);
}

template <class T>
void CLbmSolverCPU<T>::setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* src)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::setFlags() -----" << std::endl;
        std::cout << "A copy operation within the host will be performed." << std::endl;
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

	setVariable(origin, size, flags, src, 1);

    if (doLogging)
    {
        std::cout << "A copy operation within the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::setFlags(Flag* src)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setFlags(origin, size, src);
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::getVelocities() -----" << std::endl;
        std::cout << "id:                " << id << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "domain origin:     " << domain.getOrigin() << std::endl;
        std::cout << "domain size:       " << domain.getSize() << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "cuboid origin:     " << origin << std::endl;
        std::cout << "cuboid size:       " << size << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }

	getVariable(origin, size, velocities, dst, 3);

    if (doLogging)
    {
        std::cout << "A copy operation within the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(T* dst)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getVelocities(origin, size, dst);
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::setVelocities() -----" << std::endl;
        std::cout << "A copy operation within the host will be performed." << std::endl;
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

	setVariable(origin, size, velocities, src, 3);

    if (doLogging)
    {
        std::cout << "A copy operation with in the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(T* src)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setVelocities(origin, size, src);
}

template <class T>
void CLbmSolverCPU<T>::getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::getDensities() -----" << std::endl;
        std::cout << "A copy operation within the host will be performed." << std::endl;
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

	getVariable(origin, size, densities, dst, 1);

    if (doLogging)
    {
        std::cout << "A copy operation within the host was performed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::getDensities(T* dst)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    getDensities(origin, size, dst);
}

template <class T>
void CLbmSolverCPU<T>::setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
#if !TOP_DOWN_DECOMPOSITION
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
	assert(!(origin[0] > hollowCPULeftLimit[0] && origin[0] < hollowCPURightLimit[0] && origin[1] > hollowCPULeftLimit[1] && origin[1] < hollowCPURightLimit[1] && origin[2] > hollowCPULeftLimit[2] && origin[2] < hollowCPURightLimit[2]));
#else
    assert(origin[1] + size[1] <= domainSizeCPUWithHalo[1]);
#endif

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::setDensities() -----" << std::endl;
        std::cout << "A copy operation within the host will be performed." << std::endl;
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

	setVariable(origin, size, densities, src, 1);

    if (doLogging)
    {
        std::cout << "A copy operation within the host was performed." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::setDensities(T* src)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setDensities(origin, size, src);
}

#if !TOP_DOWN_DECOMPOSITION
template <class T>
CVector<3, int> CLbmSolverCPU<T>::getHollowCPULeftLimits()
{
    return hollowCPULeftLimit;
}

template <class T>
CVector<3, int> CLbmSolverCPU<T>::getHollowCPURightLimits()
{
    return hollowCPURightLimit;
}
#endif

template class CLbmSolverCPU<float>;
template class CLbmSolverCPU<double>;
