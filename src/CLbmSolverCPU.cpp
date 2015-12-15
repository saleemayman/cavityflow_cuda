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

    CVector<3, int> domainSizeGPU;
    domainSizeGPU = (solverGPU->getDomain())->getSize();


    /*
     * Set the limits of the CPU domain w.r.t to the inner GPU domain.
     */
    hollowCPULeftLimit[0] = (this->domain.getSize()[0] - domainSizeGPU[0]) / 2 - 1;
    hollowCPURightLimit[0] = hollowCPULeftLimit[0] + domainSizeGPU[0] + 1;
    hollowCPULeftLimit[1] = (this->domain.getSize()[1] - domainSizeGPU[1]) / 2 - 1;
    hollowCPURightLimit[1] = hollowCPULeftLimit[1] + domainSizeGPU[1] + 1;
    hollowCPULeftLimit[2] = (this->domain.getSize()[2] - domainSizeGPU[2]) / 2 - 1;
    hollowCPURightLimit[2] = hollowCPULeftLimit[2] + domainSizeGPU[2] + 1;


    /*
     * Allocate memory for density distributions, density, flags and velocities.
     */
    int domainCellsCPUWithHalo = this->domain.getNumOfCellsWithHalo() -((domainSizeGPU[0] - 2) * (domainSizeGPU[1] - 2) * (domainSizeGPU[2] - 2));

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
    initLbmCPU = new CLbmInitCPU<T>(
                        this->domain.getSize(), 
                        domainSizeGPU, 
                        hollowCPULeftLimit, 
                        hollowCPURightLimit, 
                        boundaryConditions);
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
    alphaLbmCPU = new CLbmAlphaCPU<T>(
                            this->domain.getSize(), 
                            domainSizeGPU, 
                            hollowCPULeftLimit, 
                            hollowCPURightLimit, 
                            acceleration);
    betaLbmCPU = new CLbmBetaCPU<T>(
                            this->domain.getSize(), 
                            domainSizeGPU, 
                            hollowCPULeftLimit, 
                            hollowCPURightLimit, 
                            acceleration, 
                            initLbmCPU->getCellIndexMap());
}

template <class T>
CLbmSolverCPU<T>::~CLbmSolverCPU()
{
    delete initLbmCPU;
    delete alphaLbmCPU;
    delete betaLbmCPU;
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
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
    //TODO: assert(origin[0] is NOT between the left[0] and right[0] CPU inner limits);
    //TODO: assert(origin[1] is NOT between the left[1] and right[1] CPU inner limits);
    //TODO: assert(origin[2] is NOT between the left[2] and right[2] CPU inner limits)
    
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
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
    //TODO: assert(origin[0] is NOT between the left[0] and right[0] CPU inner limits);
    //TODO: assert(origin[1] is NOT between the left[1] and right[1] CPU inner limits);
    //TODO: assert(origin[2] is NOT between the left[2] and right[2] CPU inner limits)

    if (doLogging)
    {
        std::cout << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
        std::cout << "id:                 " << id << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }

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
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
    //TODO: assert(origin[0] is NOT between the left[0] and right[0] CPU inner limits);
    //TODO: assert(origin[1] is NOT between the left[1] and right[1] CPU inner limits);
    //TODO: assert(origin[2] is NOT between the left[2] and right[2] CPU inner limits);

    // get global linear index of starting cell
    int dstIndex;
    int srcIndex;
    int localIndex; 
    int globalIndex;
    int indexOffset = origin[0] + origin[1] * domain.getSizeWithHalo()[0] + origin[2] * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);

    /* 
     * loop over the number of cell widths to copy in z- and y-dimensions, 
     * and then do a memcpy of all the size[0] cells in the x-dimension.
     */
    for(int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        for (int widthZ = 0; widthZ < size[2]; widthZ++)
        {
            for (int widthY = 0; widthY < size[1]; widthY++)
            {
                globalIndex = indexOffset + widthY * domain.getSizeWithHalo()[0] + widthZ * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
                localIndex = initLbmCPU->getLocalIndex(globalIndex);
                
                dstIndex = widthY * size[0] + widthZ * (size[0] * size[1]) + latticeVector * size.elements();
                srcIndex = localIndex + latticeVector * domain.getNumOfCellsWithHalo(); // TODO: domain.getNumOfCellsWithHalo() should not be used. We need the domainCellsCPU

                std::memcpy(dst + dstIndex, densityDistributions.data() + srcIndex, size[0] * sizeof(T));
            }
        }
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
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);
    // TODO: check if origin lies within the hollow region - es muss nicht!
    // TODO: check if (origin+size) overlaps with the hollow region - wieder, es musst nicht sein!

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

    // get global linear index of starting cell
    int dstIndex;
    int srcIndex;
    int localIndex; 
    int globalIndex;
    int indexOffset = origin[0] + origin[1] * domain.getSizeWithHalo()[0] + origin[2] * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);

    /* 
     * loop over the number of cell widths to copy in z- and y-dimensions, 
     * and then do a memcpy of all the size[0] cells in the x-dimension.
     */
    for(int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        for (int widthZ = 0; widthZ < size[2]; widthZ++)
        {
            for (int widthY = 0; widthY < size[1]; widthY++)
            {
                globalIndex = indexOffset + widthY * domain.getSizeWithHalo()[0] + widthZ * (domain.getSizeWithHalo()[0] * domain.getSizeWithHalo()[1]);
                localIndex = initLbmCPU->getLocalIndex(globalIndex);
                
                srcIndex = widthY * size[0] + widthZ * (size[0] * size[1]) + latticeVector * size.elements();
                dstIndex = localIndex + latticeVector * domain.getNumOfCellsWithHalo(); // TODO: we need the domainCellsCPU here.

                std::memcpy(densityDistributions.data() + dstIndex, src + srcIndex, size[0] * sizeof(T));
            }
        }
    } 

    if (doLogging)
    {
        std::cout << "A copy operation with in the host was performed." << std::endl;
        std::cout << "-------------------------------------------------------" << std::endl;
    }
}

template <class T>
void CLbmSolverCPU<T>::setDensityDistributions(T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getFlags(Flag* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setFlags(CVector<3, int> &origin, CVector<3, int> &size, Flag* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setFlags(Flag* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getVelocities(T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setVelocities(T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensities(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::getDensities(T* src)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensities(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
}

template <class T>
void CLbmSolverCPU<T>::setDensities(T* dst)
{
}

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


template class CLbmSolverCPU<float>;
template class CLbmSolverCPU<double>;
