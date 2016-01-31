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
#include <fstream>
#include <sstream>

#include "CLbmSolverCPU.hpp"

template <class T>
CLbmSolverCPU<T>::CLbmSolverCPU(
        int id,
        CDomain<T> &domain,
        CLbmSolverGPU<T>* solverGPU,
        std::vector<Flag> boundaryConditions,
        CConfiguration<T>* configuration) :
        CLbmSolver<T>(id,
                domain,
                boundaryConditions,
                configuration),
        solverGPU(solverGPU)
{
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::CLbmSolverCPU() -----" << std::endl;
            loggingFile << "id:                                                 " << this->id << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::CLbmSolverCPU() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    densityDistributions = new T[NUM_LATTICE_VECTORS * this->domain.getNumOfCellsWithHalo()];
    flags = new Flag[this->domain.getNumOfCellsWithHalo()];
    if (storeVelocities)
        velocities = new T[3 * this->domain.getNumOfCellsWithHalo()];
    if (storeDensities)
        densities = new T[this->domain.getNumOfCellsWithHalo()];

    if (configuration->doLogging) {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "size of allocated memory for density distributions: " << ((T)(NUM_LATTICE_VECTORS * this->domain.getNumOfCellsWithHalo() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
            loggingFile << "size of allocated memory for flags:                 " << ((T)(this->domain.getNumOfCellsWithHalo() * sizeof(Flag)) / (T)(1<<20)) << " MBytes" << std::endl;
            if (storeVelocities)
                loggingFile << "size of allocated memory for velocities:            " << ((T)(3 * this->domain.getNumOfCellsWithHalo() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
            if (storeDensities)
                loggingFile << "size of allocated memory for densities:             " << ((T)(this->domain.getNumOfCellsWithHalo() * sizeof(T)) / (T)(1<<20)) << " MBytes" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::CLbmSolverCPU() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    initLbmCPU = new CLbmInitCPU<T>(
            domain.getSizeWithHalo(),
            boundaryConditions);
    initLbmCPU->initLbm(
            densityDistributions,
            flags,
            densities,
            velocities,
            storeDensities,
            storeVelocities);
    alphaLbmCPU = new CLbmAlphaCPU<T>(
            domain.getSizeWithHalo(),
            accelerationDimLess);
    betaLbmCPU = new CLbmBetaCPU<T>(
            domain.getSizeWithHalo(),
            accelerationDimLess);

    if (configuration->doLogging) {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "CPU domain successfully initialized." << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::CLbmSolverCPU() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }
}

template <class T>
CLbmSolverCPU<T>::~CLbmSolverCPU()
{
    delete betaLbmCPU;
    delete alphaLbmCPU;
    delete initLbmCPU;

    if(storeDensities)
        delete[] densities;
    if(storeVelocities)
        delete[] velocities;
    delete[] flags;
    delete[] densityDistributions;
}

template <class T>
void CLbmSolverCPU<T>::getVariable(CVector<3, int> &origin, CVector<3, int> &size, T* src, T* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    // domain location and size
    params.srcPtr = make_cudaPitchedPtr(src, domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.srcPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
    // destination location and size
    params.dstPtr = make_cudaPitchedPtr(dst, size[0] * sizeof(T), size[0], size[1]);
    // destination origin
    params.dstPos = make_cudaPos(0, 0, 0);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyHostToHost;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))
}

template <class T>
void CLbmSolverCPU<T>::getVariable(CVector<3, int> &origin, CVector<3, int> &size, Flag* src, Flag* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    // domain location and size
    params.srcPtr = make_cudaPitchedPtr(src, domain.getSizeWithHalo()[0] * sizeof(Flag), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.srcPos = make_cudaPos(origin[0] * (sizeof(Flag) / sizeof(unsigned char)), origin[1], origin[2]);
    // destination location and size
    params.dstPtr = make_cudaPitchedPtr(dst, size[0] * sizeof(Flag), size[0], size[1]);
    // destination origin
    params.dstPos = make_cudaPos(0, 0, 0);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyHostToHost;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))
}

template <class T>
void CLbmSolverCPU<T>::setVariable(CVector<3, int> &origin, CVector<3, int> &size, T* src, T* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    // source location and size
    params.srcPtr = make_cudaPitchedPtr(src, size[0] * sizeof(T), size[0], size[1]);
    // source origin
    params.srcPos = make_cudaPos(0, 0, 0);
    // domain location and size
    params.dstPtr = make_cudaPitchedPtr(dst, domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.dstPos = make_cudaPos(origin[0] * (sizeof(T) / sizeof(unsigned char)), origin[1], origin[2]);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(T) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyHostToHost;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))
}

template <class T>
void CLbmSolverCPU<T>::setVariable(CVector<3, int> &origin, CVector<3, int> &size, Flag* src, Flag* dst)
{
    assert(origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0);
    assert(size[0] > 0 && size[1] > 0 && size[2] > 0);
    assert(origin[0] + size[0] <= domain.getSizeWithHalo()[0]);
    assert(origin[1] + size[1] <= domain.getSizeWithHalo()[1]);
    assert(origin[2] + size[2] <= domain.getSizeWithHalo()[2]);

    cudaMemcpy3DParms params = {0};

    // source location and size
    params.srcPtr = make_cudaPitchedPtr(src, size[0] * sizeof(Flag), size[0], size[1]);
    // source origin
    params.srcPos = make_cudaPos(0, 0, 0);
    // domain location and size
    params.dstPtr = make_cudaPitchedPtr(dst, domain.getSizeWithHalo()[0] * sizeof(T), domain.getSizeWithHalo()[0], domain.getSizeWithHalo()[1]);
    // cuboid origin
    params.dstPos = make_cudaPos(origin[0] * (sizeof(Flag) / sizeof(unsigned char)), origin[1], origin[2]);
    // cuboid size
    params.extent = make_cudaExtent(size[0] * (sizeof(Flag) / sizeof(unsigned char)), size[1], size[2]);
    params.kind = cudaMemcpyHostToHost;

    GPU_ERROR_CHECK(cudaMemcpy3D(&params))
}

template <class T>
void CLbmSolverCPU<T>::simulationStepAlpha()
{
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "---------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    alphaLbmCPU->alphaKernelCPU(
                    densityDistributions,
                    flags,
                    densities,
                    velocities,
                    tauInv,
                    velocityDimLess[0],
                    CVector<3, int>(0, 0, 0),
                    domain.getSizeWithHalo(),
                    configuration->elementsPerBlock[0],
                    storeDensities,
                    storeVelocities);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "Alpha kernel was successfully executed on the whole CPU subdomain." << std::endl;
            loggingFile << "---------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverGCU<T>::simulationStepAlpha() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "---------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    alphaLbmCPU->alphaKernelCPU(
            densityDistributions,
            flags,
            densities,
            velocities,
            tauInv,
            velocityDimLess[0],
            origin,
            size,
            configuration->elementsPerBlock[0],
            storeDensities,
            storeVelocities);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "Alpha kernel was successfully executed on the following CPU subdomain:" << std::endl;
            loggingFile << "origin:            " << origin << std::endl;
            loggingFile << "size:              " << size << std::endl;
            loggingFile << "---------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepAlpha() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }
}

template <class T>
void CLbmSolverCPU<T>::simulationStepBeta()
{
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
            loggingFile << "id:                 " << id << std::endl;
            loggingFile << "--------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    betaLbmCPU->betaKernelCPU(
            densityDistributions,
            flags,
            densities,
            velocities,
            tauInv,
            velocityDimLess[0],
            CVector<3, int>(0, 0, 0),
            domain.getSizeWithHalo(),
            configuration->elementsPerBlock[1],
            isPowerOfTwo(domain.getNumOfCellsWithHalo()),
            storeDensities,
            storeVelocities);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "Beta kernel was successfully executed on the whole CPU subdomain." << std::endl;
            loggingFile << "--------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
            loggingFile << "id:                 " << id << std::endl;
            loggingFile << "--------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    betaLbmCPU->betaKernelCPU(
            densityDistributions,
            flags,
            densities,
            velocities,
            tauInv,
            velocityDimLess[0],
            origin,
            size,
            configuration->elementsPerBlock[1],
            isPowerOfTwo(domain.getNumOfCellsWithHalo()),
            storeDensities,
            storeVelocities);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "Beta kernel was successfully executed on the following CPU subdomain." << std::endl;
            loggingFile << "origin:             " << origin << std::endl;
            loggingFile << "size:               " << size << std::endl;
            loggingFile << "--------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::simulationStepBeta() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }
}

template <class T>
void CLbmSolverCPU<T>::getDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* dst)
{
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::getDensityDistributions() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getDensityDistributions() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-------------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    for(int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        getVariable(origin, size, &densityDistributions[latticeVector * domain.getNumOfCellsWithHalo()], &dst[latticeVector * size.elements()]);
    }

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getDensityDistributions() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-------------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
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
void CLbmSolverCPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, Direction direction, T* src)
{
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

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "direction:         " << norm << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-------------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        if(norm.dotProd(lbm_units[latticeVector]) > 0)
            setVariable(origin, size, &src[latticeVector * size.elements()], &densityDistributions[latticeVector * domain.getNumOfCellsWithHalo()]);
    }

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host for lattice vectors in direction " << direction << " was performed." << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-------------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }
}

template <class T>
void CLbmSolverCPU<T>::setDensityDistributions(CVector<3, int> &origin, CVector<3, int> &size, T* src)
{
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-------------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    for (int latticeVector = 0; latticeVector < NUM_LATTICE_VECTORS; latticeVector++)
    {
        setVariable(origin, size, &src[latticeVector * size.elements()], &densityDistributions[latticeVector * domain.getNumOfCellsWithHalo()]);
    }

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "-------------------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setDensityDistributions() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "-------------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::getFlags() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getFlags() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "----------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    getVariable(origin, size, flags, dst);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getFlags() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "----------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::setFlags() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setFlags() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "----------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    setVariable(origin, size, src, flags);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "----------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setFlags() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "----------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::getVelocities() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getVelocities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    for(int i = 0; i < 3; i++)
    {
        getVariable(origin, size, &velocities[i * domain.getNumOfCellsWithHalo()], &dst[i * size.elements()]);
    }

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getVelocities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::setVelocities() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setVelocities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    for (int i = 0; i < NUM_LATTICE_VECTORS; i++)
    {
        setVariable(origin, size, &src[i * size.elements()], &velocities[i * domain.getNumOfCellsWithHalo()]);
    }

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "---------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setVelocities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "---------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::getDensities() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getDensities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    getVariable(origin, size, densities, dst);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::getDensities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
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
    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "----- CLbmSolverCPU<T>::setDensities() -----" << std::endl;
            loggingFile << "id:                " << id << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile << "domain origin:     " << domain.getOrigin() << std::endl;
            loggingFile << "domain size:       " << domain.getSize() << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile << "cuboid origin:     " << origin << std::endl;
            loggingFile << "cuboid size:       " << size << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setDensities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    setVariable(origin, size, src, densities);

    if (configuration->doLogging)
    {
        std::stringstream loggingFileName;
        loggingFileName << configuration->loggingOutputDir << "/log_" << id << ".txt";
        std::ofstream loggingFile(loggingFileName.str().c_str(), std::ios::out | std::ios::app);
        if (loggingFile.is_open())
        {
            loggingFile << "A copy operation from host to host was performed." << std::endl;
            loggingFile << "--------------------------------------------" << std::endl;
            loggingFile.close();
        } else {
            std::cerr << "----- CLbmSolverCPU<T>::setDensities() -----" << std::endl;
            std::cerr << "There is no open file to write logs." << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "--------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }
}

template <class T>
void CLbmSolverCPU<T>::setDensities(T* src)
{
    CVector<3, int> origin(1);
    CVector<3, int> size(domain.getSize());

    setDensities(origin, size, src);
}

template class CLbmSolverCPU<float>;
template class CLbmSolverCPU<double>;
