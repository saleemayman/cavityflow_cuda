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

#include "CConfiguration.hpp"

#include <cassert>
#include <limits>

#include "libmath/CMath.hpp"
#include "common.h"

template <class T>
CConfiguration<T>::CConfiguration(std::string fileName)
{
    int status;

    status = doc.LoadFile(fileName.c_str());

    if (status != tinyxml2::XML_NO_ERROR)
    {
        std::cerr << "----- CConfiguration<T>::CConfiguration() -----" << std::endl;
        std::cerr << "Loading XML configuration file \"" << fileName << "\" failed" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "-----------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
    }

    interpretXMLDoc();
    checkParameters();
}

template <class T>
CConfiguration<T>::~CConfiguration()
{
}

template <class T>
void CConfiguration<T>::interpretPhysiscsData(const tinyxml2::XMLNode* root)
{
    const tinyxml2::XMLNode* physicsChild = root->FirstChildElement(TAG_NAME_CHILD_PHYSICS);

    velocity[0] = atof(physicsChild->FirstChildElement("velocity")->FirstChildElement("x")->GetText());
    velocity[1] = atof(physicsChild->FirstChildElement("velocity")->FirstChildElement("y")->GetText());
    velocity[2] = atof(physicsChild->FirstChildElement("velocity")->FirstChildElement("z")->GetText());

    acceleration[0] = atof(physicsChild->FirstChildElement("acceleration")->FirstChildElement("x")->GetText());
    acceleration[1] = atof(physicsChild->FirstChildElement("acceleration")->FirstChildElement("y")->GetText());
    acceleration[2] = atof(physicsChild->FirstChildElement("acceleration")->FirstChildElement("z")->GetText());

    viscosity = atof(physicsChild->FirstChildElement("viscosity")->GetText());

    maxVelocityDimLess = atof(physicsChild->FirstChildElement("max-velocity")->GetText());
    maxAccelerationDimLess = atof(physicsChild->FirstChildElement("max-acceleration")->GetText());
}

template <class T>
void CConfiguration<T>::interpretGridData(const tinyxml2::XMLNode* root)
{
    const tinyxml2::XMLNode* gridChild = root->FirstChildElement(TAG_NAME_CHILD_GRID);

    domainSize[0] = atoi(gridChild->FirstChildElement("domain-size")->FirstChildElement("x")->GetText());
    domainSize[1] = atoi(gridChild->FirstChildElement("domain-size")->FirstChildElement("y")->GetText());
    domainSize[2] = atoi(gridChild->FirstChildElement("domain-size")->FirstChildElement("z")->GetText());

    domainLength[0] = atof(gridChild->FirstChildElement("domian-length")->FirstChildElement("x")->GetText());
    domainLength[1] = atof(gridChild->FirstChildElement("domian-length")->FirstChildElement("y")->GetText());
    domainLength[2] = atof(gridChild->FirstChildElement("domian-length")->FirstChildElement("z")->GetText());

    numOfSubdomains[0] = atoi(gridChild->FirstChildElement("subdomain-num")->FirstChildElement("x")->GetText());
    numOfSubdomains[1] = atoi(gridChild->FirstChildElement("subdomain-num")->FirstChildElement("y")->GetText());
    numOfSubdomains[2] = atoi(gridChild->FirstChildElement("subdomain-num")->FirstChildElement("z")->GetText());

    CPUSubdomainRatio[0] = atof(gridChild->FirstChildElement("cpu-subdomain-ratio")->FirstChildElement("x")->GetText());
    CPUSubdomainRatio[1] = atof(gridChild->FirstChildElement("cpu-subdomain-ratio")->FirstChildElement("y")->GetText());
    CPUSubdomainRatio[2] = atof(gridChild->FirstChildElement("cpu-subdomain-ratio")->FirstChildElement("z")->GetText());
}

template <class T>
void CConfiguration<T>::interpretSimulationData(const tinyxml2::XMLNode* root)
{
	const tinyxml2::XMLNode* simulationChild = root->FirstChildElement(TAG_NAME_CHILD_SIMULATION);

    loops = atoi(simulationChild->FirstChildElement("loops")->GetText());
    timestep = atof(simulationChild->FirstChildElement("timestep")->GetText());

    doBenchmark = atoi(simulationChild->FirstChildElement("benchmark")->FirstChildElement("do")->GetText());
    doLogging = atoi(simulationChild->FirstChildElement("logging")->FirstChildElement("do")->GetText());
    doValidation = atoi(simulationChild->FirstChildElement("validation")->FirstChildElement("do")->GetText());
    doVisualization = atoi(simulationChild->FirstChildElement("visualization")->FirstChildElement("do")->GetText());

    const char* benchmarkOutputDirChar = simulationChild->FirstChildElement("benchmark")->FirstChildElement("output-dir")->GetText();
    const char* validationOutputDirChar = simulationChild->FirstChildElement("validation")->FirstChildElement("output-dir")->GetText();
    const char* visualizationOutputDirChar = simulationChild->FirstChildElement("visualization")->FirstChildElement("output-dir")->GetText();
    benchmarkOutputDir.assign(benchmarkOutputDirChar);
    validationOutputDir.assign(validationOutputDirChar);
    visualizationOutputDir.assign(visualizationOutputDirChar);

    visualizationRate = atoi(simulationChild->FirstChildElement("visualization")->FirstChildElement("rate")->GetText());
}

template <class T>
void CConfiguration<T>::interpretDeviceData(const tinyxml2::XMLNode* root)
{
	const tinyxml2::XMLNode* deviceChild = root->FirstChildElement(TAG_NAME_CHILD_DEVICE);

	dim3 configuration;
    configuration.x = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("init-grid-configuration")->FirstChildElement("x")->GetText());
    configuration.y = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("init-grid-configuration")->FirstChildElement("y")->GetText());
    configuration.z = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("init-grid-configuration")->FirstChildElement("z")->GetText());
    threadsPerBlock.push_back(configuration);

    configuration.x = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("alpha-grid-configuration")->FirstChildElement("x")->GetText());
    configuration.y = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("alpha-grid-configuration")->FirstChildElement("y")->GetText());
    configuration.z = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("alpha-grid-configuration")->FirstChildElement("z")->GetText());
    threadsPerBlock.push_back(configuration);

    configuration.x = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("beta-grid-configuration")->FirstChildElement("x")->GetText());
    configuration.y = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("beta-grid-configuration")->FirstChildElement("y")->GetText());
    configuration.z = atoi(deviceChild->FirstChildElement("grid-configuration")->FirstChildElement("beta-grid-configuration")->FirstChildElement("z")->GetText());
    threadsPerBlock.push_back(configuration);
}

template <class T>
void CConfiguration<T>::interpretXMLDoc()
{
    const tinyxml2::XMLNode* root = doc.FirstChildElement(TAG_NAME_ROOT);

    interpretPhysiscsData(root);
    interpretGridData(root);
    interpretSimulationData(root);
    interpretDeviceData(root);
}

template <class T>
void CConfiguration<T>::checkParameters()
{
    assert(timestep > (T)0);
    assert(maxVelocityDimLess > (T)0);
    assert(maxAccelerationDimLess > (T)0);
    assert(domainSize[0] > 0 && domainSize[1] > 0 && domainSize[2] > 0);
    assert(domainLength[0] > (T)0 && domainLength[1] > (T)0 && domainLength[2] > (T)0);
    assert(numOfSubdomains[0] > 0 && numOfSubdomains[1] > 0 && numOfSubdomains[2] > 0);
    assert(CPUSubdomainRatio[0] >= (T)0 && CPUSubdomainRatio[1] >= (T)0 && CPUSubdomainRatio[2] >= (T)0);
    assert(CPUSubdomainRatio[0] <= (T)1 && CPUSubdomainRatio[1] <= (T)1 && CPUSubdomainRatio[2] <= (T)1);
    assert(domainSize[0] % numOfSubdomains[0] == 0 && domainSize[1] % numOfSubdomains[1] == 0 && domainSize[2] % numOfSubdomains[2] == 0);
    assert(CMath<T>::abs(domainLength[0] / (T)domainSize[0] - domainLength[1] / (T)domainSize[1]) < std::numeric_limits<T>::epsilon());
    assert(CMath<T>::abs(domainLength[1] / (T)domainSize[1] - domainLength[2] / (T)domainSize[2]) < std::numeric_limits<T>::epsilon());
    assert(CMath<T>::abs(domainLength[2] / (T)domainSize[2] - domainLength[0] / (T)domainSize[0]) < std::numeric_limits<T>::epsilon());
    assert(visualizationRate > 0);
    assert(threadsPerBlock[0].x > 0 && threadsPerBlock[0].y > 0 && threadsPerBlock[0].z > 0);
    assert(threadsPerBlock[1].x > 0 && threadsPerBlock[1].y > 0 && threadsPerBlock[1].z > 0);
    assert(threadsPerBlock[2].x > 0 && threadsPerBlock[2].y > 0 && threadsPerBlock[2].z > 0);
    assert(getSize(threadsPerBlock[0]) <= 1024);
    assert(getSize(threadsPerBlock[1]) <= 1024);
    assert(getSize(threadsPerBlock[2]) <= 1024);
}

template <class T>
void CConfiguration<T>::print()
{
    std::cout << "----- CConfiguration<T>::print() -----" << std::endl;
    std::cout << "acceleration:                      " << acceleration << std::endl;
    std::cout << "velocity:                          " << velocity << std::endl;
    std::cout << "viscosity:                         " << viscosity << std::endl;
    std::cout << "max velocity (dimension less):     " << maxVelocityDimLess << std::endl;
    std::cout << "max acceleration (dimension less): " << maxAccelerationDimLess << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "domain size:                       " << domainSize << std::endl;
    std::cout << "domain length:                     " << domainLength << std::endl;
    std::cout << "number of subdomains:              " << numOfSubdomains << std::endl;
    std::cout << "cpu/subdomain ratio:               " << CPUSubdomainRatio << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "loops:                             " << loops << std::endl;
    std::cout << "timestep:                          " << timestep << std::endl;
    std::cout << "do benchmark:                      " << doBenchmark << std::endl;
    std::cout << "do logging:                        " << doLogging << std::endl;
    std::cout << "do validation:                     " << doValidation << std::endl;
    std::cout << "do visualization:                  " << doVisualization << std::endl;
    std::cout << "benchmark directory:               " << benchmarkOutputDir << std::endl;
    std::cout << "validation directory:              " << validationOutputDir << std::endl;
    std::cout << "visualization directory:           " << visualizationOutputDir << std::endl;
    std::cout << "visualization rate:                " << visualizationRate << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "init grid configuration:           [" << threadsPerBlock[0].x << ", " << threadsPerBlock[0].y << ", " << threadsPerBlock[0].z << "]" << std::endl;
    std::cout << "alpha grid configuration:          [" << threadsPerBlock[1].x << ", " << threadsPerBlock[1].y << ", " << threadsPerBlock[1].z << "]" << std::endl;
    std::cout << "beta grid configuration:           [" << threadsPerBlock[2].x << ", " << threadsPerBlock[2].y << ", " << threadsPerBlock[2].z << "]" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
}

template <class T>
void CConfiguration<T>::serialize()
{
	numOfSubdomains[0] = numOfSubdomains[1] = numOfSubdomains[2] = 1;
}

template class CConfiguration<double>;
template class CConfiguration<float>;
