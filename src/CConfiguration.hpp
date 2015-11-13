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

#ifndef CCONFIGURATION_HPP
#define CCONFIGURATION_HPP

#include <iostream>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "../external/tinyxml2/tinyxml2.h"

#include "libmath/CVector.hpp"

/*
 * Class CConfiguration stores the necessary information for the simulation process.
 */
#define TAG_NAME_ROOT             "lbm-configuration"
#define TAG_NAME_CHILD_PHYSICS    "physics"
#define TAG_NAME_CHILD_GRID       "grid"
#define TAG_NAME_CHILD_SIMULATION "simulation"
#define TAG_NAME_CHILD_DEVICE     "device"

template <typename T>
class CConfiguration
{
private:
    tinyxml2::XMLDocument doc;

    void interpretPhysiscsData(const tinyxml2::XMLNode* root);
    void interpretGridData(const tinyxml2::XMLNode* root);
    void interpretSimulationData(const tinyxml2::XMLNode* root);
    void interpretDeviceData(const tinyxml2::XMLNode* root);
    void interpretXMLDoc();
    void checkParameters();

public:
    /*
     * physics configuration data
     */
    CVector<3, T> velocity;
    CVector<3, T> acceleration;
    T viscosity;
    T maxVelocityDimLess;
    T maxAccelerationDimLess;

    /*
     * grid configuration data
     */
    CVector<3, int> domainSize;
    CVector<3, T> domainLength;
    CVector<3, int> numOfSubdomains;
    CVector<3, T> CPUSubdomainRatio;

    /*
     * simulation configuration data
     */
    int loops;
    T timestep;
    bool doBenchmark;
    bool doLogging;
    bool doValidation;
    bool doVisualization;
    std::string benchmarkOutputDir;
    std::string validationOutputDir;
    std::string visualizationOutputDir;
    int visualizationRate;

    /*
     * Device configuration data
     */
    std::vector<dim3> threadsPerBlock;

    CConfiguration(std::string fileName);
    ~CConfiguration();

    void print();
    void serialize();
};

#endif
