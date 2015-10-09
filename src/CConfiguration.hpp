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

#include <cstdlib>
#include <iostream>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include <tinyxml2.h>

#include "libmath/CVector.hpp"

/*
 * Class CConfiguration stores the necessary information for the simulation process.
 */
namespace txml = tinyxml2;

#define TAG_NAME_ROOT           "lbm-configuration"
#define TAG_NAME_CHILD_ONE      "physics"
#define TAG_NAME_CHILD_TWO      "grid"
#define TAG_NAME_CHILD_THREE    "simulation"
#define TAG_NAME_CHILD_FOUR     "device"

template <typename T>
class CConfiguration
{
private:
    txml::XMLDocument doc;

    int load_xml(std::string file_name);
    void interpret_physiscs_data(const txml::XMLNode* root);
    void interpret_grid_data(const txml::XMLNode* root);
    void interpret_simulation_data(const txml::XMLNode* root);
    void interpret_device_data(const txml::XMLNode* root);
    void interpret_xml_doc();

public:
    // grid data
    CVector<3, int> domain_size;
    CVector<3, int> subdomain_num;
    CVector<3, T> domain_length;

    // physics configuration data
    CVector<3, T> gravitation;       ///< Specify the gravitation vector
    T viscosity;

    CVector<4, T> drivenCavityVelocity;
    // device configuration data
    size_t computation_kernel_count;
    size_t threads_per_dimension;
    int device_nr;

    // simulation configuration data
    bool do_visualization;
    T timestep;
    int loops;
    bool do_validate;

    // TODO: lbm_opencl_number_of_registers_list, lbm_opencl_number_of_threads_list
    std::vector<dim3> lbm_opencl_number_of_registers_list;
    std::vector<dim3> lbm_opencl_number_of_threads_list;

    // domain configuration data

    bool debug_mode;

    CConfiguration();
    CConfiguration(std::string file_name);
    ~CConfiguration();

    void loadFile(std::string file_name);
    void printMe();
};

#endif
