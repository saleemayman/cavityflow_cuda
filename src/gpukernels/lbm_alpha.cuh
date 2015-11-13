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

#ifndef LBM_ALPHA_CUH
#define LBM_ALPHA_CUH

#include "../common.h"

#define GRAVITATION 0

template<typename T>
__global__ void lbm_kernel_alpha(
        T* global_dd,                 // density distributions
        const Flag* flag_array,       // flags
        T* velocity,                  // velocities
        T* density,                   // densities
        const T inv_tau,
        const T gravitation_x,
        const T gravitation_y,
        const T gravitation_z,
        const T drivenCavityVelocity, // velocity parameters for modification of density distributions
        const int domainCellsX,
        const int domainCellsY,
        const int domainCellsZ,
        const bool storeDensities,
        const bool storeVelocities);

#endif
