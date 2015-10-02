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

#ifndef LBM_INIT_CUH
#define LBM_INIT_CUH

#include "../common.h"

template<typename T>
__global__ void lbm_init(
		T *global_dd,           // density distributions
        Flag *flags,            // flags
        T *velocity_array,      // velocity array (first all x components, then all y components, then z...)
        T *density,             // densities
        Flag boundaryConditionRight,
        Flag boundaryConditionLeft,
        Flag boundaryConditionTop,
        Flag boundaryConditionBottom,
        Flag boundaryConditionFront,
        Flag boundaryConditionBack,
        T drivenCavityVelocity, // velocity parameters for modification of density distributions
        const int domainCells_x,
        const int domainCells_y,
        const int domainCells_z);

#endif
