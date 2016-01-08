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

#ifndef LBM_INITCPU_HPP
#define LBM_INITCPU_HPP

#include <vector>

#include "../libmath/CVector.hpp"
#include "../common.h"

template<typename T>
class CLbmInitCPU
{
private:
    CVector<3, int> domainSize;
    std::vector<Flag> boundaryConditions;
    int numOfCells;

    T dd_param;
    T rho;
    T vela2, vela_velb, vela_velb_2;
    T velocity_x, velocity_y, velocity_z;

    void setFlags(Flag* flags);

public:
    CLbmInitCPU(
            CVector<3, int> domainSize,
            std::vector<Flag> boundaryConditions);
    ~CLbmInitCPU();

    void initLbm(
        T* densityDistributions,
        Flag* flags,
        T* densities,
        T* velocities,
        const bool storeDensities,
        const bool storeVelocities);
};
#endif
