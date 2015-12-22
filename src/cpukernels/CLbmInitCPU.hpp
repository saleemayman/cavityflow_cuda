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
#include <algorithm>

#include "../libmath/CVector.hpp"
#include "../CDomain.hpp"
#include "../common.h"

template<typename T>
class CLbmInitCPU
{
private:
    int domainCellsCPUWithHalo;
    CVector<3, int> domainSizeWithHalo;
#if !TOP_DOWN_DECOMPOSITION
    CVector<3, int> hollowCPULeftLimit;
    CVector<3, int> hollowCPURightLimit;
    std::vector<int> *localToGlobalIndexMap;
#endif

    // declare simulation related variables
    Flag flag;
    Flag boundaryConditionRight;
    Flag boundaryConditionLeft;
    Flag boundaryConditionTop;
    Flag boundaryConditionBottom;
    Flag boundaryConditionFront;
    Flag boundaryConditionBack;

    T velocity_x, velocity_y, velocity_z;
    T dd_param;
    T rho;
    T vela2;
    T vela_velb;
    T vela_velb_2;

    void setFlags(std::vector<Flag> &flags);

public:
#if !TOP_DOWN_DECOMPOSITION
    CLbmInitCPU(
            int domainCellsCPUWithHalo,
            CVector<3, int> domainSizeWithHalo,
            CVector<3, int> hollowCPULeftLimit,
            CVector<3, int> hollowCPURightLimit,
            std::vector<Flag>& boundaryConditions);
#else
    CLbmInitCPU(
            int domainCellsCPUWithHalo,
            CVector<3, int> domainSizeWithHalo,
            std::vector<Flag>& boundaryConditions);
#endif
    ~CLbmInitCPU();

    void initLbm(
        std::vector<T> &densityDistributions,
        std::vector<Flag> &flags,
        std::vector<T> &velocities,
        std::vector<T> &densities,
        const bool storeDensities,
        const bool storeVelocities);

#if !TOP_DOWN_DECOMPOSITION
    int getLocalIndex(int globalId);
    int getGlobalIndex(int localId);

    std::vector<int>* getCellIndexMap();
#endif
};
#endif
