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

#ifndef LBM_BETACPU_HPP
#define LBM_BETACPU_HPP

#include "../common.h"
#include "../libmath/CVector.hpp"

template<class T>
class CLbmBetaCPU
{
private:
    CVector<3, int> domainSize;
    CVector<3, T> gravitation;
    int numOfCells;

    int deltaPosX, deltaNegX;
    int deltaPosY, deltaNegY;
    int deltaPosZ, deltaNegZ;

    T dd_param;
    T dd0, dd1, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9, dd10, dd11, dd12, dd13, dd14, dd15, dd16, dd17, dd18;
    T rho;
    T vel2, vela2, vela_velb, vela_velb_2;
    T velocity_x, velocity_y, velocity_z;

    int domainWrap(int A, int domainCells, bool isPowTwo);
public:
    CLbmBetaCPU(
            CVector<3, int> domainSize,
            CVector<3, T> gravitation);
    ~CLbmBetaCPU();

    void betaKernelCPU(
            T* densityDistributions,
            Flag* flags,
            T* densities,
            T* velocities,
            const T inv_tau,
            const T drivenCavityVelocity,
            CVector<3, int> origin,
            CVector<3, int> size,
            const bool isDomainPowOfTwo,
            const bool storeDensities,
            const bool storeVelocities);
};
#endif
