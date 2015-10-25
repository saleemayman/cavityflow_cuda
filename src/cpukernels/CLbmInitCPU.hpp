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

#ifndef LBM_INITCPU_CUH
#define LBM_INITCPU_CUH

#include <vector>

#include "../libmath/CVector.hpp"
#include "../CDomain.hpp"
#include "../common.h"

template<typename T>
class CLbmInitCPU
{
private:
    int domainCellsCPU;
    CVector<3, int> domainSize;
    CVector<3, int> domainSizeGPU; 
    CVector<3, int> innerGPUCoordinates;
    CVector<3, int> outerGPUCoordinates;

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

    Flag setFlags(int xPosition, int yPosition, int zPosition);

public:
	CLbmInitCPU(
            CVector<3, int> domainSize,
            CVector<3, int> domainSizeGPU,
            std::vector<Flag>& boundaryConditions);
	~CLbmInitCPU();

	void initLbm(
		T *global_dd,
        Flag *flags,
        T *velocityArray,
        T *density);
};
#endif
