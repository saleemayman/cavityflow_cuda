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

#include "CLbmInitCPU.hpp"

template<class T>
CLbmInitCPU<T>::CLbmInitCPU(
                    CVector<3, int> domainSize,
                    CVector<3, int> domainSizeGPU,
                    std::vector<Flag>& boundaryConditions) : 
                        domainSize(domainSize),
                        domainSizeGPU(domainSizeGPU),
                        boundaryConditionRight(boundaryConditions[0]),
                        boundaryConditionLeft(boundaryConditions[1]),
                        boundaryConditionTop(boundaryConditions[2]),
                        boundaryConditionBottom(boundaryConditions[3]),
                        boundaryConditionFront(boundaryConditions[4]),
                        boundaryConditionBack(boundaryConditions[5])
{
    domainCellsCPU = domainSize.elements() - domainSizeGPU.elements();  

    // limits of the inner GPU domain w.r.t the CPU domain
    innerCPULimit[0] = (domainSize[0] - domainSizeGPU[0]) / 2;
    outerCPULimit[0] = innerCPULimit[0] + domainSizeGPU[0] + 1;
    innerCPULimit[1] = (domainSize[1] - domainSizeGPU[1]) / 2;
    outerCPULimit[1] = innerCPULimit[1] + domainSizeGPU[1] + 1;
    innerCPULimit[2] = (domainSize[2] - domainSizeGPU[2]) / 2;
    outerCPULimit[2] = innerCPULimit[2] + domainSizeGPU[2] + 1;
}

template<class T>
CLbmInitCPU<T>::~CLbmInitCPU()
{
}


template<class T>
void CLbmInitCPU<T>::setFlags(std::vector<Flag> &flags)
{
    // linear cell index
    int cellID = 0;
    Flag flag;

    /*
     * Iterate over all CPU domain cells in the following order:
     * x-cells, y-cells, z-cells.
     */
    for (int z = 0; z < domainSize[2]; z++)
    {
        for (int y = 0; y < domainSize[1]; y++)
        {
            for (int x = 0; x < domainSize[0]; x++)
            {
                // initialize flag field
                flag = FLUID;

                // if inside GPU domain then skip cell
                if ((x > innerCPULimit[0]-1 && x < outerCPULimit[0]-1) &&
                    (y > innerCPULimit[1]-1 && y < outerCPULimit[1]-1) &&
                    (z > innerCPULimit[2]-1 && z < outerCPULimit[2]-1))
                {
                    continue;
                }

                // set the flags for x-direction                
                if (x == 0)
                    flag = boundaryConditionRight;
                if (x == innerCPULimit[0]-1 &&
                        (y > innerCPULimit[1]-1 && y < outerCPULimit[1]-1) &&
                        (z > innerCPULimit[2]-1 && z < outerCPULimit[2]-1) &&
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (x == outerCPULimit[0]-1 &&
                        (y > innerCPULimit[1]-1 && y < outerCPULimit[1]-1) && 
                        (z > innerCPULimit[2]-1 && z < outerCPULimit[2]-1) &&
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (x == domainSize[0]-1 && flag != OBSTACLE)
                    flag = boundaryConditionLeft;
                
                // set the flags for y-direction                
                if (y == 0 && flag != OBSTACLE)
                    flag = boundaryConditionTop;
                if (y == innerCPULimit[1]-1 && 
                        (x > innerCPULimit[0]-1 && x < outerCPULimit[0]-1) && 
                        (z > innerCPULimit[2]-1 && z < outerCPULimit[2]-1) && 
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (y == outerCPULimit[1]-1 && 
                        (x > innerCPULimit[0]-1 && x < outerCPULimit[0]-1) && 
                        (z > innerCPULimit[2]-1 && z < outerCPULimit[2]-1) && 
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (y == domainSize[1]-1 && flag != OBSTACLE)
                    flag = boundaryConditionBottom;
                
                // set the flags for z-direction                
                if (z == 0 && flag != OBSTACLE)
                    flag = boundaryConditionFront;
                if (z == innerCPULimit[2]-1 && 
                        (x > innerCPULimit[0]-1 && x < outerCPULimit[0]-1) &&
                        (y > innerCPULimit[1]-1 && y < outerCPULimit[1]-1) &&  
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (z == outerCPULimit[2]-1 && 
                        (x > innerCPULimit[0]-1 && x < outerCPULimit[0]-1) &&
                        (y > innerCPULimit[1]-1 && y < outerCPULimit[1]-1) && 
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (z == domainSize[2]-1 && flag != OBSTACLE)
                    flag = boundaryConditionBack;

                flags[cellID] = flag;

                // increment cell linear id
                cellID++;
            }
        }
    }
}

template<class T>
void CLbmInitCPU<T>::initLbm(
            std::vector<T> &densityDistributions,
            std::vector<Flag> &flags,
            std::vector<T> &velocities,
            std::vector<T> &densities,
            const bool storeDensities,
            const bool storeVelocities)
{
    // initialize the flag field
    setFlags(flags);

    /*
     * Iterate over all CPU domain cells.
     * 
     */
    for (int i = 0; i < domainCellsCPU; i++)
    {
        // initialize the cell density distributions
        velocity_x = 0;
        velocity_y = 0;
        velocity_z = 0;
        rho = 1.0f;

        // compute and store velocity
        vela2 = velocity_x*velocity_x;
        dd_param = rho - (T)(3.0f/2.0f)*(vela2);

        densityDistributions[i + 0*domainCellsCPU] = eq_dd_a0(velocity_x, vela2, dd_param);
        densityDistributions[i + 1*domainCellsCPU] = eq_dd_a1(velocity_x, vela2, dd_param);

        vela2 = velocity_y*velocity_y;
        densityDistributions[i + 2*domainCellsCPU] = eq_dd_a0(velocity_y, vela2, dd_param);
        densityDistributions[i + 3*domainCellsCPU] = eq_dd_a1(velocity_y, vela2, dd_param);

        /***********************
         * DD1
         ***********************/
        vela_velb = velocity_x+velocity_y;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 4*domainCellsCPU] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 5*domainCellsCPU] = eq_dd5(vela_velb, vela_velb_2, dd_param);


        vela_velb = velocity_x-velocity_y;
        vela_velb_2 = vela_velb*vela_velb;
        densityDistributions[i + 6*domainCellsCPU] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 7*domainCellsCPU] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD2
         ***********************/
        vela_velb = velocity_x+velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 8*domainCellsCPU] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 9*domainCellsCPU] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        vela_velb = velocity_x-velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 10*domainCellsCPU] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 11*domainCellsCPU] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD3
         ***********************/
        vela_velb = velocity_y+velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 12*domainCellsCPU] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 13*domainCellsCPU] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        vela_velb = velocity_y-velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 14*domainCellsCPU] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 15*domainCellsCPU] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD4
         ***********************/
        vela2 = velocity_z*velocity_z;

        densityDistributions[i + 16*domainCellsCPU] = eq_dd_a0(velocity_z, vela2, dd_param);
        densityDistributions[i + 17*domainCellsCPU] = eq_dd_a1(velocity_z, vela2, dd_param);
        densityDistributions[i + 18*domainCellsCPU] = eq_dd18(dd_param);

        if (storeVelocities)
        {
            velocities[i] = velocity_x;
            velocities[i + domainCellsCPU] = velocity_y;
            velocities[i + 2*domainCellsCPU] = velocity_z;
        }

        if (storeDensities)
        {
            densities[i] = rho;
        }
    }
}

template class CLbmInitCPU<double>;
template class CLbmInitCPU<float>;
