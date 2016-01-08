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
        std::vector<Flag> boundaryConditions) :
        domainSize(domainSize),
        boundaryConditions(boundaryConditions),
        numOfCells(domainSize.elements())
{
}

template<class T>
CLbmInitCPU<T>::~CLbmInitCPU()
{
}

template<class T>
void CLbmInitCPU<T>::setFlags(Flag* flags)
{
    int id = 0;
    Flag flag;

    for (int x = 0; x < domainSize[0]; x++)
    {
        for (int y = 0; y < domainSize[1]; y++)
        {
        	for (int z = 0; z < domainSize[2]; z++)
            {
                // initialize flag field
                flag = FLUID;
                // set the flags for x-direction                
                if (x == 0)
                    flag = boundaryConditions[0];
                if (x == domainSize[0]-1 && flag != OBSTACLE)
                    flag = boundaryConditions[1];
                
                // set the flags for y-direction                
                if (y == 0 && flag != OBSTACLE)
                    flag = boundaryConditions[2];
                if (y == domainSize[1]-1 && flag != OBSTACLE)
                    flag = boundaryConditions[3];
                
                // set the flags for z-direction                
                if (z == 0 && flag != OBSTACLE)
                    flag = boundaryConditions[4];
                if (z == domainSize[2]-1 && flag != OBSTACLE)
                    flag = boundaryConditions[5];

                flags[id] = flag;

                id++;
            }
        }
    }
}

template<class T>
void CLbmInitCPU<T>::initLbm(
        T* densityDistributions,
        Flag* flags,
        T* densities,
        T* velocities,
        const bool storeDensities,
        const bool storeVelocities)
{
    setFlags(flags);

    for (int id = 0; id < numOfCells; id++)
    {
        // initialize the cell density distributions
        velocity_x = (T)0;
        velocity_y = (T)0;
        velocity_z = (T)0;
        rho = (T)1;

        // compute and store velocity
        vela2 = velocity_x*velocity_x;
        dd_param = rho - ((T)3/(T)2)*(vela2);

        densityDistributions[id + 0*numOfCells] = eq_dd_a0(velocity_x, vela2, dd_param);
        densityDistributions[id + 1*numOfCells] = eq_dd_a1(velocity_x, vela2, dd_param);

        vela2 = velocity_y*velocity_y;
        densityDistributions[id + 2*numOfCells] = eq_dd_a0(velocity_y, vela2, dd_param);
        densityDistributions[id + 3*numOfCells] = eq_dd_a1(velocity_y, vela2, dd_param);

        /***********************
         * DD1
         ***********************/
        vela_velb = velocity_x+velocity_y;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[id + 4*numOfCells] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[id + 5*numOfCells] = eq_dd5(vela_velb, vela_velb_2, dd_param);


        vela_velb = velocity_x-velocity_y;
        vela_velb_2 = vela_velb*vela_velb;
        densityDistributions[id + 6*numOfCells] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[id + 7*numOfCells] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD2
         ***********************/
        vela_velb = velocity_x+velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[id + 8*numOfCells] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[id + 9*numOfCells] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        vela_velb = velocity_x-velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[id + 10*numOfCells] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[id + 11*numOfCells] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD3
         ***********************/
        vela_velb = velocity_y+velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[id + 12*numOfCells] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[id + 13*numOfCells] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        vela_velb = velocity_y-velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[id + 14*numOfCells] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[id + 15*numOfCells] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD4
         ***********************/
        vela2 = velocity_z*velocity_z;

        densityDistributions[id + 16*numOfCells] = eq_dd_a0(velocity_z, vela2, dd_param);
        densityDistributions[id + 17*numOfCells] = eq_dd_a1(velocity_z, vela2, dd_param);
        densityDistributions[id + 18*numOfCells] = eq_dd18(dd_param);

        if (storeVelocities)
        {
            velocities[id] = velocity_x;
            velocities[id + numOfCells] = velocity_y;
            velocities[id + 2*numOfCells] = velocity_z;
        }

        if (storeDensities)
            densities[id] = rho;
    }
}

template class CLbmInitCPU<double>;
template class CLbmInitCPU<float>;
