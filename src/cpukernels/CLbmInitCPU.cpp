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
                    int domainCellsCPUWithHalo,
                    CVector<3, int> domainSizeWithHalo,
                    CVector<3, int> hollowCPULeftLimit,
                    CVector<3, int> hollowCPURightLimit,
                    std::vector<Flag>& boundaryConditions) : 
                        domainCellsCPUWithHalo(domainCellsCPUWithHalo),
                        domainSizeWithHalo(domainSizeWithHalo),
                        hollowCPULeftLimit(hollowCPULeftLimit),
                        hollowCPURightLimit(hollowCPURightLimit),
                        boundaryConditionRight(boundaryConditions[0]),
                        boundaryConditionLeft(boundaryConditions[1]),
                        boundaryConditionTop(boundaryConditions[2]),
                        boundaryConditionBottom(boundaryConditions[3]),
                        boundaryConditionFront(boundaryConditions[4]),
                        boundaryConditionBack(boundaryConditions[5])
{
    localToGlobalIndexMap = new std::vector<int>(domainCellsCPUWithHalo);
}

template<class T>
CLbmInitCPU<T>::~CLbmInitCPU()
{
    delete localToGlobalIndexMap;
}


template<class T>
void CLbmInitCPU<T>::setFlags(std::vector<Flag> &flags)
{
    // linear cell index
    int cellID = 0;
    int globalID = 0;
    Flag flag;

    /*
     * Iterate over all CPU domain cells in the following order:
     * x-cells, y-cells, z-cells.
     */
    for (int z = 0; z < domainSizeWithHalo[2]; z++)
    {
        for (int y = 0; y < domainSizeWithHalo[1]; y++)
        {
            for (int x = 0; x < domainSizeWithHalo[0]; x++)
            {
                // initialize flag field
                flag = FLUID;

                // if inside GPU domain then skip cell
                if ((x > (hollowCPULeftLimit[0] + 1) && x < (hollowCPURightLimit[0] - 1)) &&
                    (y > (hollowCPULeftLimit[1] + 1) && y < (hollowCPURightLimit[1] - 1)) &&
                    (z > (hollowCPULeftLimit[2] + 1) && z < (hollowCPURightLimit[2] - 1)))
                {
                    globalID++;
                    continue;
                }

                // set the flags for x-direction                
                if (x == 0)
                    flag = boundaryConditionRight;
                if (x == hollowCPULeftLimit[0] + 1 &&
                        (y > hollowCPULeftLimit[1] && y < hollowCPURightLimit[1]) &&
                        (z > hollowCPULeftLimit[2] && z < hollowCPURightLimit[2]) &&
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (x == hollowCPURightLimit[0] - 1 &&
                        (y > hollowCPULeftLimit[1] && y < hollowCPURightLimit[1]) && 
                        (z > hollowCPULeftLimit[2] && z < hollowCPURightLimit[2]) &&
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (x == domainSizeWithHalo[0]-1 && flag != OBSTACLE)
                    flag = boundaryConditionLeft;
                
                // set the flags for y-direction                
                if (y == 0 && flag != OBSTACLE)
                    flag = boundaryConditionTop;
                if (y == hollowCPULeftLimit[1] + 1 && 
                        (x > hollowCPULeftLimit[0] && x < hollowCPURightLimit[0]) && 
                        (z > hollowCPULeftLimit[2] && z < hollowCPURightLimit[2]) && 
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (y == hollowCPURightLimit[1] - 1 && 
                        (x > hollowCPULeftLimit[0] && x < hollowCPURightLimit[0]) && 
                        (z > hollowCPULeftLimit[2] && z < hollowCPURightLimit[2]) && 
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (y == domainSizeWithHalo[1]-1 && flag != OBSTACLE)
                    flag = boundaryConditionBottom;
                
                // set the flags for z-direction                
                if (z == 0 && flag != OBSTACLE)
                    flag = boundaryConditionFront;
                if (z == hollowCPULeftLimit[2] + 1 && 
                        (x > hollowCPULeftLimit[0] && x < hollowCPURightLimit[0]) &&
                        (y > hollowCPULeftLimit[1] && y < hollowCPURightLimit[1]) &&  
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (z == hollowCPURightLimit[2] - 1 && 
                        (x > hollowCPULeftLimit[0] && x < hollowCPURightLimit[0]) &&
                        (y > hollowCPULeftLimit[1] && y < hollowCPURightLimit[1]) && 
                        flag != OBSTACLE)
                    flag = GHOST_LAYER;
                if (z == domainSizeWithHalo[2]-1 && flag != OBSTACLE)
                    flag = boundaryConditionBack;

                flags[cellID] = flag;

                // add the global cell index to an array
                localToGlobalIndexMap->operator[](cellID) = globalID;

                // increment cell linear id
                cellID++;
                globalID++;
            }
        }
    }
}

template<class T>
int CLbmInitCPU<T>::getLocalIndex(int globalId)
{
    std::vector<int>::iterator localIndex;
    localIndex = std::lower_bound(localToGlobalIndexMap->begin(), localToGlobalIndexMap->end(), globalId);

    if (localIndex != localToGlobalIndexMap->end() && *localIndex == globalId)
    {
        return (int)(localIndex - localToGlobalIndexMap->begin());
    }
    else
    {
        return -1;
    }
}

template<class T>
int CLbmInitCPU<T>::getGlobalIndex(int localId)
{
    return localToGlobalIndexMap->operator[](localId);
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
    for (int i = 0; i < domainCellsCPUWithHalo; i++)
    {
        // initialize the cell density distributions
        velocity_x = 0;
        velocity_y = 0;
        velocity_z = 0;
        rho = 1.0f;

        // compute and store velocity
        vela2 = velocity_x*velocity_x;
        dd_param = rho - (T)(3.0f/2.0f)*(vela2);

        densityDistributions[i + 0*domainCellsCPUWithHalo] = eq_dd_a0(velocity_x, vela2, dd_param);
        densityDistributions[i + 1*domainCellsCPUWithHalo] = eq_dd_a1(velocity_x, vela2, dd_param);

        vela2 = velocity_y*velocity_y;
        densityDistributions[i + 2*domainCellsCPUWithHalo] = eq_dd_a0(velocity_y, vela2, dd_param);
        densityDistributions[i + 3*domainCellsCPUWithHalo] = eq_dd_a1(velocity_y, vela2, dd_param);

        /***********************
         * DD1
         ***********************/
        vela_velb = velocity_x+velocity_y;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 4*domainCellsCPUWithHalo] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 5*domainCellsCPUWithHalo] = eq_dd5(vela_velb, vela_velb_2, dd_param);


        vela_velb = velocity_x-velocity_y;
        vela_velb_2 = vela_velb*vela_velb;
        densityDistributions[i + 6*domainCellsCPUWithHalo] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 7*domainCellsCPUWithHalo] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD2
         ***********************/
        vela_velb = velocity_x+velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 8*domainCellsCPUWithHalo] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 9*domainCellsCPUWithHalo] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        vela_velb = velocity_x-velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 10*domainCellsCPUWithHalo] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 11*domainCellsCPUWithHalo] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD3
         ***********************/
        vela_velb = velocity_y+velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 12*domainCellsCPUWithHalo] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 13*domainCellsCPUWithHalo] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        vela_velb = velocity_y-velocity_z;
        vela_velb_2 = vela_velb*vela_velb;

        densityDistributions[i + 14*domainCellsCPUWithHalo] = eq_dd4(vela_velb, vela_velb_2, dd_param);
        densityDistributions[i + 15*domainCellsCPUWithHalo] = eq_dd5(vela_velb, vela_velb_2, dd_param);

        /***********************
         * DD4
         ***********************/
        vela2 = velocity_z*velocity_z;

        densityDistributions[i + 16*domainCellsCPUWithHalo] = eq_dd_a0(velocity_z, vela2, dd_param);
        densityDistributions[i + 17*domainCellsCPUWithHalo] = eq_dd_a1(velocity_z, vela2, dd_param);
        densityDistributions[i + 18*domainCellsCPUWithHalo] = eq_dd18(dd_param);

        if (storeVelocities)
        {
            velocities[i] = velocity_x;
            velocities[i + domainCellsCPUWithHalo] = velocity_y;
            velocities[i + 2*domainCellsCPUWithHalo] = velocity_z;
        }

        if (storeDensities)
        {
            densities[i] = rho;
        }
    }
}


template<class T>
std::vector<int>* CLbmInitCPU<T>::getCellIndexMap()
{
    return localToGlobalIndexMap;
}


template class CLbmInitCPU<double>;
template class CLbmInitCPU<float>;
