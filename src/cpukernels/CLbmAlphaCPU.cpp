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

#include "CLbmAlphaCPU.hpp"

template<class T>
CLbmAlphaCPU<T>::CLbmAlphaCPU(
                    int domainCellsCPUWithHalo,
                    CVector<3, int> domainSizeWithHalo,
                    CVector<3, int> hollowCPULeftLimit,
                    CVector<3, int> hollowCPURightLimit,
                    //std::vector<int> *localToGlobalIndexMap,
                    CLbmInitCPU<T> *initLbmCPU,
                    CVector<3, T> gravitation) :
                        domainCellsCPUWithHalo(domainCellsCPUWithHalo),
                        domainSizeWithHalo(domainSizeWithHalo),
                        hollowCPULeftLimit(hollowCPULeftLimit),
                        hollowCPURightLimit(hollowCPURightLimit),
                        //localToGlobalIndexMap(localToGlobalIndexMap),
                        initLbmCPU(initLbmCPU),
                        gravitation(gravitation)
{
}

template<class T>
CLbmAlphaCPU<T>::~CLbmAlphaCPU()
{
	if (isSubRegion)
		delete[] localIndices;
}


template<class T>
void CLbmAlphaCPU<T>::alphaKernelCPU(
                std::vector<T> &densityDistributions,
                std::vector<Flag> &flags,
                std::vector<T> &velocities,
                std::vector<T> &densities,
                const T inv_tau,
                const T drivenCavityVelocity,
                CVector<3, int> origin,
                CVector<3, int> size,
                const bool storeDensities,
                const bool storeVelocities)
{
    int i;  // cell local linear id
	int startIndex, endIndex;

    /*
     * Check if computation has to be done for a sub-region, if yes then create a local vector
     * containing the local linear indices of all the linear global cells for which we need to 
     * do the computation (the global cells are specified in size parameter).
     */
    if (size[0] == domainSizeWithHalo[0] && size[1] == domainSizeWithHalo[1] && size[2] == domainSizeWithHalo[2])
    {
        startIndex = 0; 
        endIndex = domainCellsCPUWithHalo;
        isSubRegion = (bool)0;
    }
    else
        isSubRegion = (bool)1;

    if (isSubRegion)
    {
        startIndex = 0; 
        endIndex = size[0] * size[1] * size[2]; 
        localIndices = new std::vector<int>(endIndex);

        for (int i = 0; i < size[2]; i++)
        {
            for (int j = 0; j < size[1]; j++)
            {
                for (int k = 0; k < size[0]; k++)
                {
                    localIndices->operator[](k + j * size[0] + i * size[0] * size[1]) = initLbmCPU->getLocalIndex((origin[0] +  k) + (origin[1] + j) * domainSizeWithHalo[0] + (origin[2] + i) * (domainSizeWithHalo[1] * domainSizeWithHalo[2]));
                }
            }
        }
    }

    /*
     * Iterate over all CPU domain cells.
     */
    for (int id = startIndex; id < endIndex; id++)
    {
        if (isSubRegion)
        {
            i = localIndices->operator[](id);
        }
        else
        {
            i = id;
        }
            
        /*
         * skip cell if it is a ghost cell. Note: a ghost cell also means
         * that the cell lies at the boundary of the GPU domain and we 
         * have to skip them as well.
         */
        if  (flags[i] == GHOST_LAYER)
            continue;

        /*
         * we have to sum the densities up in a specific order.
         * otherwise it seems that we run into numerical errors.
         */
    
        /* +++++++++++
         * +++ DD0 +++
         * +++++++++++
         * 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
         */
        dd0 = densityDistributions[i + 0*domainCellsCPUWithHalo];
        rho = dd0;
        velocity_x = dd0;
       
        dd1 = densityDistributions[i + 1*domainCellsCPUWithHalo]; 
        rho += dd1;
        velocity_x -= dd1;

        dd2 = densityDistributions[i + 2*domainCellsCPUWithHalo]; 
        rho += dd2;
        velocity_y = dd2;
    
        dd3 = densityDistributions[i + 3*domainCellsCPUWithHalo]; 
        rho += dd3;
        velocity_y -= dd3;
    
        /* +++++++++++
         * +++ DD1 +++
         * +++++++++++
         * 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
         */
        dd4 = densityDistributions[i + 4*domainCellsCPUWithHalo]; 
        rho += dd4;
        velocity_x += dd4;
        velocity_y += dd4;

        dd5 = densityDistributions[i + 5*domainCellsCPUWithHalo]; 
        rho += dd5;
        velocity_x -= dd5;
        velocity_y -= dd5;
    
        dd6 = densityDistributions[i + 6*domainCellsCPUWithHalo]; 
        rho += dd6;
        velocity_x += dd6;
        velocity_y -= dd6;

        dd7 = densityDistributions[i + 7*domainCellsCPUWithHalo]; 
        rho += dd7;
        velocity_x -= dd7;
        velocity_y += dd7;
    
        /* +++++++++++
         * +++ DD2 +++
         * +++++++++++
         * 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
         */
        dd8 = densityDistributions[i + 8*domainCellsCPUWithHalo]; 
        rho += dd8;
        velocity_x += dd8;
        velocity_z = dd8;

        dd9 = densityDistributions[i + 9*domainCellsCPUWithHalo]; 
        rho += dd9;
        velocity_x -= dd9;
        velocity_z -= dd9;
    
        dd10 = densityDistributions[i + 10*domainCellsCPUWithHalo]; 
        rho += dd10;
        velocity_x += dd10;
        velocity_z -= dd10;

        dd11 = densityDistributions[i + 11*domainCellsCPUWithHalo]; 
        rho += dd11;
        velocity_x -= dd11;
        velocity_z += dd11;
    
        /* +++++++++++
         * +++ DD3 +++
         * +++++++++++
         * dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
         */ 
        dd12 = densityDistributions[i + 12*domainCellsCPUWithHalo]; 
        rho += dd12;
        velocity_y += dd12;
        velocity_z += dd12;

        dd13 = densityDistributions[i + 13*domainCellsCPUWithHalo]; 
        rho += dd13;
        velocity_y -= dd13;
        velocity_z -= dd13;
    
        dd14 = densityDistributions[i + 14*domainCellsCPUWithHalo]; 
        rho += dd14;
        velocity_y += dd14;
        velocity_z -= dd14;

        dd15 = densityDistributions[i + 15*domainCellsCPUWithHalo]; 
        rho += dd15;
        velocity_y -= dd15;
        velocity_z += dd15;
    
        /* +++++++++++
         * +++ DD4 +++
         * +++++++++++
         * dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
         */
        dd16 = densityDistributions[i + 16*domainCellsCPUWithHalo]; 
        rho += dd16;
        velocity_z += dd16;

        dd17 = densityDistributions[i + 17*domainCellsCPUWithHalo]; 
        rho += dd17;
        velocity_z -= dd17;

        dd18 = densityDistributions[i + 18*domainCellsCPUWithHalo]; 
        rho += dd18;
    
    
        /**
         * instead of storing the density distributions after modification,
         * we store it during the modifications to hide memory waiting stalls
         */
    
#define tmp rho
        T dd_param; // modified rho as temporary variable
        switch(flags[i])
        {
            case FLUID:    // this is the whole collision operator
                vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
                dd_param = rho - (T)(3.0f/2.0f)*(vel2);
    
                tmp = gravitation[0]*(T)(1.0f/18.0f)*rho;
                vela2 = velocity_x*velocity_x;
                dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
                dd1 -= tmp;
                densityDistributions[i + 0*domainCellsCPUWithHalo] = dd1;
    
                dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
                dd0 += tmp;
                densityDistributions[i + 1*domainCellsCPUWithHalo] = dd0;
    
                tmp = gravitation[1]*(T)(-1.0f/18.0f)*rho;
                vela2 = velocity_y*velocity_y;
                dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
                dd3 -= tmp;
                densityDistributions[i + 2*domainCellsCPUWithHalo] = dd3;
    
                dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
                dd2 += tmp;
                densityDistributions[i + 3*domainCellsCPUWithHalo] = dd2;
    
#define vela_velb_2 vela2
                /***********************
                 * DD1
                 ***********************/
                vela_velb = velocity_x+velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
    
                tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;
    
                dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
                dd5 -= tmp;
                densityDistributions[i + 4*domainCellsCPUWithHalo] = dd5;
    
                dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
                dd4 += tmp;
                densityDistributions[i + 5*domainCellsCPUWithHalo] = dd4;
    
                vela_velb = velocity_x-velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;
                dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
                dd7 -= tmp;
                densityDistributions[i + 6*domainCellsCPUWithHalo] = dd7;
    
                dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
                dd6 += tmp;
                densityDistributions[i + 7*domainCellsCPUWithHalo] = dd6;
    
                /***********************
                 * DD2
                 ***********************/
                vela_velb = velocity_x+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;
    
                dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
                dd9 -= tmp;
                densityDistributions[i + 8*domainCellsCPUWithHalo] = dd9;
    
                dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
                dd8 += tmp;
                densityDistributions[i + 9*domainCellsCPUWithHalo] = dd8;
    
                tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;
                vela_velb = velocity_x-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
                dd11 -= tmp;
                densityDistributions[i + 10*domainCellsCPUWithHalo] = dd11;
    
                dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
                dd10 += tmp;
                densityDistributions[i + 11*domainCellsCPUWithHalo] = dd10;
    
                /***********************
                 * DD3
                 ***********************/
                vela_velb = velocity_y+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
    
                tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;
                dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
                dd13 -= tmp;
                densityDistributions[i + 12*domainCellsCPUWithHalo] = dd13;
    
                dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
                dd12 += tmp;
                densityDistributions[i + 13*domainCellsCPUWithHalo] = dd12;
    
                vela_velb = velocity_y-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;
    
                dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
                dd15 -= tmp;
                densityDistributions[i + 14*domainCellsCPUWithHalo] = dd15;
    
                dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
                dd14 += tmp;
                densityDistributions[i + 15*domainCellsCPUWithHalo] = dd14;
    
#undef vela_velb_2
                /***********************
                 * DD4
                 ***********************/
                vela2 = velocity_z*velocity_z;
    
                tmp = gravitation[2]*(T)(1.0f/18.0f)*rho;
                dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
                dd17 -= tmp;
                densityDistributions[i + 16*domainCellsCPUWithHalo] = dd17;
    
                dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
                dd16 += tmp;
                densityDistributions[i + 17*domainCellsCPUWithHalo] = dd16;
    
                dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
                densityDistributions[i + 18*domainCellsCPUWithHalo] = dd18;
    
                break;
    
            case OBSTACLE: // in case of an obstacle, we bounce back the values
    
                /**
                 * if we are using only bounce back method, it's not necessary to write back the results.
                 * we even wouldn't need to read the density distribution values.
                 */
    
                if (storeVelocities)
                {
                    velocity_x = (T)0;
                    velocity_y = (T)0;
                    velocity_z = (T)0;
                }
                break;
    
            case VELOCITY_INJECTION:   // this flag specifies the injection area of the fluid
                velocity_x = drivenCavityVelocity;
                velocity_y = (T)0;
                velocity_z = (T)0;
    
                rho = (T)1;
    
                vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
                dd_param = rho - ((T)3/(T)2)*(vel2);
    
                /***********************
                 * DD0
                 ***********************/
                vela2 = velocity_x*velocity_x;
                tmp = gravitation[0]*((T)1/(T)18)*rho;
    
                dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
                dd1 -= tmp;
                densityDistributions[i + 0*domainCellsCPUWithHalo] = dd1;
    
                dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
                dd0 += tmp;
                densityDistributions[i + 1*domainCellsCPUWithHalo] = dd0;
    
                vela2 = velocity_y*velocity_y;
                tmp = gravitation[1]*((T)-1/(T)18)*rho;
    
                dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
                dd3 -= tmp;
                densityDistributions[i + 2*domainCellsCPUWithHalo] = dd3;
    
                dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
                dd2 += tmp;
                densityDistributions[i + 3*domainCellsCPUWithHalo] = dd2;
    
#define vela_velb_2 vela2
                /***********************
                 * DD1
                 ***********************/
                vela_velb = velocity_x+velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;
    
                dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                dd5 -= tmp;
                densityDistributions[i + 4*domainCellsCPUWithHalo] = dd5;
    
                dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd4 += tmp;
                densityDistributions[i + 5*domainCellsCPUWithHalo] = dd4;
    
                vela_velb = velocity_x-velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;
    
                dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                dd7 -= tmp;
                densityDistributions[i + 6*domainCellsCPUWithHalo] = dd7;
    
                dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd6 += tmp;
                densityDistributions[i + 7*domainCellsCPUWithHalo] = dd6;
    
                /***********************
                 * DD2
                 ***********************/
                vela_velb = velocity_x+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;
    
                dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                dd9 -= tmp;
                densityDistributions[i + 8*domainCellsCPUWithHalo] = dd9;
    
                dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd8 += tmp;
                densityDistributions[i + 9*domainCellsCPUWithHalo] = dd8;
    
                vela_velb = velocity_x-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;
    
                dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                dd11 -= tmp;
                densityDistributions[i + 10*domainCellsCPUWithHalo] = dd11;
    
                dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd10 += tmp;
                densityDistributions[i + 11*domainCellsCPUWithHalo] = dd10;
    
                /***********************
                 * DD3
                 ***********************/
                vela_velb = velocity_y+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;
    
                dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                dd13 -= tmp;
                densityDistributions[i + 12*domainCellsCPUWithHalo] = dd13;
    
                dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd12 += tmp;
                densityDistributions[i + 13*domainCellsCPUWithHalo] = dd12;
    
                vela_velb = velocity_y-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;
    
                dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                dd15 -= tmp;
                densityDistributions[i + 14*domainCellsCPUWithHalo] = dd15;
    
                dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd14 += tmp;
                densityDistributions[i + 15*domainCellsCPUWithHalo] = dd14;
#undef vela_velb_2
                /***********************
                 * DD4
                 ***********************/
                vela2 = velocity_z*velocity_z;
    
                tmp = gravitation[2]*(T)(1.0f/18.0f)*rho;
                dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
                dd17 -= tmp;
                densityDistributions[i + 16*domainCellsCPUWithHalo] = dd17;
    
                dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
                dd16 += tmp;
                densityDistributions[i + 17*domainCellsCPUWithHalo] = dd16;
    
                dd18 = eq_dd18(dd_param);
                densityDistributions[i + 18*domainCellsCPUWithHalo] = dd18;
                break;
    
            case (GHOST_LAYER):
                    break;
        }
    
        if (storeVelocities)
        {
            // store velocity
            velocities[i + 0*domainCellsCPUWithHalo] = velocity_x;
            velocities[i + 1*domainCellsCPUWithHalo] = velocity_y;
            velocities[i + 2*domainCellsCPUWithHalo] = velocity_z;
        }
    
        if (storeDensities)
        {
            // store density (not necessary)
            densities[i] = rho;
        }
    }
}

template class CLbmAlphaCPU<double>;
template class CLbmAlphaCPU<float>;
