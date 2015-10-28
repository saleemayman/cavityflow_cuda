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
                    CVector<3, int> domainSize,
                    CVector<3, int> domainSizeGPU,
                    CVector<3, T> gravitation) :
                        domainSize(domainSize),
                        domainSizeGPU(domainSizeGPU),
                        gravitation(gravitation)
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
CLbmAlphaCPU<T>::~CLbmAlphaCPU()
{
}


template<class T>
void CLbmAlphaCPU<T>::alphaKernelCPU(
                T* global_dd,                 // density distributions
                const Flag* flag_array,        // flags
                T* velocity,                  // velocities
                T* density,                   // densities
                const T inv_tau,
                const T drivenCavityVelocity, // velocity parameters for modification of density distributions
                const bool storeDensities,
                const bool storeVelocities)
{
    // linear cell index
    int cellID = 0;

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
                /*
                 * skip cell if it is a ghost cell. Note: a ghost cell also means
                 * that the cell lies at the boundary of the GPU domain and we 
                 * have to skip them as well.
                 */
                if  (flag_array[cellID] == GHOST_LAYER)
                    break;

                /*
                 * we have to sum the densities up in a specific order.
                 * otherwise it seems that we run into numerical errors.
                 */
            
                /* +++++++++++
                 * +++ DD0 +++
                 * +++++++++++
                 * 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
                 */
                dd0 = global_dd[cellID + 0*domainCellsCPU];
                rho = dd0;
                velocity_x = dd0;
               
                dd1 = global_dd[cellID + 1*domainCellsCPU]; 
                rho += dd1;
                velocity_x -= dd1;

                dd2 = global_dd[cellID + 2*domainCellsCPU]; 
                rho += dd2;
                velocity_y = dd2;
            
                dd3 = global_dd[cellID + 3*domainCellsCPU]; 
                rho += dd3;
                velocity_y -= dd3;
            
                /* +++++++++++
                 * +++ DD1 +++
                 * +++++++++++
                 * 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
                 */
                dd4 = global_dd[cellID + 4*domainCellsCPU]; 
                rho += dd4;
                velocity_x += dd4;
                velocity_y += dd4;

                dd5 = global_dd[cellID + 5*domainCellsCPU]; 
                rho += dd5;
                velocity_x -= dd5;
                velocity_y -= dd5;
            
                dd6 = global_dd[cellID + 6*domainCellsCPU]; 
                rho += dd6;
                velocity_x += dd6;
                velocity_y -= dd6;

                dd7 = global_dd[cellID + 7*domainCellsCPU]; 
                rho += dd7;
                velocity_x -= dd7;
                velocity_y += dd7;
            
                /* +++++++++++
                 * +++ DD2 +++
                 * +++++++++++
                 * 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
                 */
                dd8 = global_dd[cellID + 8*domainCellsCPU]; 
                rho += dd8;
                velocity_x += dd8;
                velocity_z = dd8;

                dd9 = global_dd[cellID + 9*domainCellsCPU]; 
                rho += dd9;
                velocity_x -= dd9;
                velocity_z -= dd9;
            
                dd10 = global_dd[cellID + 10*domainCellsCPU]; 
                rho += dd10;
                velocity_x += dd10;
                velocity_z -= dd10;

                dd11 = global_dd[cellID + 11*domainCellsCPU]; 
                rho += dd11;
                velocity_x -= dd11;
                velocity_z += dd11;
            
                /* +++++++++++
                 * +++ DD3 +++
                 * +++++++++++
                 * dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
                 */ 
                dd12 = global_dd[cellID + 12*domainCellsCPU]; 
                rho += dd12;
                velocity_y += dd12;
                velocity_z += dd12;

                dd13 = global_dd[cellID + 13*domainCellsCPU]; 
                rho += dd13;
                velocity_y -= dd13;
                velocity_z -= dd13;
            
                dd14 = global_dd[cellID + 14*domainCellsCPU]; 
                rho += dd14;
                velocity_y += dd14;
                velocity_z -= dd14;

                dd15 = global_dd[cellID + 15*domainCellsCPU]; 
                rho += dd15;
                velocity_y -= dd15;
                velocity_z += dd15;
            
                /* +++++++++++
                 * +++ DD4 +++
                 * +++++++++++
                 * dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
                 */
                dd16 = global_dd[cellID + 16*domainCellsCPU]; 
                rho += dd16;
                velocity_z += dd16;

                dd17 = global_dd[cellID + 17*domainCellsCPU]; 
                rho += dd17;
                velocity_z -= dd17;

                dd18 = global_dd[cellID + 18*domainCellsCPU]; 
                rho += dd18;
            
            
                /**
                 * instead of storing the density distributions after modification,
                 * we store it during the modifications to hide memory waiting stalls
                 */
            
#define tmp rho
                T dd_param; // modified rho as temporary variable
                switch(flag_array[cellID])
                {
                    case FLUID:    // this is the whole collision operator
                        vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
                        dd_param = rho - (T)(3.0f/2.0f)*(vel2);
            
                        tmp = gravitation[0]*(T)(1.0f/18.0f)*rho;
                        vela2 = velocity_x*velocity_x;
                        dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
                        dd1 -= tmp;
                        global_dd[cellID + 0*domainCellsCPU] = dd1;
            
                        dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
                        dd0 += tmp;
                        global_dd[cellID + 1*domainCellsCPU] = dd0;
            
                        tmp = gravitation[1]*(T)(-1.0f/18.0f)*rho;
                        vela2 = velocity_y*velocity_y;
                        dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
                        dd3 -= tmp;
                        global_dd[cellID + 2*domainCellsCPU] = dd3;
            
                        dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
                        dd2 += tmp;
                        global_dd[cellID + 3*domainCellsCPU] = dd2;
            
#define vela_velb_2 vela2
                        /***********************
                         * DD1
                         ***********************/
                        vela_velb = velocity_x+velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
            
                        tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;
            
                        dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
                        dd5 -= tmp;
                        global_dd[cellID + 4*domainCellsCPU] = dd5;
            
                        dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
                        dd4 += tmp;
                        global_dd[cellID + 5*domainCellsCPU] = dd4;
            
                        vela_velb = velocity_x-velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;
                        dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
                        dd7 -= tmp;
                        global_dd[cellID + 6*domainCellsCPU] = dd7;
            
                        dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
                        dd6 += tmp;
                        global_dd[cellID + 7*domainCellsCPU] = dd6;
            
                        /***********************
                         * DD2
                         ***********************/
                        vela_velb = velocity_x+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;
            
                        dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
                        dd9 -= tmp;
                        global_dd[cellID + 8*domainCellsCPU] = dd9;
            
                        dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
                        dd8 += tmp;
                        global_dd[cellID + 9*domainCellsCPU] = dd8;
            
                        tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;
                        vela_velb = velocity_x-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
                        dd11 -= tmp;
                        global_dd[cellID + 10*domainCellsCPU] = dd11;
            
                        dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
                        dd10 += tmp;
                        global_dd[cellID + 11*domainCellsCPU] = dd10;
            
                        /***********************
                         * DD3
                         ***********************/
                        vela_velb = velocity_y+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
            
                        tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;
                        dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
                        dd13 -= tmp;
                        global_dd[cellID + 12*domainCellsCPU] = dd13;
            
                        dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
                        dd12 += tmp;
                        global_dd[cellID + 13*domainCellsCPU] = dd12;
            
                        vela_velb = velocity_y-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;
            
                        dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
                        dd15 -= tmp;
                        global_dd[cellID + 14*domainCellsCPU] = dd15;
            
                        dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
                        dd14 += tmp;
                        global_dd[cellID + 15*domainCellsCPU] = dd14;
            
#undef vela_velb_2
                        /***********************
                         * DD4
                         ***********************/
                        vela2 = velocity_z*velocity_z;
            
                        tmp = gravitation[2]*(T)(1.0f/18.0f)*rho;
                        dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
                        dd17 -= tmp;
                        global_dd[cellID + 16*domainCellsCPU] = dd17;
            
                        dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
                        dd16 += tmp;
                        global_dd[cellID + 17*domainCellsCPU] = dd16;
            
                        dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
                        global_dd[cellID + 18*domainCellsCPU] = dd18;
            
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
                        global_dd[cellID + 0*domainCellsCPU] = dd1;
            
                        dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
                        dd0 += tmp;
                        global_dd[cellID + 1*domainCellsCPU] = dd0;
            
                        vela2 = velocity_y*velocity_y;
                        tmp = gravitation[1]*((T)-1/(T)18)*rho;
            
                        dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
                        dd3 -= tmp;
                        global_dd[cellID + 2*domainCellsCPU] = dd3;
            
                        dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
                        dd2 += tmp;
                        global_dd[cellID + 3*domainCellsCPU] = dd2;
            
#define vela_velb_2 vela2
                        /***********************
                         * DD1
                         ***********************/
                        vela_velb = velocity_x+velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;
            
                        dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd5 -= tmp;
                        global_dd[cellID + 4*domainCellsCPU] = dd5;
            
                        dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd4 += tmp;
                        global_dd[cellID + 5*domainCellsCPU] = dd4;
            
                        vela_velb = velocity_x-velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;
            
                        dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd7 -= tmp;
                        global_dd[cellID + 6*domainCellsCPU] = dd7;
            
                        dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd6 += tmp;
                        global_dd[cellID + 7*domainCellsCPU] = dd6;
            
                        /***********************
                         * DD2
                         ***********************/
                        vela_velb = velocity_x+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;
            
                        dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd9 -= tmp;
                        global_dd[cellID + 8*domainCellsCPU] = dd9;
            
                        dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd8 += tmp;
                        global_dd[cellID + 9*domainCellsCPU] = dd8;
            
                        vela_velb = velocity_x-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;
            
                        dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd11 -= tmp;
                        global_dd[cellID + 10*domainCellsCPU] = dd11;
            
                        dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd10 += tmp;
                        global_dd[cellID + 11*domainCellsCPU] = dd10;
            
                        /***********************
                         * DD3
                         ***********************/
                        vela_velb = velocity_y+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;
            
                        dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd13 -= tmp;
                        global_dd[cellID + 12*domainCellsCPU] = dd13;
            
                        dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd12 += tmp;
                        global_dd[cellID + 13*domainCellsCPU] = dd12;
            
                        vela_velb = velocity_y-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;
            
                        dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd15 -= tmp;
                        global_dd[cellID + 14*domainCellsCPU] = dd15;
            
                        dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd14 += tmp;
                        global_dd[cellID + 15*domainCellsCPU] = dd14;
#undef vela_velb_2
                        /***********************
                         * DD4
                         ***********************/
                        vela2 = velocity_z*velocity_z;
            
                        tmp = gravitation[2]*(T)(1.0f/18.0f)*rho;
                        dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
                        dd17 -= tmp;
                        global_dd[cellID + 16*domainCellsCPU] = dd17;
            
                        dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
                        dd16 += tmp;
                        global_dd[cellID + 17*domainCellsCPU] = dd16;
            
                        dd18 = eq_dd18(dd_param);
                        global_dd[cellID + 18*domainCellsCPU] = dd18;
                        break;
            
                    case (GHOST_LAYER):
                            break;
                }
            
                if (storeVelocities)
                {
                    // store velocity
                    velocity[cellID + 0*domainCellsCPU] = velocity_x;
                    velocity[cellID + 1*domainCellsCPU] = velocity_y;
                    velocity[cellID + 2*domainCellsCPU] = velocity_z;
                }
            
                if (storeDensities)
                {
                    // store density (not necessary)
                    density[cellID] = rho;
                }

                // increment the linear index of domain cells
                cellID++;
            }
        }
    }
}

template class CLbmAlphaCPU<double>;
template class CLbmAlphaCPU<float>;
