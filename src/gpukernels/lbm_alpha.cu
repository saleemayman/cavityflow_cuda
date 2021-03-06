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

#include "lbm_alpha.cuh"

template<typename T>
__global__ void lbm_kernel_alpha(
        T* global_dd,                 // density distributions
        const Flag* flag_array,       // flags
        T* velocity,                  // velocities
        T* density,                   // densities
        const T inv_tau,
        const T gravitation_x,
        const T gravitation_y,
        const T gravitation_z,
        const T drivenCavityVelocity, // velocity parameters for modification of density distributions
        const int originX,
        const int originY,
        const int originZ,
        const int sizeX,
        const int sizeY,
        const int sizeZ,
        const int domainCellsX,
        const int domainCellsY,
        const int domainCellsZ,
        const bool storeDensities,
        const bool storeVelocities)
{
    const unsigned int DOMAIN_CELLS = domainCellsX * domainCellsY * domainCellsZ;

    const unsigned int X = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int Y = blockIdx.y * blockDim.y + threadIdx.y;
    const unsigned int Z = blockIdx.z * blockDim.z + threadIdx.z;
    const unsigned int gid = (originZ + Z) * (domainCellsX * domainCellsY) + (originY + Y) * domainCellsX + (originX + X);

    if (gid >= DOMAIN_CELLS)
        return;

    if (X >= sizeX || Y >= sizeY || Z >= sizeZ)
        return;

    Flag flag = flag_array[gid];
    if  (flag == GHOST_LAYER)
        return;

    /**
     * we use a pointer instead of accessing the array directly
     * first this reduces the number of use registers (according to profiling information)
     * secondly the program runs faster and we can use more threads
     */

    // pointer to density distributions
    //__global T *current_dds = &global_dd[gid];
    T *current_dds = &global_dd[gid];

    // velocity
    T velocity_x, velocity_y, velocity_z;

    // density distributions
    T dd0, dd1, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9, dd10, dd11, dd12, dd13, dd14, dd15, dd16, dd17, dd18;

    // density
    T rho;

    /*
     * we have to sum the densities up in a specific order.
     * otherwise it seems that we run into numerical errors.
     */

    // +++++++++++
    // +++ DD0 +++
    // +++++++++++
    //
    // 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
    dd0 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho = dd0;
    velocity_x = dd0;
    dd1 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd1;
    velocity_x -= dd1;

    dd2 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd2;
    velocity_y = dd2;
    dd3 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd3;
    velocity_y -= dd3;

    // +++++++++++
    // +++ DD1 +++
    // +++++++++++
    //
    // 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
    dd4 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd4;
    velocity_x += dd4;
    velocity_y += dd4;
    dd5 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd5;
    velocity_x -= dd5;
    velocity_y -= dd5;

    dd6 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd6;
    velocity_x += dd6;
    velocity_y -= dd6;
    dd7 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd7;
    velocity_x -= dd7;
    velocity_y += dd7;

    // +++++++++++
    // +++ DD2 +++
    // +++++++++++
    //
    // 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
    dd8 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd8;
    velocity_x += dd8;
    velocity_z = dd8;
    dd9 = *current_dds;     current_dds += DOMAIN_CELLS;
    rho += dd9;
    velocity_x -= dd9;
    velocity_z -= dd9;

    dd10 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd10;
    velocity_x += dd10;
    velocity_z -= dd10;
    dd11 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd11;
    velocity_x -= dd11;
    velocity_z += dd11;

    // +++++++++++
    // +++ DD3 +++
    // +++++++++++

    // dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
    dd12 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd12;
    velocity_y += dd12;
    velocity_z += dd12;
    dd13 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd13;
    velocity_y -= dd13;
    velocity_z -= dd13;

    dd14 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd14;
    velocity_y += dd14;
    velocity_z -= dd14;
    dd15 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd15;
    velocity_y -= dd15;
    velocity_z += dd15;

    // +++++++++++
    // +++ DD4 +++
    // +++++++++++
    //
    // dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
    dd16 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd16;
    velocity_z += dd16;
    dd17 = *current_dds;        current_dds += DOMAIN_CELLS;
    rho += dd17;
    velocity_z -= dd17;
    dd18 = *current_dds;
    rho += dd18;

//  rho *= (float)(flag != FLAG_OBSTACLE);


    // to add something to a pointer is faster than subtracting it.
    // thus we restart here with pointer to dd0
    current_dds = &global_dd[gid];

    /**
     * instead of storing the density distributions after modification,
     * we store it during the modifications to hide memory waiting stalls
     */

    T vel2;     // vel*vel
    T vela2;
    T vela_velb;
#define tmp rho


    T dd_param; // modified rho as temporary variable
    switch(flag)
    {
        case FLUID:    // this is the whole collision operator
            vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
            dd_param = rho - ((T)3/(T)2)*vel2;

            tmp = gravitation_x*((T)1/(T)18)*rho;
            vela2 = velocity_x*velocity_x;
            dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
            dd1 -= tmp;
            *current_dds = dd1;     current_dds += DOMAIN_CELLS;

            dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
            dd0 += tmp;
            *current_dds = dd0;     current_dds += DOMAIN_CELLS;

            tmp = gravitation_y*((T)-1/(T)18)*rho;
            vela2 = velocity_y*velocity_y;
            dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
            dd3 -= tmp;
            *current_dds = dd3;     current_dds += DOMAIN_CELLS;

            dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
            dd2 += tmp;
            *current_dds = dd2;     current_dds += DOMAIN_CELLS;


#define vela_velb_2 vela2
            /***********************
             * DD1
             ***********************/
            vela_velb = velocity_x+velocity_y;
            vela_velb_2 = vela_velb*vela_velb;

            tmp = (gravitation_x - gravitation_y)*((T)1/(T)36)*rho;

            dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
            dd5 -= tmp;
            *current_dds = dd5;     current_dds += DOMAIN_CELLS;

            dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
            dd4 += tmp;
            *current_dds = dd4;     current_dds += DOMAIN_CELLS;

            vela_velb = velocity_x-velocity_y;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_x + gravitation_y)*((T)1/(T)36)*rho;
            dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
            dd7 -= tmp;
            *current_dds = dd7;     current_dds += DOMAIN_CELLS;

            dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
            dd6 += tmp;
            *current_dds = dd6;     current_dds += DOMAIN_CELLS;

            /***********************
             * DD2
             ***********************/
            vela_velb = velocity_x+velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_x + gravitation_z)*((T)1/(T)36)*rho;

            dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
            dd9 -= tmp;
            *current_dds = dd9;     current_dds += DOMAIN_CELLS;

            dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
            dd8 += tmp;
            *current_dds = dd8;     current_dds += DOMAIN_CELLS;

            tmp = (gravitation_x - gravitation_z)*((T)1/(T)36)*rho;
            vela_velb = velocity_x-velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
            dd11 -= tmp;
            *current_dds = dd11;        current_dds += DOMAIN_CELLS;

            dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
            dd10 += tmp;
            *current_dds = dd10;        current_dds += DOMAIN_CELLS;

            /***********************
             * DD3
             ***********************/
            vela_velb = velocity_y+velocity_z;
            vela_velb_2 = vela_velb*vela_velb;

            tmp = (gravitation_z - gravitation_y)*((T)1/(T)36)*rho;
            dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
            dd13 -= tmp;
            *current_dds = dd13;        current_dds += DOMAIN_CELLS;

            dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
            dd12 += tmp;
            *current_dds = dd12;        current_dds += DOMAIN_CELLS;

            vela_velb = velocity_y-velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_z + gravitation_y)*((T)-1/(T)36)*rho;

            dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
            dd15 -= tmp;
            *current_dds = dd15;        current_dds += DOMAIN_CELLS;

            dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
            dd14 += tmp;
            *current_dds = dd14;        current_dds += DOMAIN_CELLS;

#undef vela_velb_2
            /***********************
             * DD4
             ***********************/
            vela2 = velocity_z*velocity_z;

            tmp = gravitation_z*((T)1/(T)18)*rho;
            dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
            dd17 -= tmp;
            *current_dds = dd17;        current_dds += DOMAIN_CELLS;

            dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
            dd16 += tmp;
            *current_dds = dd16;        current_dds += DOMAIN_CELLS;

            dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
            *current_dds = dd18;

            break;

        case OBSTACLE: // in case of an obstacle, we bounce back the values

            /**
             * if we are using only bounce back method, it's not necessary to write back the results.
             * we even wouldn't need to read the density distribution values.
             */
#if 0
            // use simple bounce back
            *current_dds = dd0;     current_dds += DOMAIN_CELLS;
            *current_dds = dd1;     current_dds += DOMAIN_CELLS;
            *current_dds = dd2;     current_dds += DOMAIN_CELLS;
            *current_dds = dd3;     current_dds += DOMAIN_CELLS;

            *current_dds = dd4;     current_dds += DOMAIN_CELLS;
            *current_dds = dd5;     current_dds += DOMAIN_CELLS;
            *current_dds = dd6;     current_dds += DOMAIN_CELLS;
            *current_dds = dd7;     current_dds += DOMAIN_CELLS;

            *current_dds = dd8;     current_dds += DOMAIN_CELLS;
            *current_dds = dd9;     current_dds += DOMAIN_CELLS;
            *current_dds = dd10;        current_dds += DOMAIN_CELLS;
            *current_dds = dd11;        current_dds += DOMAIN_CELLS;

            *current_dds = dd12;        current_dds += DOMAIN_CELLS;
            *current_dds = dd13;        current_dds += DOMAIN_CELLS;
            *current_dds = dd14;        current_dds += DOMAIN_CELLS;
            *current_dds = dd15;        current_dds += DOMAIN_CELLS;

            *current_dds = dd16;        current_dds += DOMAIN_CELLS;
            *current_dds = dd17;        current_dds += DOMAIN_CELLS;
            *current_dds = dd18;
#endif

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
            tmp = gravitation_x*((T)1/(T)18)*rho;

            dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
            dd1 -= tmp;
            *current_dds = dd1;     current_dds += DOMAIN_CELLS;

            dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
            dd0 += tmp;
            *current_dds = dd0;     current_dds += DOMAIN_CELLS;

            vela2 = velocity_y*velocity_y;
            tmp = gravitation_y*((T)-1/(T)18)*rho;

            dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
            dd3 -= tmp;
            *current_dds = dd3;     current_dds += DOMAIN_CELLS;

            dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
            dd2 += tmp;
            *current_dds = dd2;     current_dds += DOMAIN_CELLS;


#define vela_velb_2 vela2
            /***********************
             * DD1
             ***********************/
            vela_velb = velocity_x+velocity_y;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_x - gravitation_y)*((T)1/(T)36)*rho;

            dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
            dd5 -= tmp;
            *current_dds = dd5;     current_dds += DOMAIN_CELLS;

            dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
            dd4 += tmp;
            *current_dds = dd4;     current_dds += DOMAIN_CELLS;

            vela_velb = velocity_x-velocity_y;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_x + gravitation_y)*((T)1/(T)36)*rho;

            dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
            dd7 -= tmp;
            *current_dds = dd7;     current_dds += DOMAIN_CELLS;

            dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
            dd6 += tmp;
            *current_dds = dd6;     current_dds += DOMAIN_CELLS;

            /***********************
             * DD2
             ***********************/
            vela_velb = velocity_x+velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_x + gravitation_z)*((T)1/(T)36)*rho;

            dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
            dd9 -= tmp;
            *current_dds = dd9;     current_dds += DOMAIN_CELLS;

            dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
            dd8 += tmp;
            *current_dds = dd8;     current_dds += DOMAIN_CELLS;

            vela_velb = velocity_x-velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_x - gravitation_z)*((T)1/(T)36)*rho;

            dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
            dd11 -= tmp;
            *current_dds = dd11;        current_dds += DOMAIN_CELLS;

            dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
            dd10 += tmp;
            *current_dds = dd10;        current_dds += DOMAIN_CELLS;

            /***********************
             * DD3
             ***********************/
            vela_velb = velocity_y+velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_z - gravitation_y)*((T)1/(T)36)*rho;

            dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
            dd13 -= tmp;
            *current_dds = dd13;        current_dds += DOMAIN_CELLS;

            dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
            dd12 += tmp;
            *current_dds = dd12;        current_dds += DOMAIN_CELLS;

            vela_velb = velocity_y-velocity_z;
            vela_velb_2 = vela_velb*vela_velb;
            tmp = (gravitation_z + gravitation_y)*((T)-1/(T)36)*rho;

            dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
            dd15 -= tmp;
            *current_dds = dd15;        current_dds += DOMAIN_CELLS;

            dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
            dd14 += tmp;
            *current_dds = dd14;        current_dds += DOMAIN_CELLS;
#undef vela_velb_2

            /***********************
             * DD4
             ***********************/
            vela2 = velocity_z*velocity_z;

            tmp = gravitation_z*((T)1/(T)18)*rho;
            dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
            dd17 -= tmp;
            *current_dds = dd17;        current_dds += DOMAIN_CELLS;

            dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
            dd16 += tmp;
            *current_dds = dd16;        current_dds += DOMAIN_CELLS;

            dd18 = eq_dd18(dd_param);
            *current_dds = dd18;
            break;

        case (GHOST_LAYER):
                break;
    }

    if (storeVelocities)
    {
        // store velocity
        current_dds = &velocity[gid];
        *current_dds += (T)1;  current_dds += DOMAIN_CELLS;
        *current_dds += (T)1;  current_dds += DOMAIN_CELLS;
        *current_dds += (T)1;
    }

    if (storeDensities)
    {
        // store density (not necessary)
        density[gid] = rho;
        // density[gid] = flag;
    }
}

template __global__ void lbm_kernel_alpha<float>(
        float *global_dd,
        const Flag *flag_array,
        float *velocity,
        float *density,
        const float inv_tau,
        const float gravitation_x,
        const float gravitation_y,
        const float gravitation_z,
        const float drivenCavityVelocity,
        const int originX,
        const int originY,
        const int originZ,
        const int sizeX,
        const int sizeY,
        const int sizeZ,
        const int domainCellsX,
        const int domainCellsY,
        const int domainCellsZ,
        const bool storeDensities,
        const bool storeVelocities);
template __global__ void lbm_kernel_alpha<double>(
        double *global_dd,
        const Flag *flag_array,
        double *velocity,
        double *density,
        const double inv_tau,
        const double gravitation_x,
        const double gravitation_y,
        const double gravitation_z,
        const double drivenCavityVelocity,
        const int originX,
        const int originY,
        const int originZ,
        const int sizeX,
        const int sizeY,
        const int sizeZ,
        const int domainCellsX,
        const int domainCellsY,
        const int domainCellsZ,
        const bool storeDensities,
        const bool storeVelocities);
