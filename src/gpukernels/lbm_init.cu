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

#include "lbm_init.cuh"

__device__ int3 getCubePosition(int linear_position, const int domainCells_x, const int domainCells_y)
{
    int3 pos;

    /*
     * TODO
     * use AND operation to speed up (but this function is only used during initialization)
     */
    pos.x = linear_position % domainCells_x;
    linear_position /= domainCells_x;

    pos.y = linear_position % domainCells_y;
    linear_position /= domainCells_y;

    pos.z = linear_position;// % CUBE_SIZE_Z;

    return pos;
}

template<typename T>
__global__ void lbm_init(
        T *global_dd,           // density distributions
        Flag *flags,            // flags
        T *velocity_array,      // velocity array (first all x components, then all y components, then z...)
        T *density,             // densities
        Flag boundaryConditionRight,
        Flag boundaryConditionLeft,
        Flag boundaryConditionTop,
        Flag boundaryConditionBottom,
        Flag boundaryConditionFront,
        Flag boundaryConditionBack,
        T drivenCavityVelocity, // velocity parameters for modification of density distributions
        const int domainCells_x,
        const int domainCells_y,
        const int domainCells_z,
        const bool storeDensities,
        const bool storeVelocities)
{
    size_t DOMAIN_CELLS = domainCells_x * domainCells_y * domainCells_z;

    // get unique global ID
    size_t blockId = blockIdx.x + (size_t)(blockIdx.y * gridDim.x) + (size_t)(gridDim.x * gridDim.y * blockIdx.z);
    size_t gid = blockId * (size_t)(blockDim.x * blockDim.y * blockDim.z) + (size_t)(threadIdx.z * (blockDim.x * blockDim.y)) + (size_t)(threadIdx.y * blockDim.x) + threadIdx.x;

    if (gid >= DOMAIN_CELLS)
        return;

    //__global T *current_dds = &global_dd[gid];
    T *current_dds = &global_dd[gid];

    // initialize flag field
    int3 pos;
    pos = getCubePosition(gid, domainCells_x, domainCells_y);

    T velocity_x = 0;
    T velocity_y = 0;
    T velocity_z = 0;

    Flag flag = FLUID;

    if(pos.x == 0)
        flag = boundaryConditionRight;
    if(pos.x == domainCells_x-1 && flag != OBSTACLE)
        flag = boundaryConditionLeft;
    if(pos.y == 0 && flag != OBSTACLE)
        flag = boundaryConditionTop;
    if(pos.y == domainCells_y-1 && flag != OBSTACLE)
        flag = boundaryConditionBottom;
    if(pos.z == 0 && flag != OBSTACLE)
        flag = boundaryConditionFront;
    if(pos.z == domainCells_z-1 && flag != OBSTACLE)
        flag = boundaryConditionBack;

    // density distributions
    T dd0, dd1, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9, dd10, dd11, dd12, dd13, dd14, dd15, dd16, dd17, dd18;

    T dd_param;
    T vela2;
    T vela_velb;
    T rho = 1.0f;

    // compute and store velocity

    vela2 = velocity_x*velocity_x;
    dd_param = rho - (T)(3.0f/2.0f)*(vela2);

    dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
    *current_dds = dd0;     current_dds += DOMAIN_CELLS;
    dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
    *current_dds = dd1;     current_dds += DOMAIN_CELLS;

    vela2 = velocity_y*velocity_y;

    dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
    *current_dds = dd2;     current_dds += DOMAIN_CELLS;
    dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
    *current_dds = dd3;     current_dds += DOMAIN_CELLS;


#define vela_velb_2 vela2
    /***********************
     * DD1
     ***********************/
    vela_velb = velocity_x+velocity_y;
    vela_velb_2 = vela_velb*vela_velb;

    dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd4;     current_dds += DOMAIN_CELLS;
    dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd5;     current_dds += DOMAIN_CELLS;

    vela_velb = velocity_x-velocity_y;
    vela_velb_2 = vela_velb*vela_velb;

    dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd6;     current_dds += DOMAIN_CELLS;
    dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd7;     current_dds += DOMAIN_CELLS;

    /***********************
     * DD2
     ***********************/
    vela_velb = velocity_x+velocity_z;
    vela_velb_2 = vela_velb*vela_velb;

    dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd8;     current_dds += DOMAIN_CELLS;
    dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd9;     current_dds += DOMAIN_CELLS;

    vela_velb = velocity_x-velocity_z;
    vela_velb_2 = vela_velb*vela_velb;

    dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd10;        current_dds += DOMAIN_CELLS;
    dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd11;        current_dds += DOMAIN_CELLS;

    /***********************
     * DD3
     ***********************/
    vela_velb = velocity_y+velocity_z;
    vela_velb_2 = vela_velb*vela_velb;


    dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd12;        current_dds += DOMAIN_CELLS;
    dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd13;        current_dds += DOMAIN_CELLS;

    vela_velb = velocity_y-velocity_z;
    vela_velb_2 = vela_velb*vela_velb;

    dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd14;        current_dds += DOMAIN_CELLS;
    dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    *current_dds = dd15;        current_dds += DOMAIN_CELLS;


#undef vela_velb_2
    /***********************
     * DD4
     ***********************/
    vela2 = velocity_z*velocity_z;

    dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
    *current_dds = dd16;        current_dds += DOMAIN_CELLS;
    dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
    *current_dds = dd17;        current_dds += DOMAIN_CELLS;

    dd18 = eq_dd18(dd_param);
    *current_dds = dd18;

    // flag
    flags[gid] = flag;

    if (storeVelocities)
    {
		current_dds = &velocity_array[gid];
		*current_dds = velocity_x;  current_dds += DOMAIN_CELLS;
		*current_dds = velocity_y;  current_dds += DOMAIN_CELLS;
		*current_dds = velocity_z;
    }

    if (storeDensities)
        density[gid] = rho;
}

template __global__ void lbm_init<double>(
        double *global_dd,
        Flag *flags,
        double *velocity_array,
        double *density,
        Flag boundaryConditionRight,
        Flag boundaryConditionLeft,
        Flag boundaryConditionTop,
        Flag boundaryConditionBottom,
        Flag boundaryConditionFront,
        Flag boundaryConditionBack,
        double drivenCavityVelocity,
        const int domainCells_x,
        const int domainCells_y,
        const int domainCells_z,
        const bool storeDensities,
        const bool storeVelocities);
template __global__ void lbm_init<float>(
        float *global_dd,
        Flag *flags,
        float *velocity_array,
        float *density,
        Flag boundaryConditionRight,
        Flag boundaryConditionLeft,
        Flag boundaryConditionTop,
        Flag boundaryConditionBottom,
        Flag boundaryConditionFront,
        Flag boundaryConditionBack,
        float drivenCavityVelocity,
        const int domainCells_x,
        const int domainCells_y,
        const int domainCells_z,
        const bool storeDensities,
        const bool storeVelocities);
