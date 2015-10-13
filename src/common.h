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

#ifndef COMMON_H
#define COMMON_H

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include <cuda.h>
#include <cuda_runtime.h>

#include <mpi.h>

#include "libmath/CVector.hpp"
#include "constants.h"

#define E0  CVector<3, int>( 1,  0,  0)
#define E1  CVector<3, int>(-1,  0,  0)
#define E2  CVector<3, int>( 0,  1,  0)
#define E3  CVector<3, int>( 0, -1,  0)
#define E4  CVector<3, int>( 1,  1,  0)
#define E5  CVector<3, int>(-1, -1,  0)
#define E6  CVector<3, int>( 1, -1,  0)
#define E7  CVector<3, int>(-1,  1,  0)
#define E8  CVector<3, int>( 1,  0,  1)
#define E9  CVector<3, int>(-1,  0, -1)
#define E10 CVector<3, int>( 1,  0, -1)
#define E11 CVector<3, int>(-1,  0,  1)
#define E12 CVector<3, int>( 0,  1,  1)
#define E13 CVector<3, int>( 0, -1, -1)
#define E14 CVector<3, int>( 0,  1, -1)
#define E15 CVector<3, int>( 0, -1,  1)
#define E16 CVector<3, int>( 0,  0,  1)
#define E17 CVector<3, int>( 0,  0, -1)
#define E18 CVector<3, int>( 0,  0,  0)

extern CVector<3, int> lbm_units[];

enum Flag
{
    OBSTACLE           = (1<<0),
    FLUID              = (1<<1),
    VELOCITY_INJECTION = (1<<2),
    GHOST_LAYER        = (1<<3),
    GHOST_LAYER_BETA   = (GHOST_LAYER | (1<<4))
};

enum Direction
{
    RIGHT  = 0,
    LEFT   = 1,
    TOP    = 2,
    BOTTOM = 3,
    FRONT  = 4,
    BACK   = 5
};

#define GPU_ERROR_CHECK(code) \
{ \
    gpuErrorCheck((code), __FILE__, __LINE__); \
}

inline void gpuErrorCheck(cudaError_t code, std::string file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        std::cerr << "----- !!! The following CUDA API error occurred !!! -----" << std::endl;
        std::cerr << "error code:     " << code << std::endl;
        std::cerr << "error string:   " << cudaGetErrorString(code) << std::endl;
        std::cerr << "error location: " << std::endl;
        std::cerr << "          file: " << file << std::endl;
        std::cerr << "          line: " << line << std::endl;
        std::cerr << "---------------------------------------------------------" << std::endl;

        if(abort) {
            exit (code);
        }
    }
}

inline dim3 getBlocksPerGrid(int dim, CVector<3,int> size, dim3 threadsPerBlock)
{
    dim3 blocksPerGrid(1, 1, 1);

    switch(dim)
    {
    case 1:
        blocksPerGrid.x = ((size[0] - 1) / threadsPerBlock.x) + 1;
        break;
    case 2:
        blocksPerGrid.x = ((size[0] - 1) / threadsPerBlock.x) + 1;
        blocksPerGrid.y = ((size[1] - 1) / threadsPerBlock.y) + 1;
        break;
    case 3:
        blocksPerGrid.x = ((size[0] - 1) / threadsPerBlock.x) + 1;
        blocksPerGrid.y = ((size[1] - 1) / threadsPerBlock.y) + 1;
        blocksPerGrid.z = ((size[2] - 1) / threadsPerBlock.z) + 1;
        break;
    default:
        blocksPerGrid.x = ((size[0] - 1) / threadsPerBlock.x) + 1;
        blocksPerGrid.y = ((size[1] - 1) / threadsPerBlock.y) + 1;
        blocksPerGrid.z = ((size[2] - 1) / threadsPerBlock.z) + 1;
    }

    return blocksPerGrid;
}

typedef enum {
    MPI_COMM_DIRECTION_UNKNOWN = 0,
    MPI_COMM_DIRECTION_X,
    MPI_COMM_DIRECTION_X_0,
    MPI_COMM_DIRECTION_X_1,
    MPI_COMM_DIRECTION_Y,
    MPI_COMM_DIRECTION_Y_0,
    MPI_COMM_DIRECTION_Y_1,
    MPI_COMM_DIRECTION_Z,
    MPI_COMM_DIRECTION_Z_0,
    MPI_COMM_DIRECTION_Z_1,
    MPI_COMM_DIRECTION_ALL,
} MPI_COMM_DIRECTION;

inline const char* get_string_direction(MPI_COMM_DIRECTION direction) {
    std::string dir;
    switch (direction) {
    case MPI_COMM_DIRECTION_X:
        dir = "X";
        break;
    case MPI_COMM_DIRECTION_X_0:
        dir = "X0";
        break;
    case MPI_COMM_DIRECTION_X_1:
        dir = "X1";
        break;
    case MPI_COMM_DIRECTION_Y:
        dir = "Y";
        break;
    case MPI_COMM_DIRECTION_Y_0:
        dir = "Y0";
        break;
    case MPI_COMM_DIRECTION_Y_1:
        dir = "Y1";
        break;
    case MPI_COMM_DIRECTION_Z:
        dir = "Z";
        break;
    case MPI_COMM_DIRECTION_Z_0:
        dir = "Z0";
        break;
    case MPI_COMM_DIRECTION_Z_1:
        dir = "Z1";
        break;
    default:
        dir = "UNKNOWN";
        break;
    }
    return dir.c_str();
}

/*
 * we can reuse the following function because of its symmetry
 * f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0) f(0,0,1) f(0,0,-1)
 */
#define eq_dd_a0(vela, vela2, rho_alpha) \
    ((1.0/18.0)*((rho_alpha) + (3.0)*(vela) + (9.0/2.0)*(vela2)))
#define eq_dd_a1(vela, vela2, rho_alpha) \
    ((1.0/18.0)*((rho_alpha) + (-3.0)*(vela) + (9.0/2.0)*(vela2)))

/*
 * we can reuse the following functions because of the symmetry of the density distributions!
 *
 * f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
 * f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
 * f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
 */
#define eq_dd4(velx_add_vely, velx_add_vely_2, rho_alpha) \
    ((1.0/36.0)*((rho_alpha) + (3.0)*(velx_add_vely) + (9.0/2.0)*(velx_add_vely_2)))

#define eq_dd5(velx_add_vely, velx_add_vely_2, rho_alpha) \
    ((1.0/36.0)*((rho_alpha) + (-3.0)*(velx_add_vely) + (9.0/2.0)*(velx_add_vely_2)))

#define eq_dd6(velx_sub_vely, velx_sub_vely_2, rho_alpha) \
    ((1.0/36.0)*((rho_alpha) + (3.0)*(velx_sub_vely) + (9.0/2.0)*(velx_sub_vely_2)))

#define eq_dd7(velx_sub_vely, velx_sub_vely_2, rho_alpha) \
    ((1.0/36.0)*((rho_alpha) + (-3.0)*(velx_sub_vely) + (9.0/2.0)*(velx_sub_vely_2)))

/*
 * f(0,0,0)
 */
#define eq_dd18(rho_alpha) \
    ((1.0/3.0)*(rho_alpha))

#endif
