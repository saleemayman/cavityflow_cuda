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

#include "mpi.h"

#include "libmath/CVector.hpp"
#include "constants.h"

extern CVector<3,int> E0;
extern CVector<3,int> E1;
extern CVector<3,int> E2;
extern CVector<3,int> E3;
extern CVector<3,int> E4;
extern CVector<3,int> E5;
extern CVector<3,int> E6;
extern CVector<3,int> E7;
extern CVector<3,int> E8;
extern CVector<3,int> E9;
extern CVector<3,int> E10;
extern CVector<3,int> E11;
extern CVector<3,int> E12;
extern CVector<3,int> E13;
extern CVector<3,int> E14;
extern CVector<3,int> E15;
extern CVector<3,int> E16;
extern CVector<3,int> E17;
extern CVector<3,int> E18;
extern CVector<3,int> lbm_units[];

// typedef Singleton<CConfiguration<T> > ConfigSingleton; // Global declaration
// typedef Singleton<CProfiler> ProfilerSingleton; // Global declaration

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

#ifndef LOG_TO_FILE
#define  DEBUGPRINT(...) \
{ \
    int my_rank; \
    int num_procs; \
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); \
    fprintf(stderr, "P %d: ", my_rank); \
    fprintf(stderr, __VA_ARGS__); \
}
#else
#define  DEBUGPRINT(...)            \
{ \
    int my_rank; \
    int num_procs;                           \
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); \
    std::string outputfilename = LOG_OUTPUT_DIR; \
    std::stringstream ss_file; \
    ss_file << "./"  << LOG_OUTPUT_DIR  << "/" << LOG_OUTPUT_FILE_PREFIX << "_" << my_rank << ".log";    \
    std::string outputfile = ss_file.str(); \
    FILE* file = fopen(outputfile.c_str(),"a");    \
    fprintf(file, "P %d: ", my_rank); \
    fprintf(file, __VA_ARGS__); \
    file.close() \
}
#endif

inline const char * mpiGetErrorString(int errorNum)
{
    switch (errorNum) {
    case MPI_SUCCESS:
        return "CL_SUCCESS";
    case MPI_ERR_REQUEST:
        return "MPI_ERR_REQUEST";
    case MPI_ERR_ARG:
        return "MPI_ERR_ARG";
    }
    return "UNKNOWN";
}

#define MPI_CHECK_ERROR(val_a)    if ((val_a) != MPI_SUCCESS)        \
{                                \
    std::cerr    << "MPI_ERROR: file: '" << __FILE__    \
                << "', line: " << __LINE__        \
                << " - " << mpiGetErrorString(val_a) << std::endl;    \
}

#endif
