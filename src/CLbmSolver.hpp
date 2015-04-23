/*
 * Copyright 2010 Martin Schreiber
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

#ifndef CLBMOPENCL_HH
#define CLBMOPENCL_HH

#include "CLbmSkeleton.hpp"
#include "libcl/CCL.hpp"
#include "lib/CError.hpp"
#include "libmath/CVector.hpp"
#include <typeinfo>
#include <iomanip>
#include <list>

#include "common.h"

#define LBM_FLAG_OBSTACLE           (1<<0)
#define LBM_FLAG_FLUID              (1<<1)
#define LBM_FLAG_VELOCITY_INJECTION (1<<2)
#define LBM_FLAG_GHOST_LAYER        (1<<3)
#define CU_TRUE     true
#define CU_FALSE    !CU_TRUE

/**
 * the density distributions are packed cell wise to optimize the collision operation and stored
 * linearly (first x, then y, then z) in the buffer cMemDensityDistributions
 *
 * this implementation is using the D3Q19 implementation
 *
 *
 * D3Q19 Vector enumerations + directions
 *
 * f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
 * f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
 * f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
 * f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
 * f(0,0,1), f(0,0,-1),  f(0,0,0),  X
 *
 * this enumeration makes the computation of the equilibrium distribution more efficient
 *
 * velocity(vx, vy, vz, phi)
 *
 *
 *       14,12
 *     15  |  18
 *       \ | /
 *        \|/
 * 7,9----- ----- 8,10
 *        /|\
 *       / | \
 *     17  |  16
 *       11,13
 *
 */

template<typename T>
class CLbmSolver: public CLbmSkeleton<T>
{
public:
    CVector<4, T> drivenCavityVelocity;
    static const size_t SIZE_DD_HOST = 19;
private:
    int _UID;
    // CVector<3, int> _size; ///< Size of the domain
    // CVector<3, int> _subdomain_nums; ///< total number of subdomains

    static const size_t SIZE_DD_HOST_BYTES = SIZE_DD_HOST * sizeof(T);
    int _BC[3][2]; ///< Boundary conditions. First index specifys the dimension and second the upper or the lower boundary.
    // opencl handlers
    CCL::CPlatform cudaPlatform;        ///< constructor will call cuInit
    CCL::CCommandQueue &cCommandQueue;
    CCL::CContext &cContext;
    CCL::CDevice &cDevice;
    CCL::CDeviceInfo cDeviceInfo;
//  OPENCL_VERSION _cl_version;     // check if you need to verify the CUDA driver/toolkit version

    // INITIALIZATION KERNEL
    CCL::CKernel cKernelInit;
    size_t cKernelInit_WorkGroupSize;
    size_t cKernelInit_GlobalWorkGroupSize;
    size_t cKernelInit_MaxRegisters;

    // COLLISION KERNELS
    CCL::CKernel cLbmKernelAlpha;
    size_t cLbmKernelAlpha_WorkGroupSize;
    size_t cLbmKernelAlpha_GlobalWorkGroupSize;
    size_t cLbmKernelAlpha_MaxRegisters;

    CCL::CKernel cLbmKernelBeta;
    size_t cLbmKernelBeta_WorkGroupSize;
    size_t cLbmKernelBeta_GlobalWorkGroupSize;
    size_t cLbmKernelBeta_MaxRegisters;

    // INITIALIZATION KERNEL
    CCL::CKernel cKernelCopyRect;
    size_t cKernelCopyRect_WorkGroupSize;
    size_t cKernelCopyRect_GlobalWorkGroupSize;
    size_t cKernelCopyRect_MaxRegisters;

    // density distributions (D3Q19 model)
    CCL::CMem cMemDensityDistributions;

    // flags giving e.g. obstacles, (maybe gas) or other properties
    CCL::CMem cMemCellFlags;
    CCL::CMem cMemBC;

    // Velocity (3 components) and Density (1 component) of cell
    CCL::CMem cMemVelocity;
    CCL::CMem cMemDensity;

    size_t computation_kernel_count;
    size_t threads_per_dimension;
    size_t domain_cells_count;
    size_t domain_x, domain_y, domain_z;

    bool store_velocity; // store velocity for visualization
    bool store_density;
private:

    std::vector<std::string> split(const std::string& str,
            const std::string& delimiter = " ") {
        std::vector<std::string> tokens;

        std::string::size_type lastPos = 0;
        std::string::size_type pos = str.find(delimiter, lastPos);

        while (std::string::npos != pos) {
            // Found a token, add it to the vector.
            //std::cout << str.substr(lastPos, pos - lastPos) << std::endl;
            tokens.push_back(str.substr(lastPos, pos - lastPos));
            lastPos = pos + delimiter.size();
            pos = str.find(delimiter, lastPos);
        }

        tokens.push_back(str.substr(lastPos, str.size() - lastPos));
        return tokens;
    }


    inline void enqueueCopyRectKernel(CCL::CMem& src, CCL::CMem& dst,
            const int src_offset, CVector<3, int> src_origin,
            CVector<3, int> src_size, const int dst_offset,
            CVector<3, int> dst_origin, CVector<3, int> dst_size,
            CVector<3, int> block_size, bool withBarrier = true)
    {
/*
        enqueueCopyRectKernel(cBuffer, cMemCellFlags, 0,
                CVector<3, int>(0, 0, 0), size, 0, origin,
                this->domain_cells, size);
*/
        // set kernel args
        // cKernelCopyRect.setArgSize(20);      ///< reserve a vector argument container of 20 elements
        // source args
        cKernelCopyRect.setArg(0, src);
        cKernelCopyRect.setArg(1, src_offset);      printf("src_offset: %i\n", src_offset);
        cKernelCopyRect.setArg(2, src_origin[0]);   printf("src_origin[0]: %i\n", src_origin[0]);
        cKernelCopyRect.setArg(3, src_origin[1]);   printf("src_origin[1]: %i\n", src_origin[1]);
        cKernelCopyRect.setArg(4, src_origin[2]);   printf("src_origin[2]: %i\n", src_origin[2]);
        cKernelCopyRect.setArg(5, src_size[0]);     printf("src_size[0]: %i\n", src_size[0]);
        cKernelCopyRect.setArg(6, src_size[1]);     printf("src_size[1]: %i\n", src_size[1]);
        cKernelCopyRect.setArg(7, src_size[2]);     printf("src_size[2]: %i\n", src_size[2]);

        // destination args
        cKernelCopyRect.setArg(8, dst);
        cKernelCopyRect.setArg(9, dst_offset);      printf("dst_offset: %i\n", dst_offset);
        cKernelCopyRect.setArg(10, dst_origin[0]);  printf("dst_origin[0]: %i\n", dst_origin[0]);
        cKernelCopyRect.setArg(11, dst_origin[1]);  printf("dst_origin[1]: %i\n", dst_origin[1]);
        cKernelCopyRect.setArg(12, dst_origin[2]);  printf("dst_origin[2]: %i\n", dst_origin[2]);
        cKernelCopyRect.setArg(13, dst_size[0]);    printf("dst_size[0]: %i\n", dst_size[0]);
        cKernelCopyRect.setArg(14, dst_size[1]);    printf("dst_size[1]: %i\n", dst_size[1]);
        cKernelCopyRect.setArg(15, dst_size[2]);    printf("dst_size[2]: %i\n", dst_size[2]);
        cKernelCopyRect.setArg(16, block_size[0]);  printf("block_size[0]: %i\n", block_size[0]);

        size_t lGlobalSize[2];
        lGlobalSize[0] = block_size[1];     printf("lGlobalSize[0]: %lu\n", lGlobalSize[0]);
        lGlobalSize[1] = block_size[2];     printf("lGlobalSize[1]: %lu\n", lGlobalSize[1]);
        // enqueue the CopyRect kernel
        cCommandQueue.enqueueNDRangeKernel(cKernelCopyRect, // kernel
                2, // dimensions
                NULL,   //this->domain_cells,       // subdomain size in each directions
                this->domain_cells.elements(),
                NULL, // global work offset
                lGlobalSize,
                cKernelCopyRect_WorkGroupSize, //NULL,
                cKernelCopyRect.kernelArgsVec); 
        if (withBarrier)
            cCommandQueue.enqueueBarrier();
    }

public:
    size_t simulation_step_counter;
    CError error;

    std::list<int> lbm_opencl_number_of_work_items_list; ///< list with number of threads for each successively created kernel
    std::list<int> lbm_opencl_number_of_registers_list; ///< list with number of registers for each thread threads for each successively created kernel

    // viscosity parameters:
    // http://en.wikipedia.org/wiki/Viscosity
    //
    CLbmSolver(int UID, CCL::CCommandQueue &p_cCommandQueue,
            CCL::CContext &p_cContext, CCL::CDevice &p_cDevice, int BC[3][2],
            CDomain<T> &domain, CVector<3, T> &p_d_gravitation, T p_d_viscosity,
            size_t p_computation_kernel_count,
            size_t p_threads_per_dimension,
            //bool p_debug,
            bool p_store_velocity, bool p_store_density, T p_d_timestep,
            CVector<4, T>& _drivenCavityVelocity,
            std::list<int> &p_lbm_opencl_number_of_work_items_list, ///< list with number of threads for each successively created kernel
            std::list<int> &p_lbm_opencl_number_of_registers_list ///< list with number of registers for each thread threads for each successively created kernel
            ) :
            CLbmSkeleton<T>(CDomain<T>(domain), _drivenCavityVelocity), //drivenCavityVelocity(100.0, 0,0, 1),
            drivenCavityVelocity(_drivenCavityVelocity), _UID(UID), cCommandQueue(
                    p_cCommandQueue), cContext(p_cContext), cDevice(p_cDevice), cDeviceInfo(
                    p_cDevice), computation_kernel_count(p_computation_kernel_count),
                    threads_per_dimension(p_threads_per_dimension), lbm_opencl_number_of_work_items_list(
                    p_lbm_opencl_number_of_work_items_list), lbm_opencl_number_of_registers_list(
                    p_lbm_opencl_number_of_registers_list) {
        CLbmSkeleton<T>::init(p_d_gravitation, p_d_viscosity, 1.0);

        printf(" ->CLbmSolver::CLbmSolver() ");

        if (CLbmSkeleton<T>::error())
            error << CLbmSkeleton<T>::error.getString();

        store_velocity = p_store_velocity;
        store_density = p_store_density;
        //debug = p_debug;

#if DEBUG
        std::cout << "CL_VERSION " << _cl_version << std::endl;
#endif
        // setting the boundary conditions
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 2; j++)
                _BC[i][j] = BC[i][j];

        // printf("CLbmSolver._UID = %i\n", _UID);

        reload();
    }

    void addDrivenCavityValue(T value) {
        drivenCavityVelocity[0] += value;
        CVector<4, T> paramDrivenCavityVelocity = drivenCavityVelocity;
        paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

        cKernelInit.setArg(4, paramDrivenCavityVelocity[0]);
        cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);
        cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);
    }

    void reload() {
        // printf(" ->CLbmSolver::reload()");
        /**
         * WORK GROUP SIZE
         *
         * initialize the variable (postfixed appropriately with _WorkGroupSize and _MaxRegisters)
         * with either the standard value (max_local_work_group_size) or with the value from the list
         */
#define INIT_WORK_GROUP_SIZE(variable)                                      \
        if (it != this->lbm_opencl_number_of_work_items_list.end())     \
        {                                                               \
            /*variable##_WorkGroupSize = (*it != 0 ? *it : computation_kernel_count);*/ \
            variable##_WorkGroupSize = (*it != 0 ? *it : threads_per_dimension); \
            it++;                                                       \
        }                                                               \
        else                                                            \
        {                                                               \
            /*variable##_WorkGroupSize = computation_kernel_count;*/    \
            variable##_WorkGroupSize = threads_per_dimension;           \
        }                                                               \
                                                                        \
        /* enlarge global work group size to be a multiple of collprop_work_group_size */   \
        if (domain_cells_count % variable##_WorkGroupSize != 0)         \
            variable##_GlobalWorkGroupSize = (domain_cells_count / variable##_WorkGroupSize + 1) * variable##_WorkGroupSize;    \
                                                                        \
        if (ir != this->lbm_opencl_number_of_registers_list.end())      \
        {                                                               \
            variable##_MaxRegisters = (*ir != 0 ? *ir : 0);             \
            ir++;                                                       \
        }                                                               \
        else                                                            \
        {                                                               \
            variable##_MaxRegisters  = 0;                               \
        }

        std::list<int>::iterator it =
                this->lbm_opencl_number_of_work_items_list.begin();
        std::list<int>::iterator ir =
                this->lbm_opencl_number_of_registers_list.begin();

        INIT_WORK_GROUP_SIZE(cKernelInit);
        INIT_WORK_GROUP_SIZE(cLbmKernelAlpha);
        INIT_WORK_GROUP_SIZE(cLbmKernelBeta);
        INIT_WORK_GROUP_SIZE(cKernelCopyRect);

        /**
         * program defines
         */
        domain_cells_count = this->domain_cells.elements();

        domain_x = this->domain_cells.data[0];
        domain_y = this->domain_cells.data[1];
        domain_z = this->domain_cells.data[2];

        std::ostringstream cuda_program_defines;

        cuda_program_defines << " -D PROBLEM_DEFAULTS=0";
        cuda_program_defines << " -D SIZE_DD_HOST_BYTES=" << SIZE_DD_HOST_BYTES;
        cuda_program_defines << " -D DOMAIN_CELLS_X=" << this->domain_cells[0];
        cuda_program_defines << " -D DOMAIN_CELLS_Y=" << this->domain_cells[1];
        cuda_program_defines << " -D DOMAIN_CELLS_Z=" << this->domain_cells[2];
        cuda_program_defines << " -D GLOBAL_WORK_GROUP_SIZE="
            << (this->domain_cells[0] * this->domain_cells[1] * this->domain_cells[2]);
        cuda_program_defines << " -D FLAG_OBSTACLE=" << LBM_FLAG_OBSTACLE;
        cuda_program_defines << " -D FLAG_FLUID=" << LBM_FLAG_FLUID;
        cuda_program_defines << " -D FLAG_VELOCITY_INJECTION=" << LBM_FLAG_VELOCITY_INJECTION;
        cuda_program_defines << " -D FLAG_GHOST_LAYER=" << LBM_FLAG_GHOST_LAYER;

        if (store_velocity)
            cuda_program_defines << " -D STORE_VELOCITY=1";

        if (store_density)
            cuda_program_defines << " -D STORE_DENSITY=1";

        if (typeid(T) == typeid(float))
        {
            cuda_program_defines << " -D TYPE_FLOAT=1";
        }
        else if (typeid(T) == typeid(double))
        {
            cuda_program_defines << " -D TYPE_DOUBLE=1";
        }
        else
        {
            error << "unsupported class type T" << std::endl;
            return;
        }

#if DEBUG
        std::cout << cuda_program_defines.str() << std::endl;
#endif
        /*
         * ALLOCATE BUFFERS
         */
         cMemDensityDistributions.create(cContext,
                sizeof(T) * domain_cells_count * SIZE_DD_HOST);
        cMemCellFlags.create(cContext,
                sizeof(int) * domain_cells_count);
        cMemVelocity.create(cContext,
                sizeof(T) * domain_cells_count * 3);
        cMemDensity.create(cContext, 
                sizeof(T) * domain_cells_count);

        int bc_linear[6];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 2; j++)
                bc_linear[i * 2 + j] = _BC[i][j];

#if DEBUG
        std::cout << "BOUNDARY CONDITION: "<< std::endl;
        for(int i = 0; i < 6; i++)
        std::cout << " " << bc_linear[i];
        std::cout << std::endl;
#endif
    
        cMemBC.createCopyToDevice(cContext, sizeof(int) * 6, bc_linear);

        /**
         * prepare INIT KERNEL DATA
         * nothing to do because we use COLLISION KERNEL DATA
         */
        cDeviceInfo.loadDeviceInfo(cDevice);    ///< get GPU specifications. CC needed for kernel compilation

        std::string cProgramCompileOptionsString;
        std::stringstream cProgramDefinesPostfixString;
        std::stringstream cProgramCompileOptionsStream; ///< string to append
        std::string initKernelModuleFileName = "src/cu_programs/lbm_init";  //cu_programs/
        char charbuf[255];

        /*
         * INIT kernel
         */
        cProgramDefinesPostfixString << " -D LOCAL_WORK_GROUP_SIZE=";
        cProgramDefinesPostfixString << cKernelInit_WorkGroupSize;

        cProgramCompileOptionsStream << "nvcc -x cu -keep -keep-dir src/cu_programs/ "; // -keep retains the PTX and all other intermediate compile files
        if (cKernelInit_MaxRegisters != 0)
        {
            cProgramCompileOptionsStream << "-maxrregcount ";
            cProgramCompileOptionsStream << cKernelInit_MaxRegisters;
        }

        cProgramCompileOptionsStream << cuda_program_defines.str() + cProgramDefinesPostfixString.str();
        // cProgramCompileOptionsStream << " -arch=sm_";
        // cProgramCompileOptionsStream << cDeviceInfo.execution_capabilities;      ///< GPU compute capability
        cProgramCompileOptionsStream << " -arch=sm_3";
        cProgramCompileOptionsStream << "0 -m64 -I. -dc ";
        cProgramCompileOptionsStream << initKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".cu -o ";
        cProgramCompileOptionsStream << initKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".o";

        cProgramCompileOptionsString = cProgramCompileOptionsStream.str();

        cKernelInit_GlobalWorkGroupSize = domain_cells_count;
        if (cKernelInit_GlobalWorkGroupSize % cKernelInit_WorkGroupSize != 0)
            cKernelInit_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize
                    / cKernelInit_WorkGroupSize + 1)
                    * cKernelInit_WorkGroupSize;
        CCL::CProgram cProgramInit;

        // Compile the CUDA kernels using nvcc
        if (simulation_step_counter == 0)
        {
         cProgramInit.executeCommand(cProgramCompileOptionsString.c_str(), initKernelModuleFileName.c_str());
        }
        // cProgramInit.executeCommand(cProgramCompileOptionsString.c_str(), initKernelModuleFileName.c_str());
        // MPI_Barrier(MPI_COMM_WORLD);

        // Load init_lbm.ptx file with a CUDA module(program)
        initKernelModuleFileName = initKernelModuleFileName + ".ptx";
        cProgramInit.load(cContext, initKernelModuleFileName.c_str());

        if (cProgramInit.error()) {
            error << "failed to compile lbm_init.cl" << CError::endl;
            error << cProgramInit.error.getString() << CError::endl;
            return;
        }
#if DEBUG
// #if 1
        std::cout << "KernelInit:   local_work_group_size: " << cKernelInit_WorkGroupSize << "      max_registers: " << cKernelInit_MaxRegisters << std::endl;
#endif

        /*
         * ALPHA
         */
        std::string alphaKernelModuleFileName = "src/cu_programs/lbm_alpha";
        cProgramDefinesPostfixString.str(std::string());    // reset the contents of the StringStream
        cProgramDefinesPostfixString << " -D LOCAL_WORK_GROUP_SIZE=";
        cProgramDefinesPostfixString << cLbmKernelAlpha_WorkGroupSize;

        // reset the contents of the StringStream
        cProgramCompileOptionsStream.str(std::string());

        cProgramCompileOptionsStream << "nvcc -x cu -keep -keep-dir src/cu_programs/ "; // -keep retains the PTX and all other intermediate compile files
        if (cLbmKernelAlpha_MaxRegisters != 0)
        {
            cProgramCompileOptionsStream << "-maxrregcount ";
            cProgramCompileOptionsStream << cLbmKernelAlpha_MaxRegisters;
        }

        cLbmKernelAlpha_GlobalWorkGroupSize = domain_cells_count;
        if (cLbmKernelAlpha_GlobalWorkGroupSize % cLbmKernelAlpha_WorkGroupSize
                != 0)
            cLbmKernelAlpha_GlobalWorkGroupSize =
                    (cKernelInit_GlobalWorkGroupSize
                            / cLbmKernelAlpha_WorkGroupSize + 1)
                            * cLbmKernelAlpha_WorkGroupSize;

        cProgramCompileOptionsStream << cuda_program_defines.str() + cProgramDefinesPostfixString.str();
        // cProgramCompileOptionsStream << " -arch=sm_";
        // cProgramCompileOptionsStream << cDeviceInfo.execution_capabilities;      ///< GPU compute capability
        cProgramCompileOptionsStream << " -arch=sm_3";
        cProgramCompileOptionsStream << "0 -m64 -I. -dc ";
        cProgramCompileOptionsStream << alphaKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".cu -o ";
        cProgramCompileOptionsStream << alphaKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".o";

        cProgramCompileOptionsString = cProgramCompileOptionsStream.str();

        CCL::CProgram cProgramAlpha;
        
        // Compile cuda kernels to ptx file using nvcc
        if (simulation_step_counter == 0)
        {
         cProgramAlpha.executeCommand(cProgramCompileOptionsString.c_str(), alphaKernelModuleFileName.c_str());
        }
        // cProgramAlpha.executeCommand(cProgramCompileOptionsString.c_str(), alphaKernelModuleFileName.c_str());
        // MPI_Barrier(MPI_COMM_WORLD);

        //TODO: check which program phase I need to load the module
        alphaKernelModuleFileName = alphaKernelModuleFileName + ".ptx";
        cProgramAlpha.load(cContext, alphaKernelModuleFileName.c_str());

// TODO: CProgram.build() member does not exist! Make sure the CUDA version does not need it.
//      cProgramAlpha.build(cDevice, cProgramCompileOptionsString.c_str());
        if (cProgramAlpha.error()) {
            error << "failed to compile lbm_alpha.cl" << CError::endl;
            error << cProgramAlpha.error.getString() << CError::endl;
            return;
        }
#if DEBUG
// #if 1
        std::cout << "KernelAlpha:  local_work_group_size: " << cLbmKernelAlpha_WorkGroupSize << "      max_registers: " << cLbmKernelAlpha_MaxRegisters << std::endl;
#endif
        /*
         * BETA
         */
        std::string betaKernelModuleFileName = "src/cu_programs/lbm_beta";
        cProgramDefinesPostfixString.str(std::string());    // reset the contents of the StringStream
        cProgramDefinesPostfixString << " -D LOCAL_WORK_GROUP_SIZE=";
        cProgramDefinesPostfixString << cLbmKernelBeta_WorkGroupSize;

        // reset the contents of the StringStream
        cProgramCompileOptionsStream.str(std::string());

        cProgramCompileOptionsStream << "nvcc -x cu -keep -keep-dir src/cu_programs/ "; // -keep retains the PTX and all other intermediate compile files
        if (cLbmKernelBeta_MaxRegisters != 0)
        {
            cProgramCompileOptionsStream << "-maxrregcount ";
            cProgramCompileOptionsStream << cLbmKernelBeta_MaxRegisters;
        }

        cLbmKernelBeta_GlobalWorkGroupSize = domain_cells_count;
        if (cLbmKernelBeta_GlobalWorkGroupSize % cLbmKernelBeta_WorkGroupSize
                != 0)
            cLbmKernelBeta_GlobalWorkGroupSize =
                    (cKernelInit_GlobalWorkGroupSize
                            / cLbmKernelBeta_WorkGroupSize + 1)
                            * cLbmKernelBeta_WorkGroupSize;

        cProgramCompileOptionsStream << cuda_program_defines.str() + cProgramDefinesPostfixString.str();
        // cProgramCompileOptionsStream << " -arch=sm_";
        // cProgramCompileOptionsStream << cDeviceInfo.execution_capabilities;      ///< GPU compute capability
        cProgramCompileOptionsStream << " -arch=sm_3";
        cProgramCompileOptionsStream << "0 -m64 -I. -dc ";
        cProgramCompileOptionsStream << betaKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".cu -o ";
        cProgramCompileOptionsStream << betaKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".o";

        cProgramCompileOptionsString = cProgramCompileOptionsStream.str();

        CCL::CProgram cProgramBeta;
        
        // compile cuda kernel to ptx
        if (simulation_step_counter == 0)
        {
         cProgramBeta.executeCommand(cProgramCompileOptionsString.c_str(), betaKernelModuleFileName.c_str());
        }
        // cProgramBeta.executeCommand(cProgramCompileOptionsString.c_str(), betaKernelModuleFileName.c_str());
        // MPI_Barrier(MPI_COMM_WORLD);

        //TODO: check which program phase I need to load the module
        betaKernelModuleFileName = betaKernelModuleFileName + ".ptx";
        cProgramBeta.load(cContext, betaKernelModuleFileName.c_str());

// TODO: CProgram.build() member does not exist! Make sure the CUDA version does not need it.
//      cProgramBeta.build(cDevice, cProgramCompileOptionsString.c_str());
        if (cProgramBeta.error()) {
            error << "failed to compile lbm_beta.cl" << CError::endl;
            error << cProgramBeta.error.getString() << CError::endl;

            return;
        }
#if DEBUG
// #if 1
        std::cout << "KernelBeta:   local_work_group_size: " << cLbmKernelBeta_WorkGroupSize << "       max_registers: " << cLbmKernelBeta_MaxRegisters << std::endl;
#endif

        /*
         * INIT CopyBufferRect
         */
        std::string copyRectKernelModuleFileName = "src/cu_programs/copy_buffer_rect";
        cProgramDefinesPostfixString.str(std::string());    // reset the contents of the StringStream
        cProgramDefinesPostfixString << " -D LOCAL_WORK_GROUP_SIZE=";
        cProgramDefinesPostfixString << cKernelCopyRect_WorkGroupSize;

        // reset the contents of the StringStream
        cProgramCompileOptionsStream.str(std::string());

        cProgramCompileOptionsStream << "nvcc -x cu -keep -keep-dir src/cu_programs/ "; // -keep retains the PTX and all other intermediate compile files
        if (cKernelCopyRect_MaxRegisters != 0)
        {
            cProgramCompileOptionsStream << "-maxrregcount ";
            cProgramCompileOptionsStream << cKernelCopyRect_MaxRegisters;
        }


        cKernelCopyRect_GlobalWorkGroupSize = domain_cells_count;
        if (cKernelInit_GlobalWorkGroupSize % cKernelInit_WorkGroupSize != 0)
            cKernelInit_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize
                    / cKernelInit_WorkGroupSize + 1)
                    * cKernelInit_WorkGroupSize;

        cProgramCompileOptionsStream << cuda_program_defines.str() + cProgramDefinesPostfixString.str();
        // cProgramCompileOptionsStream << " -arch=sm_";
        // cProgramCompileOptionsStream << cDeviceInfo.execution_capabilities;      ///< GPU compute capability
        cProgramCompileOptionsStream << " -arch=sm_3";
        cProgramCompileOptionsStream << "0 -m64 -I. -dc ";
        cProgramCompileOptionsStream << copyRectKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".cu -o ";
        cProgramCompileOptionsStream << copyRectKernelModuleFileName.data();
        cProgramCompileOptionsStream << ".o";

        cProgramCompileOptionsString = cProgramCompileOptionsStream.str();

        CCL::CProgram cProgramCopyRect;
        
        // compile kernel with nvcc
        if (simulation_step_counter == 0)
        {
         cProgramCopyRect.executeCommand(cProgramCompileOptionsString.c_str(), copyRectKernelModuleFileName.c_str());
        }
        // cProgramCopyRect.executeCommand(cProgramCompileOptionsString.c_str(), copyRectKernelModuleFileName.c_str());
        // MPI_Barrier(MPI_COMM_WORLD);

        //TODO: check which program phase I need to load the module
        copyRectKernelModuleFileName = copyRectKernelModuleFileName + ".ptx";
        cProgramCopyRect.load(cContext, copyRectKernelModuleFileName.c_str());

        
// TODO: CProgram.build() member does not exist! Make sure the CUDA version does not need it.
        //cProgramCopyRect.build(cDevice, cProgramCompileOptionsString.c_str());
        if (cProgramCopyRect.error()) {
            error << "failed to compile copy_buffer_rect.cl" << CError::endl;
            error << cProgramCopyRect.error.getString() << CError::endl;
            return;
        }
#if DEBUG
// #if 1
        std::cout << "KernelCopyRect:   local_work_group_size: " << cKernelCopyRect_WorkGroupSize << "      max_registers: " << cKernelCopyRect_MaxRegisters << std::endl;
#endif

        /**
         * create kernels and setup arguments
         */
        CVector<4, T> paramDrivenCavityVelocity = CLbmSkeleton<T>::drivenCavityVelocity;
//      CVector<4, T> paramDrivenCavityVelocity = drivenCavityVelocity;
//      paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

#if DEBUG
        {
            std::cout << "driven cavity velocity: " << paramDrivenCavityVelocity << std::endl;
            std::cout << "inverse tau: " << this->inv_tau << std::endl;
            std::cout << "d_timestep: " << this->d_timestep << std::endl;
            std::cout << "gravitaton: " << this->gravitation << std::endl;
        }
#endif
        /**
         * SETUP ARGUMENTS
         */
        // initialization kernel
        cKernelInit.create(cProgramInit, "init_kernel"); // attach kernel to CUDA module
        cKernelInit.setArg(0, cMemDensityDistributions);
        cKernelInit.setArg(1, cMemCellFlags);
        cKernelInit.setArg(2, cMemVelocity);
        cKernelInit.setArg(3, cMemDensity);
        cKernelInit.setArg(4, cMemBC);
        cKernelInit.setArg(5, paramDrivenCavityVelocity[0]);

        // collision and propagation kernels (alpha and beta)
        cLbmKernelAlpha.create(cProgramAlpha, "lbm_kernel_alpha");   // attach kernel to CUDA module
        cLbmKernelAlpha.setArg(0, cMemDensityDistributions);
        cLbmKernelAlpha.setArg(1, cMemCellFlags);
        cLbmKernelAlpha.setArg(2, cMemVelocity);
        cLbmKernelAlpha.setArg(3, cMemDensity);
        cLbmKernelAlpha.setArg(4, this->inv_tau);
        cLbmKernelAlpha.setArg(5, this->gravitation[0]);
        cLbmKernelAlpha.setArg(6, this->gravitation[1]);
        cLbmKernelAlpha.setArg(7, this->gravitation[2]);
        cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);

        cLbmKernelBeta.create(cProgramBeta, "lbm_kernel_beta"); // attach kernel to CUDA module
        cLbmKernelBeta.setArg(0, cMemDensityDistributions);
        cLbmKernelBeta.setArg(1, cMemCellFlags);
        cLbmKernelBeta.setArg(2, cMemVelocity);
        cLbmKernelBeta.setArg(3, cMemDensity);
        cLbmKernelBeta.setArg(4, this->inv_tau);
        cLbmKernelBeta.setArg(5, this->gravitation[0]);
        cLbmKernelBeta.setArg(6, this->gravitation[1]);
        cLbmKernelBeta.setArg(7, this->gravitation[2]);
        cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);

        cKernelCopyRect.create(cProgramCopyRect, "copy_buffer_rect"); // attach kernel to CUDA module

        printf("\n");

        reset();
    }

    void reset() {
        // printf(" ->CLbmSolver::reset()");

        simulation_step_counter = 0;

#if DEBUG
        std::cout << "Init Simulation: " << std::flush;
#endif
        // TODO: decide how grid_size_x value is set when calling enqueueNDRangeKernel
        // printf("CLbmSolver.reset() before call\n");
        cCommandQueue.enqueueNDRangeKernel(cKernelInit, ///< kernel
                1,                                      ///< dimensions
                NULL,   //this->domain_cells.size(),
                domain_cells_count,                     ///< total elements to process 
                NULL,                                   ///< global work offset
                &cKernelInit_GlobalWorkGroupSize,
                cKernelInit_WorkGroupSize,
                cKernelInit.kernelArgsVec);

        cCommandQueue.enqueueBarrier();


        // int *cell_flags = new int[domain_cells_count];
        // storeFlags(cell_flags);

        // printf("FLAG_OBSTACLE: %lu\n", (1 << 0));
        // printf("FLAG_FLUID: %lu\n", (1 << 1));
        // printf("FLAG_GHOST_LAYER: %lu\n", (1 << 3));
        // printf("cell flags: size= %lu\n", domain_cells_count);

        // for (size_t i = 0; i < domain_cells_count; i++)
        // {
        //     printf(" %i ", cell_flags[i]);
        //     if ((i + 1) % (domain_x) == 0)
        //     {
        //         printf("\n");
        //     }
        // }
        // printf("\n");

        // delete [] cell_flags;

#if DEBUG
        std::cout << "OK" << std::endl;
#endif
    }


    void simulationStepAlpha() {
        // printf(" -> CLbmSolver::simulationStepAlpha()");
#if DEBUG
        std::cout << "--> Running Alpha kernel" << std::endl;
#endif
        cCommandQueue.enqueueNDRangeKernel(
                cLbmKernelAlpha,    // kernel
                1,                  // dimensions
                NULL,   //this->domain_cells.size(),
                domain_cells_count, // total number of elements
                NULL,               // global work offset
                &cLbmKernelAlpha_GlobalWorkGroupSize,
                cLbmKernelAlpha_WorkGroupSize,
                cLbmKernelAlpha.kernelArgsVec);
    }

    void simulationStepBeta() {
        // printf(" -> CLbmSolver::simulationStepBeta()");
#if DEBUG
        std::cout << "--> Running BETA kernel" << std::endl;
#endif
        cCommandQueue.enqueueNDRangeKernel(
                cLbmKernelBeta, // kernel
                1, // dimensions
                NULL,   //this->domain_cells.size(),
                domain_cells_count,
                NULL, // global work offset
                &cLbmKernelBeta_GlobalWorkGroupSize,
                cLbmKernelBeta_WorkGroupSize,
                cLbmKernelBeta.kernelArgsVec);

        // print density values which are overwritten with FLAGS
        // printDensityMemObj(size_t elems)
    }

    /**
     * start one simulation step (enqueue kernels)
     */
    void simulationStep() {
        /*
         * collision kernels are inserted as 1d kernels because they work cell wise without neighboring information
         */
        // printf(" -> CLbmSolver::simulationStep()");
        if (simulation_step_counter & 1) {
            simulationStepAlpha();
        } else {
            simulationStepBeta();
        }
        cCommandQueue.enqueueBarrier();

        simulation_step_counter++;
    }

    /**
     * wait until all kernels have finished
     */
    void wait() {
        cCommandQueue.finish();
    }

    /**
     * store density distribution values to host memory
     */
    void storeDensityDistribution(T *dst) {
        // printf(" -> CLbmSolver::storeDensityDistribution(1)");
        
        size_t byte_size = cMemDensityDistributions.getSize();

        cCommandQueue.enqueueReadBuffer(cMemDensityDistributions, CU_TRUE, // sync reading
                0, byte_size, dst);
    }

    void storeDensityDistribution(T *dst, CVector<3, int> &origin,
            CVector<3, int> &size) {
        // printf(" -> CLbmSolver::storeDensityDistribution(3)");
        
        CCL::CMem cBuffer;
        //cBuffer.create(cContext, CL_MEM_READ_WRITE,
        //      sizeof(T) * size.elements() * SIZE_DD_HOST, NULL);
        cBuffer.create(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST);

        for (int f = 0; f < SIZE_DD_HOST; f++)
        {
            enqueueCopyRectKernel(cMemDensityDistributions, cBuffer,
                    f * this->domain_cells_count, origin,
                    this->domain_cells, f * size.elements(),
                    CVector<3, int>(0, 0, 0), size, size, false);
        }
        cCommandQueue.enqueueBarrier();

        //clEnqueueBarrier(cCommandQueue.command_queue);    This sync point not needed, same thing done in last call above
        cCommandQueue.enqueueReadBuffer(cBuffer, CU_TRUE, // sync reading
                0, cBuffer.getSize(), dst);
    }

    void setDensityDistribution(T *src, CVector<3, int> &origin,
            CVector<3, int> &size) {
        // printf(" -> CLbmSolver::setDensityDistribution(3)");
        
        CCL::CMem cBuffer;
        cBuffer.createCopyToDevice(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST, src);

        for (int f = 0; f < SIZE_DD_HOST; f++)
        {
            enqueueCopyRectKernel(cBuffer, cMemDensityDistributions,
                    f * size.elements(), CVector<3, int>(0, 0, 0), size,
                    f * this->domain_cells_count, origin,
                    this->domain_cells, size, false);
        }
        cCommandQueue.enqueueBarrier();

        // inline void enqueueWriteBuffer(  CMem &cMem, bool block_write, size_t offset, size_t buffer_size, const void *buffer_ptr)
    }

    void setDensityDistribution(T *src, CVector<3, int> &origin,
            CVector<3, int> &size, CVector<3, int> norm) {
        // printf(" -> CLbmSolver::setDensityDistribution(4)");
        CCL::CMem cBuffer;
        //cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        //      sizeof(T) * size.elements() * SIZE_DD_HOST, src);
        cBuffer.createCopyToDevice(cContext, 
                sizeof(T) * size.elements() * SIZE_DD_HOST, src);


        for (int f = 0; f < SIZE_DD_HOST; f++)
        {
            if (norm.dotProd(lbm_units[f]) > 0)
            {
                enqueueCopyRectKernel(cBuffer, cMemDensityDistributions,
                        f * size.elements(), CVector<3, int>(0, 0, 0), size,
                        f * this->domain_cells_count, origin,
                        this->domain_cells, size, false);
            }
        }
        cCommandQueue.enqueueBarrier();
    }
    /**
     * store velocity and density values to host memory
     * this is useful for a host memory based visualization
     *
     * @param dst The buffer that will contain the return values
     */
    void storeVelocity(T *dst) {
        // printf(" -> CLbmSolver::storeVelocity(1)");
        size_t byte_size = cMemVelocity.getSize();

        // printf("byte_size: %lu \n", byte_size);
        cCommandQueue.enqueueReadBuffer(cMemVelocity, CU_TRUE, // sync reading
             0, byte_size, dst);

        // cMemVelocity.printVelMemObj(byte_size/sizeof(T));
    }

    /**
     * Store a block of velocity data from device to host
     *
     * @param dst The buffer that will contain the return values
     * @param origin The origin point of data block
     * @param size The size of data block
     */
    void storeVelocity(T* dst, CVector<3, int> &origin, CVector<3, int> &size) {
        // printf(" -> CLbmSolver::storeVelocity(3)");
        CCL::CMem cBuffer;
        //cBuffer.create(cContext, CL_MEM_READ_WRITE,
        //      sizeof(T) * size.elements() * 3, NULL);
        cBuffer.create(cContext, 
                sizeof(T) * size.elements() * 3);

        for (int dim = 0; dim < 3; dim++) {
            enqueueCopyRectKernel(cMemVelocity, cBuffer,
                    dim * this->domain_cells_count, origin,
                    this->domain_cells, dim * size.elements(),
                    CVector<3, int>(0, 0, 0), size, size, false);
        }
        cCommandQueue.enqueueBarrier();

        cCommandQueue.enqueueReadBuffer(cBuffer, CU_TRUE, // sync reading
                0, cBuffer.getSize(), dst);
    }

    /**
     * Set a block of velocity data from host to device
     *
     * @param src The buffer that contain the input values
     * @param origin The origin point of data block
     * @param size The size of data block
     */
    void setVelocity(T* src, CVector<3, int> &origin, CVector<3, int> &size) {
        CCL::CMem cBuffer;
        //cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        //      sizeof(T) * size.elements() * 3, src);
        cBuffer.createCopyToDevice(cContext, 
                sizeof(T) * size.elements() * 3, src);

        // printf(" -> CLbmSolver::setVelocity(3)");
        for (int dim = 0; dim < 3; dim++) {
            enqueueCopyRectKernel(cBuffer, cMemVelocity,
                    dim * size.elements(), CVector<3, int>(0, 0, 0), size,
                    dim * this->domain_cells_count, origin,
                    this->domain_cells, size, false);
        }
        cCommandQueue.enqueueBarrier();
    }
    /**
     * Store velocity data from device to host
     *
     * @param dst The buffer that will contain the return values
     */
    void storeDensity(T *dst) {
        size_t byte_size = cMemDensity.getSize();

        // printf(" -> CLbmSolver::storeDensity(1)");
        cCommandQueue.enqueueReadBuffer(cMemDensity, CU_TRUE, // sync reading
                0, byte_size, dst);

/*        printf("density as FLAGS: size= %i\n", byte_size/sizeof(T) );
        for (size_t i = 0; i < byte_size/sizeof(T); i++)
        {
            printf(" %i ", (int)dst[i]);
            if ((i + 1) % (domain_x) == 0)
            {
                printf("\n");
            }
        }
        printf("\n");
*/
    }

    /**
     * Store a block of density data from device to host
     *
     * @param dst The buffer that will contain the return values
     * @param origin The origin point of data block
     * @param size The size of data block
     */
    void storeDensity(T *dst, CVector<3, int> &origin, CVector<3, int> &size) {
        CCL::CMem cBuffer;
        //cBuffer.create(cContext, CL_MEM_READ_WRITE, sizeof(T) * size.elements(),
        //      NULL);
        cBuffer.create(cContext, sizeof(T) * size.elements());

        // printf(" -> CLbmSolver::storeDensity(3)");
        enqueueCopyRectKernel(cMemDensity, cBuffer, 0, origin,
                this->domain_cells, 0, CVector<3, int>(0, 0, 0), size,
                size);

        cCommandQueue.enqueueReadBuffer(cBuffer, CU_TRUE, // sync reading
                0, cBuffer.getSize(), dst);
    }

    /**
     * Set a block of density data from host to device
     *
     * @param src The buffer that will contain the source values
     * @param origin The origin point of data block
     * @param size The size of data block
     */
    void setDensity(T *src, CVector<3, int> &origin, CVector<3, int> &size)
    {
        CCL::CMem cBuffer;
        cBuffer.createCopyToDevice(cContext,
                sizeof(T) * size.elements(), src);

        // printf(" -> CLbmSolver::setDensity(4)");
        enqueueCopyRectKernel(cBuffer, cMemDensity, 0,
                CVector<3, int>(0, 0, 0), size, 0, origin,
                this->domain_cells, size);
    }

    void storeFlags(int *dst) {
        size_t byte_size = cMemDensity.getSize();


        cCommandQueue.enqueueReadBuffer(cMemCellFlags, CU_TRUE, // sync reading
                0, byte_size, dst);
    }

    void storeFlags(int *dst, CVector<3, int> &origin, CVector<3, int> &size)
    {
        CCL::CMem cBuffer;
        cBuffer.create(cContext, sizeof(int) * size.elements());

        // For CUDA: using cuMemcpyDtoH instead of a CUStream, for now.
        cCommandQueue.enqueueReadBuffer(cBuffer, CU_TRUE, // sync reading
                0, cBuffer.getSize(), dst);
    }

    void setFlags(int *src, CVector<3, int> &origin, CVector<3, int> &size)
    {
        CCL::CMem cBuffer;
        cBuffer.createCopyToDevice(cContext, 
                sizeof(int) * size.elements(), src);

        printf("src [cBuffer -> cMemCellFlags]: size= %i\n", size.elements() );
        for (size_t i = 0; i < size.elements(); i++)
        {
            printf(" %i ", src[i]);
            if ((i + 1) % (domain_x) == 0)
            {
                printf("\n");
            }
        }
        printf("\n");

        enqueueCopyRectKernel(cBuffer, cMemCellFlags, 0,
                CVector<3, int>(0, 0, 0), size, 0, origin,
                this->domain_cells, size);
    }

private:
    void debugChar(CCL::CMem &cMem, size_t wrap_size = 20) {
        size_t char_size = cMem.getSize();
        char *buffer = new char[char_size];

        cCommandQueue.enqueueReadBuffer(cMem, CU_TRUE, // sync reading
                0, char_size, buffer);

        for (size_t i = 0; i < char_size; i++) {
            if (i % wrap_size == 0) {
                std::cout << std::endl;
                std::cout << (i / wrap_size) << ": ";
            }

            std::cout << (int) buffer[i] << " ";
        }

        delete[] buffer;
    }

    void debugFloat(CCL::CMem &cMem, size_t wrap_size = 20) {
        size_t byte_size = cMem.getSize();
        size_t T_size = byte_size / sizeof(T);

        T *buffer = new T[T_size];

        cCommandQueue.enqueueReadBuffer(cMem, CU_TRUE, // sync reading
                0, byte_size, buffer);

        std::streamsize ss = std::cout.precision();
        std::cout.precision(4);
        std::cout.setf(std::ios::fixed, std::ios::floatfield);

        for (size_t i = 0; i < T_size; i++) {
            if (i % wrap_size == 0) {
                std::cout << std::endl;
                std::cout << (i / wrap_size) << ": ";
            }

            std::cout << buffer[i] << " ";
        }

        std::cout.precision(ss);
        std::cout << std::resetiosflags(std::ios::fixed);

        delete[] buffer;
    }

    /**
     * debug handler to output some useful information
     */
public:
    void debug_print() {
        // read out DensityDistributions
        std::cout << "DENSITY DISTRIBUTIONS:";
        //debugFloat(cMemDensityDistributions, SIZE_DD_HOST);
        debugFloat(cMemDensityDistributions, 16);
        std::cout << std::endl;

        // read out Velocity
        std::cout << std::endl;
        std::cout << "VELOCITY:";
        debugFloat(cMemVelocity, 4 * 3);
        std::cout << std::endl;

        // read out VelocityDensity
        std::cout << std::endl;
        std::cout << "DENSITY:";
        debugFloat(cMemDensity, 4);
        std::cout << std::endl;

        // read out Flags
        std::cout << std::endl;
        std::cout << "FLAGS:";
        debugChar(cMemCellFlags, 4 * 4);
        std::cout << std::endl;
    }

    /**
     * debug density distributions (works only in linear mode!)
     * \param   dd_id   specifies the dd number
     */
    void debugDD(size_t dd_id = 0, size_t wrap_size = 16,
            size_t empty_line = 16) {
        CCL::CMem &cMem = cMemDensityDistributions;

        size_t byte_size = cMem.getSize();
        size_t T_size = byte_size / sizeof(T);

        T *buffer = new T[T_size];

        cCommandQueue.enqueueReadBuffer(cMem, CU_TRUE, // sync reading
                0, byte_size, buffer);

        std::streamsize ss = std::cout.precision();
        std::cout.precision(4);
        std::cout.setf(std::ios::fixed, std::ios::floatfield);

        size_t start = this->domain_cells.elements() * dd_id;
        size_t end = this->domain_cells.elements() * (dd_id + 1);

        for (size_t i = start; i < end; i++) {
            if (empty_line != wrap_size)
                if (i % empty_line == 0 && i != start) {
                    std::cout << std::endl;
                }

            if (i % wrap_size == 0) {
                if (i != start)
                    std::cout << std::endl;
                std::cout << (i / wrap_size) << ": ";
            }

            std::cout << buffer[i] << " ";
        }

        std::cout.precision(ss);
        std::cout << std::resetiosflags(std::ios::fixed);
        std::cout << std::endl;

        delete[] buffer;
    }

    float getVelocityChecksum() {
        T *velocity = new T[this->domain_cells.elements() * 3];
        int *flags = new int[this->domain_cells.elements() * 3];

        storeVelocity(velocity);
        storeFlags(flags);

        T *velx = velocity;
        T *vely = velocity + this->domain_cells.elements();
        T *velz = velocity + this->domain_cells.elements() * 2;

        float checksum = 0;
        for (int a = 0; a < this->domain_cells.elements(); a++)
            if (flags[a] == LBM_FLAG_FLUID)
                checksum += velx[a] + vely[a] + velz[a];

        delete[] flags;
        delete[] velocity;

        return checksum;
    }
};

#endif
