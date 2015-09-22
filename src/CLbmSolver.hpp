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
//#include "CComm.hpp" //included in CLbmSkeleton
#include "libcuda/CCL.hpp"
#include "lib/CError.hpp"
#include "libmath/CVector.hpp"
#include <typeinfo>
#include <iomanip>
#include <list>

#include "common.h"

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
    bool _isDomainSizePowOfTwo;
    bool _isLocalSizePowOfTwo;

    int _BC[3][2]; ///< Boundary conditions. First index specifys the dimension and second the upper or the lower boundary.

    CCL::CPlatform cudaPlatform;        ///< constructor will call cuInit
    CCL::CCommandQueue &cCommandQueue;
    CCL::CContext &cContext;
    CCL::CDevice &cDevice;
    CCL::CDeviceInfo cDeviceInfo;

    // INITIALIZATION KERNEL
    CCL::CKernel cKernelInit;
    size_t cKernelInit_ThreadsPerDim;
    size_t cKernelInit_GlobalWorkGroupSize;
    size_t cKernelInit_MaxRegisters;

    // COLLISION KERNELS
    CCL::CKernel cLbmKernelAlpha;
    size_t cLbmKernelAlpha_ThreadsPerDim;
    size_t cLbmKernelAlpha_GlobalWorkGroupSize;
    size_t cLbmKernelAlpha_MaxRegisters;

    CCL::CKernel cLbmKernelBeta;
    size_t cLbmKernelBeta_ThreadsPerDim;
    size_t cLbmKernelBeta_GlobalWorkGroupSize;
    size_t cLbmKernelBeta_MaxRegisters;

    // INITIALIZATION KERNEL
    CCL::CKernel cKernelCopyRect;
    size_t cKernelCopyRect_ThreadsPerDim;
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
	// initialize density buffers once
	// list of communication vector
    std::vector<CComm<T>*> _comm_container; ///< A std::Vector containing all the communcation objects for the subdomain
	std::vector< CVector<3, int> > send_size_vec;
	std::vector< CVector<3, int> > recv_size_vec;
	CCL::CMem *cStoreDensityBuffer;
	CCL::CMem *cSetDensityBuffer;
	//cStoreBuffer.create(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST);


    std::vector<std::string> split(const std::string& str,
            const std::string& delimiter = " ")
    {
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

        // set kernel args
        // source args
        cKernelCopyRect.setArg(0, src);
        cKernelCopyRect.setArg(1, src_offset); 
        cKernelCopyRect.setArg(2, src_origin[0]);
        cKernelCopyRect.setArg(3, src_origin[1]);
        cKernelCopyRect.setArg(4, src_origin[2]);
        cKernelCopyRect.setArg(5, src_size[0]);
        cKernelCopyRect.setArg(6, src_size[1]);
        cKernelCopyRect.setArg(7, src_size[2]);

        // destination args
        cKernelCopyRect.setArg(8, dst);
        cKernelCopyRect.setArg(9, dst_offset); 
        cKernelCopyRect.setArg(10, dst_origin[0]);
        cKernelCopyRect.setArg(11, dst_origin[1]);
        cKernelCopyRect.setArg(12, dst_origin[2]);
        cKernelCopyRect.setArg(13, dst_size[0]);
        cKernelCopyRect.setArg(14, dst_size[1]);
        cKernelCopyRect.setArg(15, dst_size[2]);
        cKernelCopyRect.setArg(16, block_size[0]);
        cKernelCopyRect.setArg(17, block_size[1]);
        cKernelCopyRect.setArg(18, block_size[2]);

        size_t lGlobalSize[2];
        lGlobalSize[0] = block_size[1]; 
        lGlobalSize[1] = block_size[2]; 

		cKernelCopyRect.setGridAndBlockSize(2, cKernelCopyRect_ThreadsPerDim, domain_cells_count);

        // enqueue the CopyRect kernel
        cCommandQueue.enqueueNDRangeKernel(cKernelCopyRect, // kernel
				CU_FALSE,							// shared memory use flag
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
                    p_lbm_opencl_number_of_registers_list)
    {
        CLbmSkeleton<T>::init(p_d_gravitation, p_d_viscosity, 1.0);


        if (CLbmSkeleton<T>::error())
            error << CLbmSkeleton<T>::error.getString();

        store_velocity = p_store_velocity;
        store_density = p_store_density;
        //debug = p_debug;

#if DEBUG
        // std::cout << "CL_VERSION " << _cl_version << std::endl;
#endif
        // setting the boundary conditions
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 2; j++)
                _BC[i][j] = BC[i][j];

        reload();
    }

    ~CLbmSolver()
	{
		delete[] cStoreDensityBuffer;
		delete[] cSetDensityBuffer;
	}

    void addDrivenCavityValue(T value)
    {
        drivenCavityVelocity[0] += value;
        CVector<4, T> paramDrivenCavityVelocity = drivenCavityVelocity;
        paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

        cKernelInit.setArg(4, paramDrivenCavityVelocity[0]);
        cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);
        cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);
    }

    void reload()
    {
        /**
         * WORK GROUP SIZE
         *
         * initialize the variable (postfixed appropriately with _ThreadsPerDim and _MaxRegisters)
         * with either the standard value (max_local_work_group_size) or with the value from the list
         */
#define INIT_WORK_GROUP_SIZE(variable)                                      \
        if (it != this->lbm_opencl_number_of_work_items_list.end())     \
        {                                                               \
            /*variable##_ThreadsPerDim = (*it != 0 ? *it : computation_kernel_count);*/ \
            variable##_ThreadsPerDim = (*it != 0 ? *it : threads_per_dimension); \
            it++;                                                       \
        }                                                               \
        else                                                            \
        {                                                               \
            /*variable##_ThreadsPerDim = computation_kernel_count;*/    \
            variable##_ThreadsPerDim = threads_per_dimension;           \
        }                                                               \
                                                                        \
        /* enlarge global work group size to be a multiple of collprop_work_group_size */   \
        if (domain_cells_count % variable##_ThreadsPerDim != 0)         \
            variable##_GlobalWorkGroupSize = (domain_cells_count / variable##_ThreadsPerDim + 1) * variable##_ThreadsPerDim;    \
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

        
        domain_cells_count = this->domain_cells.elements();
        _isDomainSizePowOfTwo = isDomainPowerOfTwo(domain_cells_count);
        _isLocalSizePowOfTwo = isDomainPowerOfTwo(cLbmKernelBeta_ThreadsPerDim);


        /**
         * program defines
         */
        domain_x = this->domain_cells.data[0];
        domain_y = this->domain_cells.data[1];
        domain_z = this->domain_cells.data[2];
		printf("domain: [%li, %li, %li]", domain_x, domain_y, domain_z);		

        std::ostringstream cuda_program_defines;

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
        std::cout << "RANK: " << _UID <<" BOUNDARY CONDITION: "<< std::endl;
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


        /*
         * INIT kernel
         */
        cKernelInit_GlobalWorkGroupSize = domain_cells_count;
        if (cKernelInit_GlobalWorkGroupSize % cKernelInit_ThreadsPerDim != 0)
            cKernelInit_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize
                    / cKernelInit_ThreadsPerDim + 1)
                    * cKernelInit_ThreadsPerDim;

        // Load init_lbm.ptx file with a CUDA module(program)
        CCL::CProgram cProgramInit;
        std::string initKernelModuleFileName = "src/cu_programs/lbm_init";  //cu_programs/
        initKernelModuleFileName = initKernelModuleFileName + ".ptx";
        cProgramInit.load(cContext, initKernelModuleFileName.c_str());

#if DEBUG
        std::cout << "KernelInit:   local_work_group_size: " << cKernelInit_ThreadsPerDim << "      max_registers: " << cKernelInit_MaxRegisters << std::endl;
#endif

        /*
         * ALPHA
         */
        CCL::CProgram cProgramAlpha;
        std::string alphaKernelModuleFileName = "src/cu_programs/lbm_alpha";
        alphaKernelModuleFileName = alphaKernelModuleFileName + ".ptx";
        cProgramAlpha.load(cContext, alphaKernelModuleFileName.c_str());

#if DEBUG
        std::cout << "KernelAlpha:  local_work_group_size: " << cLbmKernelAlpha_ThreadsPerDim << "      max_registers: " << cLbmKernelAlpha_MaxRegisters << std::endl;
#endif
        /*
         * BETA
         */
        CCL::CProgram cProgramBeta;
        std::string betaKernelModuleFileName = "src/cu_programs/lbm_beta";
        betaKernelModuleFileName = betaKernelModuleFileName + ".ptx";
        cProgramBeta.load(cContext, betaKernelModuleFileName.c_str());

#if DEBUG
        std::cout << "KernelBeta:   local_work_group_size: " << cLbmKernelBeta_ThreadsPerDim << "       max_registers: " << cLbmKernelBeta_MaxRegisters << std::endl;
#endif

        /*
         * INIT CopyBufferRect
         */
        CCL::CProgram cProgramCopyRect;
        std::string copyRectKernelModuleFileName = "src/cu_programs/copy_buffer_rect";
        copyRectKernelModuleFileName = copyRectKernelModuleFileName + ".ptx";
        cProgramCopyRect.load(cContext, copyRectKernelModuleFileName.c_str());

#if DEBUG
        std::cout << "KernelCopyRect:   local_work_group_size: " << cKernelCopyRect_ThreadsPerDim << "      max_registers: " << cKernelCopyRect_MaxRegisters << std::endl;
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

		// set the grid and block size for all kernels
		cKernelInit.setGridAndBlockSize(3, cKernelInit_ThreadsPerDim, domain_cells_count);
		cLbmKernelAlpha.setGridAndBlockSize(3, cLbmKernelAlpha_ThreadsPerDim, domain_cells_count);
		cLbmKernelBeta.setGridAndBlockSize(3, cLbmKernelBeta_ThreadsPerDim, domain_cells_count);
		cLbmKernelBeta_ThreadsPerDim = cLbmKernelBeta.blockSize.x * cLbmKernelBeta.blockSize.y * cLbmKernelBeta.blockSize.z;

		printf("BetaKernel --> blockSize: [%u, %u, %u] \n", cLbmKernelBeta.blockSize.x, cLbmKernelBeta.blockSize.y, cLbmKernelBeta.blockSize.z);
		printf("BetaKernel --> gridSize: [%u, %u, %u] \n", cLbmKernelBeta.gridSize.x, cLbmKernelBeta.gridSize.y, cLbmKernelBeta.gridSize.z);
 
		printf("AlphaKernel --> blockSize: [%u, %u, %u] \n", cLbmKernelAlpha.blockSize.x, cLbmKernelAlpha.blockSize.y, cLbmKernelAlpha.blockSize.z);
		printf("AlphaKernel --> gridSize: [%u, %u, %u] \n", cLbmKernelAlpha.gridSize.x, cLbmKernelAlpha.gridSize.y, cLbmKernelAlpha.gridSize.z);


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
        // cKernelInit.setArg(5, paramDrivenCavityVelocity[0]);
        cKernelInit.setArg(5, CLbmSkeleton<T>::drivenCavityVelocity[0]);
        cKernelInit.setArg(6, domain_x);
        cKernelInit.setArg(7, domain_y);
        cKernelInit.setArg(8, domain_z);

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
        // cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);
        cLbmKernelAlpha.setArg(8, CLbmSkeleton<T>::drivenCavityVelocity[0]);
        cLbmKernelAlpha.setArg(9, domain_x);
        cLbmKernelAlpha.setArg(10, domain_y);
        cLbmKernelAlpha.setArg(11, domain_z);

        cLbmKernelBeta.create(cProgramBeta, "lbm_kernel_beta"); // attach kernel to CUDA module
        cLbmKernelBeta.setArg(0, cMemDensityDistributions);
        cLbmKernelBeta.setArg(1, cMemCellFlags);
        cLbmKernelBeta.setArg(2, cMemVelocity);
        cLbmKernelBeta.setArg(3, cMemDensity);
        cLbmKernelBeta.setArg(4, this->inv_tau);
        cLbmKernelBeta.setArg(5, this->gravitation[0]);
        cLbmKernelBeta.setArg(6, this->gravitation[1]);
        cLbmKernelBeta.setArg(7, this->gravitation[2]);
        // cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);
        cLbmKernelBeta.setArg(8, CLbmSkeleton<T>::drivenCavityVelocity[0]);
        cLbmKernelBeta.setArg(9, this->domain_cells.data[0]);
        cLbmKernelBeta.setArg(10, this->domain_cells.data[1]);
        cLbmKernelBeta.setArg(11, this->domain_cells.data[2]);
        cLbmKernelBeta.setArg(12, cLbmKernelBeta_ThreadsPerDim);

//		printf("cLbmKernelBeta_ThreadsPerDim: %u\n", cLbmKernelBeta_ThreadsPerDim);
//		printf("(cLbmKernelBeta.blockSize.x * cLbmKernelBeta.blockSize.y * cLbmKernelBeta.blockSize.z): %u \n", (cLbmKernelBeta.blockSize.x * cLbmKernelBeta.blockSize.y * cLbmKernelBeta.blockSize.z));
//		printf("blockSize: [%u, %u, %u] \n", cLbmKernelBeta.blockSize.x, cLbmKernelBeta.blockSize.y, cLbmKernelBeta.blockSize.z);
//		printf("gridSize: [%u, %u, %u] \n", cLbmKernelBeta.gridSize.x, cLbmKernelBeta.gridSize.y, cLbmKernelBeta.gridSize.z);
        cLbmKernelBeta.setArg(13, _isDomainSizePowOfTwo);
        cLbmKernelBeta.setArg(14, _isLocalSizePowOfTwo);


        cKernelCopyRect.create(cProgramCopyRect, "copy_buffer_rect"); // attach kernel to CUDA module

        reset();
    }

    void reset()
    {
        simulation_step_counter = 0;

#if DEBUG
        std::cout << "Init Simulation: " << std::flush;
#endif
        cCommandQueue.enqueueNDRangeKernel(cKernelInit, ///< kernel
				CU_FALSE,								// shared memory use flag
                cKernelInit.kernelArgsVec);

//		printf("rank: %i, initKernel: after launch\n", _UID);
        cCommandQueue.enqueueBarrier();

#if DEBUG
        std::cout << "OK" << std::endl;
#endif
    }


    void simulationStepAlpha()
    {
#if DEBUG
        // std::cout << "--> Running Alpha kernel" << std::endl;
#endif
        cCommandQueue.enqueueNDRangeKernel(
                cLbmKernelAlpha,    // kernel
				CU_FALSE,			// shared memory use flag
                cLbmKernelAlpha.kernelArgsVec);
    }

    void simulationStepBeta()
    {
#if DEBUG
        // std::cout << "--> Running BETA kernel" << std::endl;
#endif
        cCommandQueue.enqueueNDRangeKernel(
                cLbmKernelBeta, 	// kernel
				CU_TRUE,			// shared memory use flag
                cLbmKernelBeta.kernelArgsVec);
	
    }

    /**
     * start one simulation step (enqueue kernels)
     */
    void simulationStep() 
    {
        /*
         * collision kernels are inserted as 1d kernels because they work cell wise without neighboring information
         */
        if (simulation_step_counter & 1) {
            simulationStepAlpha();
        } else {
            simulationStepBeta();
        }
        cCommandQueue.enqueueBarrier();

        simulation_step_counter++;
    }

	// accept communication vector from Controller
    void setGhostLayerBuffers(std::vector<CComm<T> *> &comm)
	{
		this->_comm_container = comm;
        typename std::vector<CComm<T>*>::iterator it = this->_comm_container.begin();

		// add the CVector objects to the vector container
		for(; it != this->_comm_container.end(); it++)
		{
			send_size_vec.push_back( (*it)->getSendSize() );	
			recv_size_vec.push_back( (*it)->getRecvSize() );	
		}

		// number of CMem objects to be created for send and recv buffers
		int i = 0;
		size_t storeMemObjects = send_size_vec.size(); 
		size_t setMemObjects = recv_size_vec.size();
		
		cStoreDensityBuffer = new CCL::CMem[storeMemObjects]; 
		cSetDensityBuffer = new CCL::CMem[setMemObjects]; 

		// allocate the buffer objects by accessing the CVector container objects
        typename std::vector<CVector<3, int> >::iterator it_s = this->send_size_vec.begin();
        typename std::vector<CVector<3, int> >::iterator it_r = this->recv_size_vec.begin();
		//for(; it_s != this->send_size_vec.end(); it_s++, i++)
		it = this->_comm_container.begin();
		for(; it != this->_comm_container.end(); it++, it_s++, it_r++, i++)
		{
			cStoreDensityBuffer[i] = CCL::CMem();
			cSetDensityBuffer[i] = CCL::CMem();
			cStoreDensityBuffer[i].create(cContext, sizeof(T) * (*it_s).elements() * SIZE_DD_HOST);
			cSetDensityBuffer[i].create(cContext, sizeof(T) * (*it_r).elements() * SIZE_DD_HOST);
#if DEBUG
/*			std::cout << "Solver --> rank: " << _UID ", send_size: " << (*it_s).elements() << ", recv_size: " << (*it_r).elements() << std::endl;
			std::cout << "Solver --> send_elems: " << (*it_s).data[0] << ", " << (*it_s).data[1] << ", " << (*it_s).data[2] << std::endl;
			std::cout << "Solver --> recv_elems: " << (*it_r).data[0] << ", " << (*it_r).data[1] << ", " << (*it_r).data[2] << std::endl;
*/
#endif
		}
	}

    /**
     * wait until all kernels have finished
     */
    void wait()
    {
        cCommandQueue.finish();
    }

    /**
     * store density distribution values to host memory
     */
    void storeDensityDistribution(T *dst)
    {
        size_t byte_size = cMemDensityDistributions.getSize();

        cCommandQueue.enqueueReadBuffer(cMemDensityDistributions, CU_TRUE, // sync reading
                0, byte_size, dst);
    }

    void storeDensityDistribution(T *dst, CVector<3, int> &origin,
            CVector<3, int> &size, int i)
    {
        //CCL::CMem cBuffer;
        //cBuffer.create(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST);

        for (int f = 0; f < SIZE_DD_HOST; f++)
        {
            //enqueueCopyRectKernel(cMemDensityDistributions, cBuffer,
            enqueueCopyRectKernel(cMemDensityDistributions, cStoreDensityBuffer[i],
                    f * this->domain_cells_count, origin,
                    this->domain_cells, f * size.elements(),
                    CVector<3, int>(0, 0, 0), size, size, false);
        }
        cCommandQueue.enqueueBarrier();

        //cCommandQueue.enqueueReadBuffer(cBuffer, CU_TRUE, 0, cBuffer.getSize(), dst);
        cCommandQueue.enqueueReadBuffer(cStoreDensityBuffer[i], CU_TRUE, 
			0, cStoreDensityBuffer[i].getSize(), dst);
    }

    void setDensityDistribution(T *src, CVector<3, int> &origin,
            CVector<3, int> &size, int i)
    {
        
        //CCL::CMem cBuffer;
        //cBuffer.create(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST);
		//cCommandQueue.enqueueWriteBuffer( cBuffer, CU_TRUE, 0, sizeof(T) * size.elements() * SIZE_DD_HOST, src);
		cCommandQueue.enqueueWriteBuffer( cSetDensityBuffer[i], CU_TRUE, 0,
			sizeof(T) * size.elements() * SIZE_DD_HOST, src);

        for (int f = 0; f < SIZE_DD_HOST; f++)
        {
            //enqueueCopyRectKernel(cBuffer, cMemDensityDistributions,
            enqueueCopyRectKernel(cSetDensityBuffer[i], cMemDensityDistributions,
                    f * size.elements(), CVector<3, int>(0, 0, 0), size,
                    f * this->domain_cells_count, origin,
                    this->domain_cells, size, false);
        }
        cCommandQueue.enqueueBarrier();
    }

    void setDensityDistribution(T *src, CVector<3, int> &origin,
            CVector<3, int> &size, CVector<3, int> norm, int i)
    {
        //CCL::CMem cBuffer;
        //cBuffer.create(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST);
		//cCommandQueue.enqueueWriteBuffer( cBuffer, CU_TRUE, 0, sizeof(T) * size.elements() * SIZE_DD_HOST, src);
		cCommandQueue.enqueueWriteBuffer( cSetDensityBuffer[i], CU_TRUE, 0, 
			sizeof(T) * size.elements() * SIZE_DD_HOST, src);

        for (int f = 0; f < SIZE_DD_HOST; f++)
        {
            if (norm.dotProd(lbm_units[f]) > 0)
            {
                //enqueueCopyRectKernel(cBuffer, cMemDensityDistributions,
                enqueueCopyRectKernel(cSetDensityBuffer[i], cMemDensityDistributions,
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
    void storeVelocity(T *dst)
    {
        size_t byte_size = cMemVelocity.getSize();

        cCommandQueue.enqueueReadBuffer(cMemVelocity, CU_TRUE, // sync reading
             0, byte_size, dst);
    }

    /**
     * Store a block of velocity data from device to host
     *
     * @param dst The buffer that will contain the return values
     * @param origin The origin point of data block
     * @param size The size of data block
     */
    void storeVelocity(T* dst, CVector<3, int> &origin, CVector<3, int> &size)
    {
        CCL::CMem cBuffer;

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
    void setVelocity(T* src, CVector<3, int> &origin, CVector<3, int> &size)
    {
        CCL::CMem cBuffer;

        cBuffer.createCopyToDevice(cContext, 
                sizeof(T) * size.elements() * 3, src);

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
    void storeDensity(T *dst)
    {
        size_t byte_size = cMemDensity.getSize();

        cCommandQueue.enqueueReadBuffer(cMemDensity, CU_TRUE, // sync reading
                0, byte_size, dst);

    }

    /**
     * Store a block of density data from device to host
     *
     * @param dst The buffer that will contain the return values
     * @param origin The origin point of data block
     * @param size The size of data block
     */
    void storeDensity(T *dst, CVector<3, int> &origin, CVector<3, int> &size)
    {
        CCL::CMem cBuffer;

        cBuffer.create(cContext, sizeof(T) * size.elements());

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

        enqueueCopyRectKernel(cBuffer, cMemDensity, 0,
                CVector<3, int>(0, 0, 0), size, 0, origin,
                this->domain_cells, size);
    }

    void storeFlags(int *dst)
    {
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

        // printf("src [cBuffer -> cMemCellFlags]: size= %i\n", size.elements() );
        // for (size_t i = 0; i < size.elements(); i++)
        // {
        //     printf(" %i ", src[i]);
        //     if ((i + 1) % (domain_x) == 0)
        //     {
        //         printf("\n");
        //     }
        // }
        // printf("\n");

        enqueueCopyRectKernel(cBuffer, cMemCellFlags, 0,
                CVector<3, int>(0, 0, 0), size, 0, origin,
                this->domain_cells, size);
    }

private:
    bool isDomainPowerOfTwo(size_t x)
    {
        return ( (x == (1<<0)) || (x == (1<<1)) || (x == (1<<2)) || (x == (1<<3)) || (x == (1<<4)) ||
                 (x == (1<<5)) || (x == (1<<6)) || (x == (1<<7)) || (x == (1<<8)) || (x == (1<<9)) || 
                 (x == (1<<10)) || (x == (1<<11)) || (x == (1<<12)) || (x == (1<<13)) || (x == (1<<14)) ||
                 (x == (1<<15)) || (x == (1<<16)) || (x == (1<<17)) || (x == (1<<18)) || (x == (1<<19)) ||
                 (x == (1<<20)) || (x == (1<<21)) || (x == (1<<22)) || (x == (1<<23)) || (x == (1<<24)) ||
                 (x == (1<<25)) || (x == (1<<26)) || (x == (1<<27)) || (x == (1<<28)) || (x == (1<<29)) || 
                 (x == (1<<30)) || (x == (1<<31)) );
    }

    void debugChar(CCL::CMem &cMem, size_t wrap_size = 20)
    {
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

    void debugFloat(CCL::CMem &cMem, size_t wrap_size = 20)
    {
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
    void debug_print()
    {
        // for (int i = 0; i < 2; i++)
        // {
        //     if (_UID == i)
        //     {
                // read out DensityDistributions
                // std::cout << "Rank: " << _UID << " DENSITY DISTRIBUTIONS:";
                // //debugFloat(cMemDensityDistributions, SIZE_DD_HOST);
                // debugFloat(cMemDensityDistributions, 16);
                // std::cout << std::endl;

                // read out Velocity
                std::cout << std::endl;
                std::cout << "Rank: " << _UID << " VELOCITY:";
                debugFloat(cMemVelocity, 4 * 3);
                std::cout << std::endl;

                // read out Density
                std::cout << std::endl;
                std::cout << "Rank: " << _UID << " DENSITY:";
                debugFloat(cMemDensity, 4);
                std::cout << std::endl;

                // read out Flags
                std::cout << std::endl;
                std::cout << "Rank: " << _UID << " FLAGS:";
                debugChar(cMemCellFlags, 4 * 4);
                std::cout << std::endl;
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        // }
    }

    /**
     * debug density distributions (works only in linear mode!)
     * \param   dd_id   specifies the dd number
     */
    void debugDD(size_t dd_id = 0, size_t wrap_size = 16,
            size_t empty_line = 16)
    {
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

    float getVelocityChecksum()
    {
        T *velocity = new T[this->domain_cells.elements() * 3];
        int *flags = new int[this->domain_cells.elements() * 3];

        storeVelocity(velocity);
        storeFlags(flags);

        T *velx = velocity;
        T *vely = velocity + this->domain_cells.elements();
        T *velz = velocity + this->domain_cells.elements() * 2;

        float checksum = 0;
        for (int a = 0; a < this->domain_cells.elements(); a++)
		{
//			if (_UID == 0)
//				printf("D: %i, [%f, %f, %f] \n", a, velx[a], vely[a], velz[a]);

            if (flags[a] == FLAG_FLUID)
            {
                checksum += velx[a] + vely[a] + velz[a];
            }
		}

        delete[] flags;
        delete[] velocity;

        return checksum;
    }
};

#endif
