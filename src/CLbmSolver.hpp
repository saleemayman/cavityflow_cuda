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

#ifndef CLBMSOLVER_HPP
#define CLBMSOLVER_HPP

#include <list>
#include <vector>

#include "CLbmSkeleton.hpp"
//#include "CComm.hpp" //included in CLbmSkeleton
#include "libcuda/CCL.hpp"
#include "lib/CError.hpp"
#include "libmath/CVector.hpp"

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

	// initialize density buffers once
	// list of communication vector
    std::vector<CComm<T>*> _comm_container; ///< A std::Vector containing all the communcation objects for the subdomain
	std::vector< CVector<3, int> > send_size_vec;
	std::vector< CVector<3, int> > recv_size_vec;
	CCL::CMem *cStoreDensityBuffer;
	CCL::CMem *cSetDensityBuffer;
	//cStoreBuffer.create(cContext, sizeof(T) * size.elements() * SIZE_DD_HOST);

    std::vector<std::string> split(const std::string& str,
            const std::string& delimiter = " ");
    inline void enqueueCopyRectKernel(CCL::CMem& src, CCL::CMem& dst,
            const int src_offset, CVector<3, int> src_origin,
            CVector<3, int> src_size, const int dst_offset,
            CVector<3, int> dst_origin, CVector<3, int> dst_size,
            CVector<3, int> block_size, bool withBarrier = true);
    bool isDomainPowerOfTwo(size_t x);
    void debugChar(CCL::CMem &cMem, size_t wrap_size = 20);
    void debugFloat(CCL::CMem &cMem, size_t wrap_size = 20);

public:
    CVector<4, T> drivenCavityVelocity;
    static const size_t SIZE_DD_HOST = 19;
    size_t simulation_step_counter;
    CError error;
    std::list<int> lbm_opencl_number_of_work_items_list; ///< list with number of threads for each successively created kernel
    std::list<int> lbm_opencl_number_of_registers_list; ///< list with number of registers for each thread threads for each successively created kernel

    CLbmSolver(int UID, CCL::CCommandQueue &p_cCommandQueue,
            CCL::CContext &p_cContext, CCL::CDevice &p_cDevice, int BC[3][2],
            CDomain<T> &domain, CVector<3, T> &p_d_gravitation, T p_d_viscosity,
            size_t p_computation_kernel_count,
            size_t p_threads_per_dimension,
            //bool p_debug,
            bool p_store_velocity, bool p_store_density, T p_d_timestep,
            CVector<4,T>& _drivenCavityVelocity,
            std::list<int> &p_lbm_opencl_number_of_work_items_list, ///< list with number of threads for each successively created kernel
            std::list<int> &p_lbm_opencl_number_of_registers_list);
    ~CLbmSolver();

    void addDrivenCavityValue(T value);
    void reload();
    void reset();
    void simulationStepAlpha();
    void simulationStepBeta();
    void simulationStep();
    void setGhostLayerBuffers(std::vector<CComm<T> *> &comm);
    void wait();
    void storeDensityDistribution(T *dst);
    void storeDensityDistribution(T *dst, CVector<3, int> &origin,
            CVector<3, int> &size, int i);
    void setDensityDistribution(T *src, CVector<3, int> &origin,
            CVector<3, int> &size, int i);
    void setDensityDistribution(T *src, CVector<3, int> &origin,
            CVector<3, int> &size, CVector<3, int> norm, int i);
    void storeVelocity(T *dst);
    void storeVelocity(T* dst, CVector<3, int> &origin, CVector<3, int> &size);
    void setVelocity(T* src, CVector<3, int> &origin, CVector<3, int> &size);
    void storeDensity(T *dst);
    void storeDensity(T *dst, CVector<3, int> &origin, CVector<3, int> &size);
    void setDensity(T *src, CVector<3, int> &origin, CVector<3, int> &size);
    void storeFlags(int *dst);
    void storeFlags(int *dst, CVector<3, int> &origin, CVector<3, int> &size);
    void setFlags(int *src, CVector<3, int> &origin, CVector<3, int> &size);
    void debug_print();
    void debugDD(size_t dd_id = 0, size_t wrap_size = 16, size_t empty_line = 16);
    float getVelocityChecksum();
};

#endif
