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


#ifndef CCOMPLANG_HPP
#define CCOMPLANG_HPP

#include <iostream>
#include <strings.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <sstream>
#include <cmath>
#include <vector>
#include "CCLErrors.hpp"
#include "../lib/CError.hpp"
#include "../lib/CFile.hpp"
#include "../libmath/CVector.hpp"
#include "../common.h"
#include <sys/types.h>

// load defaults kernel execution paramters
#define DEFAULTS (1)
#include "../cu_programs/lbm_defaults.h"

#ifdef C_GL_TEXTURE_HPP
    #if GL3_PROTOTYPES
        #define __gl_h_
    #endif
    #include <CL/cl_gl.h>

    #ifndef CL_GL_CONTEXT_KHR
        #define     CL_GL_CONTEXT_KHR       0x2008
        #define     CL_EGL_DISPLAY_KHR      0x2009
        #define     CL_GLX_DISPLAY_KHR      0x200A
        #define     CL_WGL_HDC_KHR          0x200B
        #define     CL_CGL_SHAREGROUP_KHR   0x200C
    #endif

    #include <GL/glx.h>
#endif

/**
 * \brief OpenCL C++ abstraction class
 */
class CCL
{
public:

    /**
     * \brief Load and get information about an OpenCL platform
     */
    class CPlatform
    {
    public:
        char *version;  ///< version string
        char *name;     ///< name string
        char *vendor;   ///< vendor string
        char *extensions;   ///< extensions string

private:
        /**
         * initialize data structure with default values (NULL)
         */
        inline void initCPlatform()
        {
            version = NULL;
            name = NULL;
            vendor = NULL;
            extensions = NULL;

            // cuda driver API initialization
            if (cuInit(0) != CUDA_SUCCESS)
            {
                std::cerr << "CUDA Driver API could not initialize, exiting!" << std::endl;
                exit(-1);
            }
        }

        /**
         * cleanup platform (free strings)
         */
        inline void cleanupCPlatform()
        {
            delete[] version;
            delete[] name;
            delete[] vendor;
            delete[] extensions;
        }

public:
        inline void load()
        {
            initCPlatform();
        }

        inline CPlatform()
        {
            initCPlatform();
        }

        inline ~CPlatform()
        {
            cleanupCPlatform();
        }
    };


    /**
     * \brief Handle an OpenCL device
     */
    class CDevice
    {
    public:
        CUdevice device_id;             ///< CUDA Device ID

        /**
         * initialize CDevice class with existing device id
         */
        inline CDevice(CUdevice p_device_id)
        {
            device_id = p_device_id;
        }
        /**
         * initialize CDevice class with existing class CDevice
         */
        inline CDevice(const CDevice &cDevice)
        {
            device_id = cDevice.device_id;
        }

        /**
         * default constructor
         */
        inline CDevice()
        {
            device_id = 0;
        }

        /**
         * initialize CDevice class with existing class CDevice
         */
        inline void initWithExistingDevice(const CDevice &cDevice)
        {
            device_id = cDevice.device_id;
        }

        /**
         * set device_id to p_device_id
         */
        inline void set(CUdevice p_device_id)
        {
            device_id = p_device_id;
        }

        /**
         * set device_id to device_id in cDevice
         */
        inline void set(const CDevice &cDevice)
        {
            device_id = cDevice.device_id;
        }
    };

    class CContext;
    /**
     * \brief Load a list of devices for different OpenCL contexts, platforms and/or types
     */
    class CDevices  : public std::vector<CDevice>
    {
    public:
        CUcontext context;
        CUdevice *device_ids;
        int device_ids_count;       ///< number of contexts in device_ids[]

        /**
         * load device list belonging to context cContext
         * \param cContext  existing context
         */
         inline void load()
        {
            // load device information
            CudaCallCheckError( cuDeviceGetCount(&device_ids_count) );
            if (device_ids_count == 0)  std::cerr << "Warning: no device found!" << std::endl;
            delete device_ids;
            device_ids = new CUdevice[device_ids_count];

            /// get a handle for each device available
            for (int i = 0; i < device_ids_count; i++)
            {
                CudaCallCheckError( cuDeviceGet(&device_ids[i], i) );   
            }

            initDeviceVector();
        }

        /**
         * initialize device list with NULL
         */
        inline void initCDevices()
        {
            device_ids = NULL;
            device_ids_count = 0;
            clear();
        }
        /**
         * load device list belonging to cContext
         */
        inline CDevices()
        {
            initCDevices();
            load();
        }

        /**
         * deconstructor: free unnecessary device data
         */
        inline ~CDevices()
        {
            if (device_ids != NULL)
                delete[] device_ids;
            clear();
        }

    private:
        inline void initDeviceVector()
        {
            (*this).resize(device_ids_count);
            for (uint i = 0; i < device_ids_count; i++)
            {
                (*this)[i].set(device_ids[i]);
            }
        }
    };



    /**
     * \brief Create and manage an OpenCL context
     */
    class CContext
    {
    public:
        CUresult error;     ///< error handler
        CUcontext context;  ///< OpenCL context id

        /**
         * load context by platform and device
         */
        inline void load(   const CDevice &cDevice      ///< CUDA device
        )
        {
            release();

            CUresult err_ret;
            err_ret = cuCtxCreate(&context, CU_CTX_SCHED_AUTO, cDevice.device_id);
            CudaCallCheckError(err_ret);
            retain();
        }

// ask Arash about creating a context from a list of devices!
#if 0
        /**
         * load context by platform and device list
         */
        inline void load(   const CPlatform &cPlatform, ///< platform for parameters
                            const CDevices &cDevices    ///< device list available for context
        )
        {
            release();

            cl_int err_ret;
            cl_context_properties context_properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id, 0, 0};
            context = clCreateContext(context_properties, cDevices.device_ids_count, cDevices.device_ids, NULL, NULL, &err_ret);
            CL_CHECK_ERROR(err_ret);

            retain();
        }
#endif


        /**
         * initialize with NULL context
         */
        inline void initCContext()
        {
            context = NULL;
        }

        /**
         * initialize with NULL context
         */
        inline CContext()
        {
            context = NULL;
        }

        /**
         * initialize and load context with device and platform
         */
        inline CContext(    CDevice &cDevice        ///< OpenCL device
        )
        {
            load(cDevice);
        }


        inline ~CContext()
        {
            release();
        }


        /**
         * initialize context handler with existing context
         */
        inline void initWithExistingContext(const CContext &cContext)
        {
            release();
            context = cContext.context;
            retain();
        }


private:
        /**
         * increment reference counter to context to avoid deletion of context if
         * multiple classes use the same context
         */
        inline void retain()
        {
            CudaCallCheckError( cuCtxAttach(&context, 0) );
        }

        /**
         * decrement reference counter to context if context is valid
         */
        inline void release()
        {
            if (context != NULL)
            {
                cuCtxDetach(context);
                context = NULL;
            }
        }
    };

    /**
     * \brief   handle multiple platforms available for computations
     */
    class CPlatforms
    {
    public:
        CUresult err_ret;

        /**
         * initialize with NULL values
         */
        inline void initCPlatforms()
        {
            err_ret = cuInit(0);
            if (err_ret != CUDA_SUCCESS)
            {
                CudaCallCheckError(err_ret);
                std::cerr << "Warning: CUDA Driver API not initialized!" << std::endl;
            }
        }

        /**
         * load all available platforms
         */
        inline void load()
        {
            initCPlatforms();
        }

        /**
         * default constructor: load all available platforms
         */
        inline CPlatforms()
        {
            initCPlatforms();
        }

        inline ~CPlatforms()
        {
        }
    };

    /**
     * \brief handler for OpenCL buffer objects
     *
     * a buffer object is initialized via the CContextDevices class
     */
    class CMem
    {
    public:
        CUdeviceptr memobj;

        inline CMem()
        {
            memobj = '\0';  //NULL;
        }

        /**
         * initialize reference of memory object with existing reference
         */
        inline CMem(const CMem &bo)
        {
            memobj = bo.memobj;
        }

        /**
         * create OpenCL memory buffer
         */
        inline CMem(    CContext&       cContext,   ///< context for buffer
                        size_t          size        ///< Size of memory buffer
            )
        {
            memobj = '\0';  //NULL;
            create(cContext, size); 
        }

        /**
         * release memory object (decrements OpenCL reference counter)
         */
        void release()
        {
            if (memobj != '\0')     //if (memobj != NULL)
                CudaCallCheckError( cuMemFree(memobj) );
            memobj = '\0';  //NULL;
        }

        /**
         * deconstructor (release)
         */
        inline ~CMem()
        {
            release();
        }

        /**
         * return the size of the memory object
         */
        inline size_t getSize()
        {
            CUdeviceptr pbase;
            size_t mem_size;
            CudaCallCheckError( cuMemGetAddressRange(&pbase, &mem_size, memobj) );
            return mem_size;        ///< returns the size of memobj in size_t
        }

        /**
         * create new memory object
         */
        inline void create( CContext&       cContext,   ///< context for buffer
                            size_t          size        ///< Size of memory buffer
        )
        {
            CUresult errcode_ret;

            release();
            errcode_ret = cuMemAlloc(&memobj, size);
            CudaCallCheckError(errcode_ret);
        }

        /*
         * same as create() but copies the data from host to device
         * OpenCL equivalent is clCreateBuffer with  CL_MEM_COPY_HOST_PTR  flag
         */
        inline void createCopyToDevice( CContext&       cContext,   ///< context for buffer
                                        size_t          size,       ///< Size of memory buffer
                                        void*           host_ptr    ///< Host pointer to initialize buffer
        )
        {
            // Allocate GPU memory
            create(cContext, size);

            // Copy data from host_ptr to GPU
            CudaCallCheckError( cuMemcpyHtoD(memobj, host_ptr, size) );
        }

        // for debugging
        inline void printVelMemObj(size_t elems)
        {
            T *vel_memobj_h = new T[elems];

            // copy data to host buffer
            CudaCallCheckError( cuMemcpyDtoH(vel_memobj_h, memobj, elems * sizeof(T)) );

            printf("Velocity memObj size: %lu \n", elems);
            for (int i = 0; i < elems; i+=10)
            {
                printf("\t %f", vel_memobj_h[i]);
            }
            printf("\n");

            delete[] vel_memobj_h;
        }

        // for debugging
        inline void printDensityMemObj(size_t elems)
        {
            T *rho_memobj_h = new T[elems];

            // copy device densities to host
            CudaCallCheckError( cuMemcpyDtoH(rho_memobj_h, memobj, elems * sizeof(T)) );

            printf("Density memObj size: %lu \n", elems);
            for (int i = 0; i < elems; i+=10)
            {
                printf("\t %f\n", rho_memobj_h[i]);
            }
            printf("\n");

            delete[] rho_memobj_h;
        }
/*
 * CGlTexture has to be included before CCL.hpp to enable the creation of
 * OpenCL memory handlers from OpenGL textures
 */
#ifdef C_GL_TEXTURE_HPP
        /**
         * create OpenCL memory object from existing 3D OpenGL texture
         */
        inline void createFromGLTexture3D(
                CContext&       cContext,   ///< context for buffer
                cl_mem_flags    flags,      ///< OpenCL flags
                CGlTexture      &texture    ///< reference to 3d texture
        )
        {
            cl_int errcode_ret;
            memobj = clCreateFromGLTexture3D(cContext.context, flags, GL_TEXTURE_3D, 0, texture.textureid, &errcode_ret);
            CL_CHECK_ERROR(errcode_ret);
        }

        /**
         * create OpenCL memory object from existing 2D OpenGL texture
         */
        inline void createFromGLTexture2D(
                CContext&       cContext,   ///< context for buffer
                cl_mem_flags    flags,      ///< OpenCL flags
                CGlTexture      &texture    ///< reference to 3d texture
        )
        {
            cl_int errcode_ret;
            memobj = clCreateFromGLTexture2D(cContext.context, flags, texture.target, 0, texture.textureid, &errcode_ret);
            CL_CHECK_ERROR(errcode_ret);
        }
#endif
    };


    /**
     * \brief Load, use and manage OpenCL programs to create kernels
     */
    class CProgram
    {
private:
        std::string command;    ///< command string to execute
        std::stringstream command_string_stream;    ///< string to append
        const char *file_name;  ///< cuda file to compile

public:
        CError error;           ///< error handler
        CUmodule program;       ///< CUDA program id equivalent
        std::string filepath;   ///< string with filename (for debugging purposes)

        inline CProgram()
        {
            program = 0;
        }

        /**
         * load program given by p_filepath
         */
        inline CProgram(    CContext &cContext, // context
                            const std::string &p_filepath   // source file
        )
        {
            CProgram();
            load(cContext, p_filepath.c_str() );
        }

        /**
         * load kernel using source given in strings
         */
        inline void load(   CContext &cContext, // context
                            const char *cudaKernelFileName  // cuda cubin/ptx source code
        )
        {
            CUresult errcode_ret;
            errcode_ret = cuModuleLoad(&program, cudaKernelFileName);
            if (errcode_ret != CUDA_SUCCESS)
            {
                CudaCallCheckError(errcode_ret);
                CudaCallCheckError( cuCtxDestroy(cContext.context) );
            }
        }


        /**
         * execute compile command on system
         */
        inline void executeCommand( const std::string &command, 
                                    const char *file_name
        )
        {
            int err_ret;
            std::cout << "Executing Compile command --> " << command << std::endl;
            err_ret = system(command.c_str());

            if (err_ret != 0)
            {
                std::cerr << "failed to compile " << file_name << std::endl;
                exit(err_ret);
            }
            else
            {
                std::cout << "Compile Successful, file --> " << file_name << std::endl;
            }
            // printf("\n");
        }
    };


    /**
     * \brief Create kernels from CProgram, use and set arguments for them
     */
    class CKernel
    {
    public:
        CUfunction kernel;  ///< CUDA funtion handle
        CUresult error;     ///< cuda error handler
        std::vector<void *> kernelArgsVec;  ///< vector to hold list of kernel input parameters

		dim3 blockSize; 		///< threads per grid dimension = dim3(32, 1, 1);
		dim3 gridSize;  		///< number of blocks to launch = dim3((size + blockSize.x - 1) / blockSize.x, 1, 1);

        /**
         * create kernel from OpenCL program
         */
        inline void create( CProgram &program,          ///< OpenCL program
                            const char *kernel_name     ///< name of kernel function to create kernel from
            )
        {
            // printf("kernel file name: %s\n", kernel_name);

            CudaCallCheckError( cuModuleGetFunction(&kernel, program.program, kernel_name) );
        }

        /**
         * create kernel from OpenCL program
         */
        inline CKernel( CProgram &program,          ///< OpenCL program
                        const char *kernel_name     ///< name of kernel function to create kernel from
        )
        {
            create(program, kernel_name);
        }

        /**
         * create NULL kernel
         */
        inline CKernel()
        {
            kernel = NULL;
            kernelArgsVec.reserve(20);
        }

        /**
         * create kernel from existing kernel
         */
        inline CKernel(const CKernel &k)
        {
            kernel = k.kernel;
        }

        /**
         * default deconstructor:
         * release kernel
         */
        inline ~CKernel()
        {
        }

        /**
         * allocate enough space for the kernel arguments in std::vector container
         */
        inline void setArgSize(int size)
        {
            kernelArgsVec.reserve(size);
        }
        /**
         * set memory object kernel argument
         */
        inline void setArg( int arg_index, ///< argument number
                            CMem &cmem ///< OpenCL memory object
        )
        {
            kernelArgsVec[arg_index] = &(cmem.memobj);
        }

        /**
         * set cl_float kernel argument
         */
        inline void setArg( int arg_index,
                            float &arg      ///< reference to float argument
        )
        {
            kernelArgsVec[arg_index] = &arg;
        }

        /**
         * set cl_double kernel argument
         */
        inline void setArg( int arg_index,
                            double &arg     ///< reference to float argument
        )
        {
            kernelArgsVec[arg_index] = &arg;
        }

        /**
         * set cl_int kernel argument
         */
        inline void setArg( int arg_index,
                            int &arg
        )
        {
            kernelArgsVec[arg_index] = &arg;
        }

        /**
         * set size_t kernel argument
         */
        inline void setArg( int arg_index,
                            bool &arg
        )
        {
            kernelArgsVec[arg_index] = &arg;
        }

        /**
         * set size_t kernel argument
         */
        inline void setArg( int arg_index,
                            size_t &arg
        )
        {
            kernelArgsVec[arg_index] = &arg;
        }


        inline void setArg( int arg_index,
                            const int &arg
        )
        {
            kernelArgsVec[arg_index] = const_cast<int *>(&arg); //const_cast < new_type > ( expression )        
        }

        /**
         * set the number of grid size and threads per block for GPU
         */
        void setGridAndBlockSize(	int work_dim_def,
									CVector<3, int> &threads_per_dim, const size_t total_elements
		)
        {
            unsigned int threadsPerBlock;
            unsigned int grids_x, grids_y, grids_z;

			blockSize = dim3(threads_per_dim[0], threads_per_dim[1], threads_per_dim[2]);
			threadsPerBlock = blockSize.x * blockSize.y * blockSize.z;

			// get the block/grid dimension implicitly from threads_per_dim
			int work_dim = 3;
			for (int i = 0; i < 3; i++)
			{
				if (threads_per_dim[i] == 1)
					work_dim--;
			}

			// use default work dimensions if specified
			if (work_dim_def != NULL)		
				work_dim = work_dim_def;

            if (work_dim == 1)
            {
                grids_x = (total_elements + threadsPerBlock - 1)/threadsPerBlock;

                gridSize = dim3(grids_x + 1, 1, 1);
            }
            else if(work_dim == 2)
            {
                grids_x = sqrt((total_elements + threadsPerBlock - 1) / threadsPerBlock);
                grids_y = grids_x;

                gridSize = dim3(grids_x + 1, grids_y + 1, 1);
            }
            else
            {
//				blockSize = dim3(threads_per_dim[0], threads_per_dim[1], threads_per_dim[2]);
                blockSize = dim3(threads_per_dim[0], threads_per_dim[1], 1);
                threadsPerBlock = blockSize.x * blockSize.y * blockSize.z;

				grids_x = pow((total_elements + threadsPerBlock - 1)/threadsPerBlock, 1./3.);
                grids_y = grids_x;
                grids_z = grids_x;

                gridSize = dim3(grids_x + 1, grids_y + 1, grids_z + 1);
            }
        }

    };

    /**
     * \brief OpenCL event handler
     */
    class CEvent
    {
    public:
        //cl_event event;   ///< openCL event
        CUevent event;      ///< CUDA event

        inline CEvent()
        {
            event = 0;
        }

        /**
         * wait for event to complete
         */
        inline void waitAndRelease()
        {
            if (event != 0)
            {
                CudaCallCheckError( cuEventSynchronize(event) );
                release();
                event = 0;
            }
        }

        /**
         * return the execution status belonging to the event
         */
        inline CUresult getExecutionStatus()
        {
            CUresult param_value;
            param_value = cuEventQuery(event);
            CudaCallCheckError(param_value);

            return param_value;
        }

        inline ~CEvent()
        {
            if (event != 0)
                CudaCallCheckError( cuEventDestroy(event) );
        }

        /**
         * release the event
         */
        inline void release()
        {
            CudaCallCheckError( cuEventDestroy(event) );
            event = 0;
        }
    };


    /**
     * \brief OpenCL command queue handler
     */
    class CCommandQueue
    {
    public:
		CError error;       	///< error handler
        CUstream cuda_stream;   ///< CUDA stream handler

        /**
         * create cuda stream
         */
        inline void create()
        {
            CudaCallCheckError( cuStreamCreate(&cuda_stream, 0) );  // flag has to be zero
        }

        /**
         * decrement OpenCL reference counter to command queue
         */
        inline void release()
        {
            if (cuda_stream != 0)   // check if this is correct
            {
                CudaCallCheckError(cuStreamDestroy(cuda_stream));
            }
        }

        inline CCommandQueue()
        {
            cuda_stream = 0;
            create();
        }

        inline ~CCommandQueue()
        {
            release();
        }

        /**
         * enqueue barrier to command queue
         */
        inline void enqueueBarrier()
        {
            CUresult errcode_ret;   // cuda error type
            errcode_ret = cuStreamSynchronize(cuda_stream);
            if (errcode_ret != CUDA_SUCCESS)
                error << getCudaDrvErrorString(errcode_ret) << std::endl;   // CUDA Driver API erro handling

        }

        /**
         * copy data from buffer in host memory to CL device
         */
        inline void enqueueWriteBuffer( CMem &cMem,
                                        bool block_write,
                                        size_t offset,
                                        size_t buffer_size,
                                        const void *buffer_ptr
        )
        {
            // TODO: copy data to device
            if (block_write)    // sync write
            {
                CudaCallCheckError( cuMemcpyHtoD(cMem.memobj, buffer_ptr, buffer_size));
            }
            else    // async write
            {
                CudaCallCheckError( cuMemcpyHtoDAsync(cMem.memobj, buffer_ptr, buffer_size, cuda_stream) );
            }
        }

        /**
         * read from CL device to buffer located in host memory
         */
        inline void enqueueReadBuffer(  CMem &cMem,             ///< memory object
                                        bool block_write,   ///< don't return until data is written
                                        size_t offset,          ///< offset to start writing from
                                        size_t buffer_size,     ///< size of buffer to write
                                        void *buffer_ptr        ///< host memory pointer
        )
        {
            if (block_write)    // sync write
            {
                CudaCallCheckError( cuMemcpyDtoH(buffer_ptr, cMem.memobj, buffer_size) );
            }
            else    // async write
            {
                CudaCallCheckError( cuMemcpyDtoHAsync(buffer_ptr, cMem.memobj, buffer_size, cuda_stream) );
            }
        }


        /**
         * enqueue nd range kernel
         */
        inline void enqueueNDRangeKernel(   CKernel &cKernel,           ///< enqueue a OpenCL kernel
											bool use_shared_mem,       			///< shared memory flag needed only for Beta kernel
                                            std::vector<void *>& kernelParams   ///< kernel input arguments
        )
        {
            //dim3 blockSize;
            //dim3 gridSize;
            unsigned int sharedMemBytes;

			// allocate shared memory only when required by kernel			
			if (use_shared_mem)
			{
				// Total amount of shared memory per block:       49152 bytes
            	sharedMemBytes = sizeof(T) * (size_t)(12 * cKernel.blockSize.x * cKernel.blockSize.y * cKernel.blockSize.z);
			}
			else
			{
				sharedMemBytes = 0;			
			}
		
//			printf("sharedMemory: %u bytes\n", sharedMemBytes);	
//			printf("blockSize: [%u, %u, %u] \n", cKernel.blockSize.x, cKernel.blockSize.y, cKernel.blockSize.z);
//			printf("gridSize: [%u, %u, %u] \n", cKernel.gridSize.x, cKernel.gridSize.y, cKernel.gridSize.z);

            CudaCallCheckError( cuLaunchKernel(cKernel.kernel,
                                                cKernel.gridSize.x,
                                                cKernel.gridSize.y,
                                                cKernel.gridSize.z,
                                                cKernel.blockSize.x,
                                                cKernel.blockSize.y,
                                                cKernel.blockSize.z,
                                                sharedMemBytes,
                                                cuda_stream,   //cuda_stream,
                                                &kernelParams[0],
                                                0        // extra     
                                                ) );
        }

        /**
         * wait until all enqueued object from command queue are finished
         */
        inline void finish()
        {
            CudaCallCheckError( cuStreamSynchronize(cuda_stream) );
        }
    };

    /**
     * \brief load information about a specific device
     */
    class CDeviceInfo : public CDevice
    {
    public:
        CUdevprop device_prop;
        int max_compute_units;      ///< maximum compute units available
        int max_work_item_dimensions;   ///< maximum number of dimensions

        size_t *max_work_item_sizes;    ///< maximum amount of work items
        int max_work_group_size;        ///< maximum group size
        int max_clock_frequency;            ///< maximum clock frequency
        int address_bits;                   ///< address bits for device
        int global_mem_cache_size;              ///< size of global memory cache
        size_t global_mem_size;                 ///< size of global memory

        int max_constant_buffer_size;           ///< maximum bytes for constant buffer
        int local_mem_size;                     ///< size of local memory
        int available;                          ///< true, if device available
        int execution_capabilities_major;             ///< kernel execution capabilities
        int execution_capabilities_minor;             ///< kernel execution capabilities
        int queue_properties;                   ///< queue properties
        
        char *name;             ///< name of device
        int driver_version;     ///< driver version of device

        /**
         * initialize device information with NULL data
         */
        inline void initCDeviceInfo()
        {
            max_work_item_sizes = NULL;

            name = NULL;
            driver_version = '\0';  //NULL;
        }

        inline CDeviceInfo()
            : CDevice()
        {
            initCDeviceInfo();
        }

        /**
         * initialize device information from existing device
         */
        inline CDeviceInfo(const CDevice &cDevice)
            : CDevice(cDevice)
        {
            initCDeviceInfo();
            loadDeviceInfo(cDevice);
        }

        inline ~CDeviceInfo()
        {
            delete[] max_work_item_sizes;
            delete[] name;

        }

        /**
         * return type of the device
         */
        inline const char* getTypeString()
        {
            return "GPU";
        }

        /**
         * load device information given by device_id
         */
        inline void loadDeviceInfo(const CDevice &cDevice)
        {
            set(cDevice.device_id); // set device id
    
            getCudaDrvErrorString( cuDeviceGetCount(&max_compute_units) );
            max_work_item_dimensions = 3;

            delete max_work_item_sizes;
            max_work_item_sizes = new size_t[max_work_item_dimensions];
            
            getCudaDrvErrorString( cuDeviceGetAttribute(&max_work_group_size, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, device_id) );
            getCudaDrvErrorString( cuDeviceGetAttribute(&max_clock_frequency, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, device_id) );
            max_clock_frequency = max_clock_frequency * 1e-3f;  // convert to MHz from kHz
            
            getCudaDrvErrorString( cuDeviceGetAttribute(&address_bits, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, device_id) );
            getCudaDrvErrorString( cuDeviceGetAttribute(&global_mem_cache_size, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, device_id) );
            getCudaDrvErrorString( cuDeviceTotalMem(&global_mem_size, device_id) );
            getCudaDrvErrorString( cuDeviceGetAttribute(&max_constant_buffer_size, CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, device_id) );
            getCudaDrvErrorString( cuDeviceGetAttribute(&local_mem_size, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, device_id) );
            getCudaDrvErrorString( cuDeviceGetAttribute(&available, CU_DEVICE_ATTRIBUTE_COMPUTE_MODE, device_id) );
            getCudaDrvErrorString( cuDeviceComputeCapability(&execution_capabilities_major, &execution_capabilities_minor, device_id) );
            getCudaDrvErrorString( cuDeviceGetAttribute(&queue_properties, CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS, device_id) );

            getCudaDrvErrorString( cuDeviceGetName(name, 256, device_id) ); 
            getCudaDrvErrorString( cuDriverGetVersion(&driver_version) );
            getCudaDrvErrorString( cuDeviceGetProperties(&device_prop, cDevice.device_id));

            // printf("Test!! Device Compute capability: %i, %i \n", execution_capabilities_major, execution_capabilities_minor);
        }

        /**
         * output previously loaded device information
         */
        inline void printDeviceInfo(std::string sLinePrefix)
        {
            std::cout << sLinePrefix << "MAX_COMPUTE_UNITS: " << max_compute_units << std::endl;
            std::cout << sLinePrefix << "MAX_WORK_ITEM_DIMENSIONS: " << max_work_item_dimensions << std::endl;

            for (int w = 0; w < max_work_item_dimensions; w++)
            {
                std::cout << sLinePrefix << "MAX_WORK_ITEM_SIZES[" << w << "]: " << max_work_item_sizes[w] << std::endl;
            }
            std::cout << sLinePrefix << "MAX_WORK_GROUP_SIZE: " << max_work_group_size << std::endl;
            std::cout << sLinePrefix << std::endl;


            std::cout << sLinePrefix << "MAX_CCLOCK_FREQUENCY: " << max_clock_frequency << std::endl;
            std::cout << sLinePrefix << "ADDRESS_BITS: " << address_bits << std::endl;
            std::cout << sLinePrefix << std::endl;

            std::cout << sLinePrefix << "GLOBAL_MEM_CACHE_SIZE: " << global_mem_cache_size << std::endl;
            std::cout << sLinePrefix << "GLOBAL_MEM_SIZE: " << global_mem_size << " (" << (global_mem_size >> 20) << "MB)" << std::endl;
            std::cout << sLinePrefix << std::endl;
            std::cout << sLinePrefix << "MAX_CONSTANT_BUFFER_SIZE: " << max_constant_buffer_size << std::endl;
            std::cout << sLinePrefix << "SHARED_MEM_PER_BLOCK: " << device_prop.sharedMemPerBlock << " bytes" << std::endl;


            std::cout << sLinePrefix << "LOCAL_MEM_SIZE: " << local_mem_size << std::endl;
            std::cout << sLinePrefix << "AVAILABLE: " << available << std::endl;
            std::cout << sLinePrefix << "COMPUTE_CAPABILITY: " << execution_capabilities_major << std::endl;
            std::cout << sLinePrefix << "CONCURRENT_KERNELS: " << queue_properties << std::endl;
            std::cout << sLinePrefix << "NAME: " << name << std::endl;
            std::cout << sLinePrefix << "DRIVER_VERSION: " << driver_version << std::endl;
        }
    };
};


#endif
