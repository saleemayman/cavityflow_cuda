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

//#include <CL/cl.h>    // remove this once all OpenCL commands are gone.
//#include <helper_cuda_drvapi.h>
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
//#include "../libtools/CProfiler.hpp"
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
        //cl_platform_id platform_id;   ///< OpenCL platform handler

        char *profile;  ///< profile string
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
            profile = NULL;
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
            delete[] profile;
            delete[] version;
            delete[] name;
            delete[] vendor;
            delete[] extensions;
        }

public:

#if 0
        /**
         * load information about platform and store to class
         */
        inline void loadPlatformInfo()
        {
            size_t size_retval;
#define loadP(flag, variable)   \
            CL_CHECK_ERROR(clGetPlatformInfo(   platform_id,    flag,   0,  NULL,   &size_retval));     \
            variable = new char[size_retval];                                           \
            CL_CHECK_ERROR(clGetPlatformInfo(   platform_id,    flag,   size_retval,    variable,   NULL)); \

            loadP(CL_PLATFORM_PROFILE, profile)
            loadP(CL_PLATFORM_VERSION, version)
            loadP(CL_PLATFORM_NAME, name)
            loadP(CL_PLATFORM_VENDOR, vendor)
            loadP(CL_PLATFORM_EXTENSIONS, extensions)
#undef load
        }
        /**
         * setup platform id to load information for platform_id
         */
        inline void load(cl_platform_id p_platform_id)
        {
            cleanupCPlatform();
            initCPlatform();
            platform_id = p_platform_id;
        }

        /**
         * load information about Platform given by platform_id
         */
        inline CPlatform(cl_platform_id p_platform_id)
        {
            initCPlatform();
            load(p_platform_id);
        }
#endif
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
        //cl_device_id  device_id;      ///< OpenCL Device ID
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
        //cl_context context;           ///< OpenCL context
        //cl_device_id *device_ids; ///< array with devices available for the context
        CUcontext context;
        CUdevice *device_ids;
        int device_ids_count;       ///< number of contexts in device_ids[]

        /**
         * load device list belonging to context cContext
         * \param cContext  existing context
         */
        //inline void load(CContext &cContext)
         inline void load()
        {
            // load device information
//          int value_size_ret;
//          CL_CHECK_ERROR(clGetContextInfo(cContext.context, CL_CONTEXT_DEVICES, 0, NULL, &value_size_ret));
            CudaCallCheckError( cuDeviceGetCount(&device_ids_count) );
            //device_ids_count = value_size_ret / sizeof(cl_device_id);
            if (device_ids_count == 0)  std::cerr << "Warning: no device found!" << std::endl;
            delete device_ids;
            device_ids = new CUdevice[device_ids_count];

            //CL_CHECK_ERROR(clGetContextInfo(cContext.context, CL_CONTEXT_DEVICES, value_size_ret, device_ids, NULL));
            /// get a handle for each device available
            for (int i = 0; i < device_ids_count; i++)
            {
                CudaCallCheckError( cuDeviceGet(&device_ids[i], i) );   
            }

            initDeviceVector();
        }


        /**
         * load device list of type 'device type' belonging to platform cPlatform
         */
#if 0        
        inline void load(   CPlatform &cPlatform,
                            cl_device_type device_type = CL_DEVICE_TYPE_ALL
        )
        {
            // read out device ids for given platform
            CL_CHECK_ERROR(clGetDeviceIDs(  cPlatform.platform_id, device_type, 0, NULL, &device_ids_count));

            delete device_ids;
            device_ids = new cl_device_id[device_ids_count];

            CL_CHECK_ERROR(clGetDeviceIDs(  cPlatform.platform_id, device_type, device_ids_count, device_ids, NULL));

            initDeviceVector();
        }
#endif
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
         * initialize devices
         */
        /*inline CDevices()
        {
            initCDevices();
        }*/

        /**
         * load device list belonging to cContext
         */
        inline CDevices()
        {
            initCDevices();
            load();
        }

        /**
         * load device list of type 'device type' belonging to platform cPlatform
         */
#if 0        
        inline CDevices(CPlatform &cPlatform, cl_device_type device_type = CL_DEVICE_TYPE_ALL)
        {
            initCDevices();
            load(cPlatform, device_type);
        }
#endif
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
            //clear();
            (*this).resize(device_ids_count);
            for (uint i = 0; i < device_ids_count; i++)
            {
                (*this)[i].set(device_ids[i]);
            }
        }
/*
        inline CDevices operator=(const CDevices& c)
        {
            return c;
        }
*/
    };



    /**
     * \brief Create and manage an OpenCL context
     */
    class CContext
    {
    public:
        CUresult error;     ///< error handler
        CUcontext context;  ///< OpenCL context id


#if 0
        /**
         * load context and device list by platform and type
         */
        inline void load(   const CPlatform &cPlatform,     ///< platform for parameters
                            cl_device_type device_type      ///< type of context and device list
        )
        {
            release();

            cl_int err_ret;
            cl_context_properties context_properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id, 0, 0};
            context = clCreateContextFromType(context_properties, device_type, NULL, NULL, &err_ret);
            CL_CHECK_ERROR(err_ret);

            retain();
        }
#endif
        /**
         * load context by platform and device
         
        inline void load(   const CPlatform &cPlatform, ///< platform for parameters
                            const CDevice &cDevice      ///< OpenCL device
        )*/
        inline void load(   const CDevice &cDevice      ///< CUDA device
        )
        {
            release();

            //cl_int err_ret;
            CUresult err_ret;
            err_ret = cuCtxCreate(&context, CU_CTX_SCHED_AUTO, cDevice.device_id);
            CudaCallCheckError(err_ret);
            //cl_context_properties context_properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id, 0, 0};
            //context = clCreateContext(context_properties, 1, &cDevice.device_id, NULL, NULL, &err_ret);
            //CL_CHECK_ERROR(err_ret);
/*
            cuCtxCreate(CUcontext *pctx, unsigned int flags, CUdevice dev);
*/
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


#if 0
#ifdef C_GL_TEXTURE_HPP

        /**
         * load context from currently running GL context
         * \param cPlatform desired CL Platform information
         */
        inline bool loadCurrentGlContext(CPlatform &cPlatform)
        {
            release();

            cl_int err_ret;

#ifdef WIN32
    #error TODO
#elifdef MACOSX
    #error TODO
#else
            GLXContext glxContext = glXGetCurrentContext();
            if (glxContext == NULL)
            {
                error << "no current glx context found" << std::endl;
                return false;
            }

            Display *display = glXGetCurrentDisplay();
            if (display == NULL)
            {
                error << "no current glx display found" << std::endl;
                return false;
            }

            cl_context_properties context_properties[] = {
                            CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id,
                            CL_GL_CONTEXT_KHR,  (cl_context_properties)glxContext,
                            CL_GLX_DISPLAY_KHR, (cl_context_properties)display,
                            NULL, NULL
                    };
#endif
            context = clCreateContextFromType(context_properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &err_ret);
//          context = clCreateContext(context_properties, 0, NULL, NULL, NULL, &err_ret);
            CL_CHECK_ERROR(err_ret);

            retain();
            return true;
        }
#endif
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

#if 0
        /**
         * initialize and load context with device_type and platform
         */
        inline CContext(    const CPlatform &cPlatform, ///< platform for parameters
                            cl_device_type device_type  ///< type of context and device list
        )
        {
            initCContext();
            load(cPlatform, device_type);
        }
#endif
        /**
         * initialize and load context with device and platform

        inline CContext(    CPlatform &cPlatform,   ///< platform for parameters
                            CDevice &cDevice        ///< OpenCL device
        )*/
        inline CContext(    CDevice &cDevice        ///< OpenCL device
        )
        {
            //load(cPlatform, cDevice);
            load(cDevice);
        }


#if 0
        /**
         * load context by platform and device list
         */
        inline CContext(    CPlatform &cPlatform,   ///< platform for parameters
                            CDevices &cDevices      ///< device list available for context
        )
        {
            initCContext();
            load(cPlatform, cDevices);
        }
#endif
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
        /*
        inline CContext(const CContext &)
        {
        }
        */
/*
        inline const CContext& operator=(const CContext &c)
        {
            return c;
        }
*/
        /**
         * increment reference counter to context to avoid deletion of context if
         * multiple classes use the same context
         */
        inline void retain()
        {
            //CL_CHECK_ERROR(clRetainContext(context));
            CudaCallCheckError( cuCtxAttach(&context, 0) );
        }

        /**
         * decrement reference counter to context if context is valid
         */
        inline void release()
        {
            if (context != NULL)
            {
                //clReleaseContext(context);
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
        /*cl_platform_id *platform_ids; ///< array with platform ids
        cl_uint platform_ids_count;     ///< number of platforms in platform array*/
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
            /*platform_ids = NULL;
            platform_ids_count = 0;*/
        }

        /**
         * load all available platforms
         */
        inline void load()
        {
            initCPlatforms();

            /*CL_CHECK_ERROR(clGetPlatformIDs(
                        0,
                        0,
                        &platform_ids_count
                    ));

            delete platform_ids;
            platform_ids = new cl_platform_id[platform_ids_count];

            CL_CHECK_ERROR(clGetPlatformIDs(
                        platform_ids_count,
                        platform_ids,
                        NULL
                    ));*/
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
            /*if (platform_ids != NULL)
                delete platform_ids;*/
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
        //cl_mem memobj;    ///< OpenCL memory object handler
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
/*      inline CMem(    CContext&       cContext,   ///< context for buffer
                        unsigned int    flags,      ///< CUDA flags
                        size_t          size,       ///< Size of memory buffer
                        void*           host_ptr    ///< Host pointer to initialize buffer
            )   */
        inline CMem(    CContext&       cContext,   ///< context for buffer
                        size_t          size        ///< Size of memory buffer
            )
        {
            memobj = '\0';  //NULL;
            //create(cContext, flags, size, host_ptr);  
            create(cContext, size); 
        }

        /**
         * release memory object (decrements OpenCL reference counter)
         */
        void release()
        {
            if (memobj != '\0')     //if (memobj != NULL)
                //clReleaseMemObject(memobj);
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
            // CudaCallCheckError( cuMemGetAddressRange(NULL, &mem_size, memobj) );
            CudaCallCheckError( cuMemGetAddressRange(&pbase, &mem_size, memobj) );
            // printf(" ->getSize() mem_size: %lu, pbase: %llu\n", mem_size, pbase);
            return mem_size;        ///< returns the size of memobj in size_t
        }

        /**
         * create new memory object
         */
        /*inline void create(
                        CContext&       cContext,   ///< context for buffer
                        unsigned int    flags,      ///< CUDA flags
                        size_t          size,       ///< Size of memory buffer
                        void*           host_ptr    ///< Host pointer to initialize buffer
            )
        {*/
        inline void create( CContext&       cContext,   ///< context for buffer
                            size_t          size        ///< Size of memory buffer
        )
        {
            CUresult errcode_ret;

            release();
            errcode_ret = cuMemAlloc(&memobj, size);
            CudaCallCheckError(errcode_ret);
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
            //program = clCreateProgramWithSource(cContext.context, count, strings, lengths, &errcode_ret);
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
/*          if (kernel != NULL)
                CL_CHECK_ERROR(clReleaseKernel(kernel));    */
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
                //CL_CHECK_ERROR(clWaitForEvents(1, &event));
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
            //cl_int param_value;
            //CL_CHECK_ERROR(clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &param_value, NULL));
            CUresult param_value;
            param_value = cuEventQuery(event);
            CudaCallCheckError(param_value);

            return param_value;
        }

        inline ~CEvent()
        {
            if (event != 0)
                //clReleaseEvent(event);
                CudaCallCheckError( cuEventDestroy(event) );
        }

        /**
         * release the event
         */
        inline void release()
        {
            //clReleaseEvent(event);
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
        CError error;       ///< error handler
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
                //CL_CHECK_ERROR(clReleaseCommandQueue(command_queue));
                CudaCallCheckError(cuStreamDestroy(cuda_stream));
            }
        }

        inline CCommandQueue()
        {
            //command_queue = 0;
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
            //cl_int errcode_ret;
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
                //CUresult cuMemcpyHtoD(CUdeviceptr dstDevice, const void *srcHost, size_t ByteCount);
            }
            else    // async write
            {
                CudaCallCheckError( cuMemcpyHtoDAsync(cMem.memobj, buffer_ptr, buffer_size, cuda_stream) );
                // CUresult cuMemcpyHtoDAsync(CUdeviceptr dstDevice, const void *srcHost, size_t ByteCount, CUstream hStream);
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
            // printf(" -> CCL::enqueueReadBuffer()\n");
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
         * set the number of grid size and threads per block for GPU
         */
        void setGridAndBlockSize(   dim3 &numBlocks, dim3 &threadsPerBlock, unsigned int work_dim,
                                    size_t &local_work_size, const size_t total_elements, size_t *global_work_size
        )
        {
            size_t totalThreadsPerBlock;
            size_t grids_x, grids_y, grids_z;

            if (work_dim == 1)
            {
                threadsPerBlock = dim3(local_work_size, 1, 1);
                totalThreadsPerBlock = threadsPerBlock.x * threadsPerBlock.y * threadsPerBlock.z;

                grids_x = (total_elements + threadsPerBlock.x - 1)/totalThreadsPerBlock;

                numBlocks = dim3(grids_x + 1, 1, 1);
            }
            else if(work_dim == 2)
            {
                threadsPerBlock = dim3(local_work_size, local_work_size, 1);
                // threadsPerBlock = dim3(1, 1, 1);
                totalThreadsPerBlock = threadsPerBlock.x * threadsPerBlock.y * threadsPerBlock.z;

                grids_x = sqrt((total_elements + threadsPerBlock.x * threadsPerBlock.y - 1) / totalThreadsPerBlock);
                // grids_x = global_work_size[0];
                // grids_y = global_work_size[1];
                grids_y = grids_x;

                numBlocks = dim3(grids_x + 1, grids_y + 1, 1);
            }
            else
            {
                threadsPerBlock = dim3(local_work_size, local_work_size, 1);

                grids_x = pow((total_elements + totalThreadsPerBlock - 1)/totalThreadsPerBlock, 1./3.);
                grids_y = grids_x;
                grids_z = grids_x;

                numBlocks = dim3(grids_x + 1, grids_y + 1, grids_z + 1);
            }
        }

        /**
         * enqueue nd range kernel
         */
        inline void enqueueNDRangeKernel(   CKernel &cKernel,           ///< enqueue a OpenCL kernel
                                            unsigned int work_dim,      ///< number of work dimensions (0, 1 or 2)
                                            int *dummy, //const CVector<3, int> &subdomain_size,        ///< number of elements in the concerned domain
                                            const size_t total_elements,        ///< total number of elements to process
                                            const size_t *global_work_offset,   ///< global work offset
                                            size_t *global_work_size,           ///< global work size
                                            size_t &local_work_size,            ///< local work size
                                            std::vector<void *>& kernelParams   ///< kernel input arguments
        )
        {
            // printf(" -> CCL::enqueueNDRangeKernel()");
#if PROFILE
          event = new cl_event();
#endif
            dim3 threadsPerBlock; //  threads per grid dimension = dim3(32, 1, 1);
            dim3 numBlocks;  //  number of blocks to launch = dim3((size + threadsPerBlock.x - 1) / threadsPerBlock.x, 1, 1);
            int sharedMemBytes;

            setGridAndBlockSize(numBlocks, threadsPerBlock, work_dim, local_work_size, total_elements, global_work_size);

            // printf("\n");
            // printf("threadsPerBlock: [%u, %u, %u] \n", threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z);
            // printf("numBlocks: [%u, %u, %u] \n", numBlocks.x, numBlocks.y, numBlocks.z);

            // sharedMemBytes = 12 * threadsPerBlock.x * threadsPerBlock.y * threadsPerBlock.z;
            // sharedMemBytes = local_work_size * 12;
            sharedMemBytes = sizeof(T) * local_work_size * 12;
            // printf("dim: %i, sharedMemBytes: %i \n", work_dim, sharedMemBytes);
            // get handle for the cuda stream
            // create();

            CudaCallCheckError( cuLaunchKernel(cKernel.kernel,
                                                numBlocks.x,
                                                numBlocks.y,
                                                numBlocks.z,
                                                threadsPerBlock.x,
                                                threadsPerBlock.y,
                                                threadsPerBlock.z,
                                                sharedMemBytes,
                                                cuda_stream,   //cuda_stream,
                                                &kernelParams[0],
                                                0        // extra     
                                                ) );
// TODO: Check if size_t is a valid argument for cuLaunchKernel() for numBlocks and threadsPerBlock variables.

#if PROFILE
            clWaitForEvents(1, event);
            char kernel_name[128];
            CL_CHECK_ERROR( clGetKernelInfo (cKernel.kernel,
                             CL_KERNEL_FUNCTION_NAME,
                             128,
                             kernel_name,
                             NULL
                             ));
            CProfilerEvent* profEvent = new CProfilerEvent(kernel_name, event);
            ProfilerSingleton::Instance()->addProfilerEvent(profEvent);
            delete event;
#endif
        }

        /**
         * wait until all enqueued object from command queue are finished
         */
        inline void finish()
        {
            //CL_CHECK_ERROR(   clFinish(command_queue));
            CudaCallCheckError( cuStreamSynchronize(cuda_stream) );
        }
    };

    /**
     * \brief load information about a specific device
     */
    class CDeviceInfo : public CDevice
    {
    public:
        //cl_device_type device_type;       ///< OpenCL device type
        //cl_uint vendor_id;                ///< OpenCL vendor id
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
//          vendor = NULL;
            driver_version = '\0';  //NULL;
//          profile = NULL;
//          version = NULL;
//          extensions = NULL;
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
//          delete[] vendor;
//          delete[] driver_version;
//          delete[] profile;
//          delete[] version;
//          delete[] extensions;

        }

        /**
         * return type of the device
         */
        inline const char* getTypeString()
        {
            return "GPU";
// CUDA Driver API only deals with device type GPU          
#if 0
            switch(device_type)
            {
                case CL_DEVICE_TYPE_CPU:    return "CPU";
                case CL_DEVICE_TYPE_GPU:    return "GPU";
                case CL_DEVICE_TYPE_ACCELERATOR:    return "ACCELERATOR";
                case CL_DEVICE_TYPE_DEFAULT:    return "DEFAULT";
                case CL_DEVICE_TYPE_ALL:    return "ALL";
                default:            return "unknown";
            }
#endif
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
