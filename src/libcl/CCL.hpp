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

//#include <CL/cl.h>	// remove this once all OpenCL commands are gone.
#include </usr/local/cuda/include/cuda.h>	// check if it actually reads the definitions from cuda.h
//#include <helper_cuda_drvapi.h>
#include <iostream>
#include <strings.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <sstream>
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
		#define		CL_GL_CONTEXT_KHR		0x2008
		#define		CL_EGL_DISPLAY_KHR		0x2009
		#define		CL_GLX_DISPLAY_KHR		0x200A
		#define		CL_WGL_HDC_KHR			0x200B
		#define		CL_CGL_SHAREGROUP_KHR	0x200C
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
		//cl_platform_id platform_id;	///< OpenCL platform handler

		char *profile;	///< profile string
		char *version;	///< version string
		char *name;		///< name string
		char *vendor;	///< vendor string
		char *extensions;	///< extensions string

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
#define loadP(flag, variable)	\
			CL_CHECK_ERROR(clGetPlatformInfo(	platform_id,	flag,	0,	NULL,	&size_retval));		\
			variable = new char[size_retval];											\
			CL_CHECK_ERROR(clGetPlatformInfo(	platform_id,	flag,	size_retval,	variable,	NULL));	\

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
		//cl_device_id	device_id;		///< OpenCL Device ID
		CUdevice device_id;				///< CUDA Device ID

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
	class CDevices	: public std::vector<CDevice>
	{
	public:
		//cl_context context;			///< OpenCL context
		//cl_device_id *device_ids;	///< array with devices available for the context
		CUcontext context;
		CUdevice *device_ids;
		int device_ids_count;		///< number of contexts in device_ids[]

		/**
		 * load device list belonging to context cContext
		 * \param cContext	existing context
		 */
		//inline void load(CContext &cContext)
		 inline void load()
		{
			// load device information
//			int value_size_ret;
//			CL_CHECK_ERROR(clGetContextInfo(cContext.context, CL_CONTEXT_DEVICES, 0, NULL, &value_size_ret));
			CudaCallCheckError( cuDeviceGetCount(&device_ids_count) );
			//device_ids_count = value_size_ret / sizeof(cl_device_id);
			if (device_ids_count == 0)	std::cerr << "Warning: no device found!" << std::endl;
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
		inline void load(	CPlatform &cPlatform,
							cl_device_type device_type = CL_DEVICE_TYPE_ALL
		)
		{
			// read out device ids for given platform
			CL_CHECK_ERROR(clGetDeviceIDs(	cPlatform.platform_id, device_type, 0, NULL, &device_ids_count));

			delete device_ids;
			device_ids = new cl_device_id[device_ids_count];

			CL_CHECK_ERROR(clGetDeviceIDs(	cPlatform.platform_id, device_type, device_ids_count, device_ids, NULL));

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
		CUresult error;		///< error handler
		CUcontext context;	///< OpenCL context id


#if 0
		/**
		 * load context and device list by platform and type
		 */
		inline void load(	const CPlatform &cPlatform,		///< platform for parameters
							cl_device_type device_type		///< type of context and device list
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
		 
		inline void load(	const CPlatform &cPlatform,	///< platform for parameters
							const CDevice &cDevice		///< OpenCL device
		)*/
		inline void load(	const CDevice &cDevice		///< CUDA device
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
		inline void load(	const CPlatform &cPlatform,	///< platform for parameters
							const CDevices &cDevices	///< device list available for context
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
		 * \param cPlatform	desired CL Platform information
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
							CL_GL_CONTEXT_KHR,	(cl_context_properties)glxContext,
							CL_GLX_DISPLAY_KHR,	(cl_context_properties)display,
							NULL, NULL
					};
#endif
			context = clCreateContextFromType(context_properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &err_ret);
//			context = clCreateContext(context_properties, 0, NULL, NULL, NULL, &err_ret);
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
		inline CContext(	const CPlatform &cPlatform,	///< platform for parameters
							cl_device_type device_type	///< type of context and device list
		)
		{
			initCContext();
			load(cPlatform, device_type);
		}
#endif
		/**
		 * initialize and load context with device and platform

		inline CContext(	CPlatform &cPlatform,	///< platform for parameters
							CDevice &cDevice		///< OpenCL device
		)*/
		inline CContext(	CDevice &cDevice		///< OpenCL device
		)
		{
			//load(cPlatform, cDevice);
			load(cDevice);
		}


#if 0
		/**
		 * load context by platform and device list
		 */
		inline CContext(	CPlatform &cPlatform,	///< platform for parameters
							CDevices &cDevices		///< device list available for context
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
	 * \brief	handle multiple platforms available for computations
	 */
	class CPlatforms
	{
	public:
		/*cl_platform_id *platform_ids;	///< array with platform ids
		cl_uint platform_ids_count;		///< number of platforms in platform array*/
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
		//cl_mem memobj;	///< OpenCL memory object handler
		CUdeviceptr memobj;

		inline CMem()
		{
			memobj = '\0';	//NULL;
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
/*		inline CMem(	CContext&		cContext,	///< context for buffer
						unsigned int	flags,		///< CUDA flags
						size_t			size,		///< Size of memory buffer
						void*			host_ptr	///< Host pointer to initialize buffer
			)	*/
		inline CMem(	CContext&		cContext,	///< context for buffer
						size_t			size		///< Size of memory buffer
			)
		{
			memobj = '\0';	//NULL;
			//create(cContext, flags, size, host_ptr);	
			create(cContext, size);	
		}

		/**
		 * release memory object (decrements OpenCL reference counter)
		 */
		void release()
		{
			if (memobj != '\0')		//if (memobj != NULL)
				//clReleaseMemObject(memobj);
				CudaCallCheckError( cuMemFree(memobj) );
			memobj = '\0';	//NULL;
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
			size_t mem_size;
			CudaCallCheckError( cuMemGetAddressRange(NULL, &mem_size, memobj) );
			//clGetMemObjectInfo(memobj, CL_MEM_SIZE, sizeof(size_t), &mem_size, NULL);
			return mem_size;		///< returns the size of memobj in size_t
		}

#if 0
		/**
		 * create OpenCL 2D image memory object
		 */
		// TODO: Replace this with CUDA texture reference
		inline void createImage2D(
									CContext&		cContext,	///< context for buffer
									cl_mem_flags	flags,		///< OpenCL flags
									const cl_image_format	*image_format,	///< image format
									size_t	image_width,		///< width of image
									size_t	image_height,		///< height of image
									size_t	image_row_pitch = 0,	///< row pitch to improve aligned access
									void 	*host_ptr = NULL	///< host pointer to initalize with existing data
			)
		{
			release();

			cl_int errcode_ret;
			memobj = clCreateImage2D(	cContext.context, flags, image_format,
										image_width, image_height, image_row_pitch,
										host_ptr, &errcode_ret);

			CL_CHECK_ERROR(errcode_ret);
		}
#endif
		/**
		 * create new memory object
		 */
		/*inline void create(
						CContext&		cContext,	///< context for buffer
						unsigned int	flags,		///< CUDA flags
						size_t			size,		///< Size of memory buffer
						void*			host_ptr	///< Host pointer to initialize buffer
			)
		{*/
		inline void create(	CContext&		cContext,	///< context for buffer
							size_t			size		///< Size of memory buffer
		)
		{
			release();
			/*
			cl_int errcode_ret;
			memobj = clCreateBuffer(cContext.context, flags, size, host_ptr, &errcode_ret);
			CL_CHECK_ERROR(errcode_ret);
			*/

			CUresult errcode_ret;

			errcode_ret = cuMemAlloc(&memobj, size);
			CudaCallCheckError(errcode_ret);

		}

		/*
		 * same as create() but copies the data from host to device
		 * OpenCL equivalent is clCreateBuffer with  CL_MEM_COPY_HOST_PTR  flag
		 */
		inline void createCopyToDevice(	CContext&		cContext,	///< context for buffer
										size_t			size,		///< Size of memory buffer
										void*			host_ptr	///< Host pointer to initialize buffer
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
				CContext&		cContext,	///< context for buffer
				cl_mem_flags	flags,		///< OpenCL flags
				CGlTexture 		&texture	///< reference to 3d texture
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
				CContext&		cContext,	///< context for buffer
				cl_mem_flags	flags,		///< OpenCL flags
				CGlTexture 		&texture	///< reference to 3d texture
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
		std::string command;	///< command string to execute
		std::stringstream command_string_stream;	///< string to append
		const char *file_name;	///< cuda file to compile

public:
		/*CError error;			///< error handler
		cl_program program;		///< OpenCL program id*/
		CError error;			///< error handler
		CUmodule program;		///< CUDA program id equivalent
		std::string filepath;	///< string with filename (for debugging purposes)

		inline CProgram()
		{
			program = 0;
		}

		/**
		 * load program given by p_filepath
		 */
		inline CProgram(	CContext &cContext,	// context
					const std::string &p_filepath	// source file
		)
		{
			CProgram();
			load(cContext, p_filepath.c_str() );
		}

#if 0
		/**
		 * load program given by p_filepath and place prefix_string before the filecontent
		 */
		inline CProgram(	CContext &cContext,	// context
					const std::string &p_filepath,	// source file
					const std::string &prefix_string	// prefix placed before source file
		)
		{
			CProgram();
			load(cContext, p_filepath, prefix_string);
		}

		/**
		 * initialize program with source code
		 */
		inline CProgram(	CContext &cContext,		// context
							uint count,			// number of strings
							const char **strings,	// source code
							const size_t *lengths	// length of source code strings. if NULL, the strings are \0 terminated
		)
		{
			CProgram();
			load(cContext, count, strings, lengths);
		}
#endif
		/**
		 * load kernel using source given in strings
		 */
/*		inline void load(	CContext &cContext,	// context
							uint count,		// number of strings
							const char **strings,	// source code
							const size_t *lengths	// length of source code strings. if NULL, the strings are \0 terminated
		)	*/
		inline void load(	CContext &cContext,	// context
							const char *cudaKernelFileName	// cuda cubin/ptx source code
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
		 * load kernel source from file and prepent prefix_string
		 */
		inline void load(	CContext &cContext,					// context
							const std::string &sPrefixString,	// prefix placed before source file
							const std::string &p_filepath		// source file
		)
		{
			filepath = p_filepath;

			std::string fileContent;
			std::string errorLog;
#if 1
			std::string source = sPrefixString;
			source += "\n";
			source += "#include \"";
			source += filepath;
			source += "\"";

			const char* strings[] = {source.c_str()};
			size_t lengths[] = {source.size()};
			//load(cContext, 1, strings, lengths);
			load(cContext, p_filepath.c_str() );		// for cuda module load
#else
			if (!CFile::fileContents(filepath, fileContent, errorLog))
			{
				error << "CProgram: Failed to open file '" << filepath << "' - " << strerror(errno) << std::endl;
				return;
			}

			const char* strings[] = {sPrefixString.c_str(), fileContent.c_str()};
			size_t lengths[] = {sPrefixString.size(), fileContent.size()};
			load(cContext, 2, strings, lengths);
#endif
		}

#if 0
		/**
		 * load kernel source from file
		 */
		inline void load(	CContext &cContext,				// context
							const std::string &p_filepath	// source file
		)
		{
			filepath = p_filepath;

			std::string fileContent;
			std::string errorLog;

			CFile cFile;
			cFile.fileContents(filepath, fileContent, errorLog);

			if (cFile.error())
			{
				std::cerr << "CProgram: Failed to open file '" << filepath << "' - " << cFile.error.getString() << std::endl;
				return;
			}

			const char* strings[] = {fileContent.c_str()};
			size_t lengths[] = {fileContent.size()};
			load(cContext, 1, strings, lengths);
		}

		/**
		 * return build information string
		 */
		std::string getBuildInfo(	CDevice &device		///< device on which the code was build for
								)
		{
			char *build_log;
			size_t ret_val_size;
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
			build_log = new char[ret_val_size+1];
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);

			// to be carefully, terminate with \0
			// there's no information in the reference wheter the string is 0 terminated or not
			build_log[ret_val_size] = '\0';

			std::string ret_string = build_log;
			delete build_log;
			return ret_string;
		}
#endif
		/**
		 * execute compile command on system
		 */
		inline void executeCommand(	const std::string &command, 
									const char *file_name
		)
		{
			int err_ret;
			std::cout << "Executing Compile command --> " << command << std::endl;
			err_ret = system(command.c_str());
			std::cout << "Compile Successful, file --> " << file_name << std::endl;

			if (err_ret != 0)
			{
				std::cerr << "failed to compile " << file_name << std::endl;
				exit(err_ret);
			}
		}
#if 0
		/**
		 * output OpenCL binary code
		 */
		void printBinaries()
		{
			cl_uint program_num_devices;
			CL_CHECK_ERROR(clGetProgramInfo(	program,
								CL_PROGRAM_NUM_DEVICES,
								sizeof(cl_uint),
								&program_num_devices,
								NULL
					));

			if (program_num_devices == 0)
			{
				std::cerr << "no valid binary was found" << std::endl;
				return;
			}

			size_t *binaries_sizes = new size_t[program_num_devices];

			CL_CHECK_ERROR(clGetProgramInfo(	program,
								CL_PROGRAM_BINARY_SIZES,
								program_num_devices*sizeof(size_t),
								binaries_sizes,
								NULL
					));

			char **binaries = new char*[program_num_devices];

			for (size_t i = 0; i < program_num_devices; i++)
				binaries[i] = new char[binaries_sizes[i]+1];

			CL_CHECK_ERROR(clGetProgramInfo(program, CL_PROGRAM_BINARIES, program_num_devices*sizeof(size_t), binaries, NULL));

			for (size_t i = 0; i < program_num_devices; i++)
			{
				binaries[i][binaries_sizes[i]] = '\0';
				std::cout << "Program " << i << ":" << std::endl;
				std::cout << binaries[i];
			}

			for (size_t i = 0; i < program_num_devices; i++)
				delete [] binaries[i];

			delete [] binaries;
			delete [] binaries_sizes;
		}
#endif
	};

/*
TODO: change all setArg(...) to include the variable index.
*/

	/**
	 * \brief Create kernels from CProgram, use and set arguments for them
	 */
	class CKernel
	{
	public:
/*		cl_kernel kernel;	///< OpenCL kernel handler
		CError error;		///< error handler 	*/
		CUfunction kernel;	///< CUDA funtion handle
		CUresult error;		///< cuda error handler
		std::vector<void *> kernelArgsVec;	///< vector to hold list of kernel input parameters
		//void **kernelArgs;	///< array to send to culaunchKernel

		/**
		 * create kernel from OpenCL program
		 */
		inline void create(	CProgram &program,			///< OpenCL program
							const char *kernel_name		///< name of kernel function to create kernel from
			)
		{
			// check how to include this. Now it gives a compiler error
			//if (error())
			//	std::cout << "FUCK" << std::endl;

			/*cl_int errcode_ret;
			kernel = clCreateKernel(program.program, kernel_name, &errcode_ret);
			if (errcode_ret != CL_SUCCESS)
				error << cclGetErrorString(errcode_ret) << std::endl;	*/

			CudaCallCheckError( cuModuleGetFunction(&kernel, program.program, kernel_name) );
		}

		/**
		 * create kernel from OpenCL program
		 */
		inline CKernel(	CProgram &program,			///< OpenCL program
						const char *kernel_name	///< name of kernel function to create kernel from
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
/*			if (kernel != NULL)
				CL_CHECK_ERROR(clReleaseKernel(kernel));	*/
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
/*			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(cl_mem), &(cmem.memobj));

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;	*/
			//kernelArgsVec.push_back(cmem);
			kernelArgsVec[arg_index] = &cmem;
		}

		/**
		 * set cl_float kernel argument
		 */
		inline void setArg(	int arg_index,
							float &arg		///< reference to float argument
		)
		{
/*			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(float), &arg);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;	*/
			kernelArgsVec[arg_index] = &arg;
		}

		/**
		 * set cl_double kernel argument
		 */
		inline void setArg(	int arg_index,
							double &arg		///< reference to float argument
		)
		{
/*			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(double), &arg);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;	*/
			kernelArgsVec[arg_index] = &arg;
		}

		/**
		 * set cl_int kernel argument
		 */
		inline void setArg(	int arg_index,
							int &arg
		)
		{
/*			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(cl_int), &size);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;	*/
			kernelArgsVec[arg_index] = &arg;
		}

		inline void setArg(	int arg_index,
							const int &arg
		)
		{
/*			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(cl_int), &size);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;	*/
			kernelArgsVec[arg_index] = const_cast<int *>(&arg);	//const_cast < new_type > ( expression ) 		
		}

// this member is not called in CLbmSolver so I am commenting it out for now
#if 0
		/**
		 * return the maximum work group size
		 * \param cDevice	device where the kernel should be executed
		 */
		inline size_t getMaxWorkGroupSize(CDevice &cDevice)
		{
			size_t param_value;
			clGetKernelWorkGroupInfo(kernel, cDevice.device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &param_value, NULL);
			return param_value;
		}
#endif
	};
/**
 * 00:07 Uhr 05.11.2014
 * CUDA Drv API only requires the cubin/ptx filename of the kernel to load
 * the kernel; OpenCL clCreateProgram... uses all the different options
 * given in the **strings argument(as displayed when the exec was run)
 * TODO:
 * 	- load the kernel using CUDA Drv API cuModuleLoad (Done: load(CContext, cudaKernelFileName))
 */

	/**
	 * \brief OpenCL event handler
	 */
	class CEvent
	{
	public:
		//cl_event event;	///< openCL event
		CUevent event;		///< CUDA event

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

#if 0
	/**
	 * \brief OpenCL event list management (up to 10 events)
	 */
	class CEvents
	{
	public:
		cl_event events[10];	///< events (max: 10)
		cl_uint events_count;	///< number of events in events list

		/**
		 * initialize empty event list
		 */
		inline CEvents()
		{
			events_count = 0;
		}

		/**
		 * initialize with 1 event
		 */
		inline CEvents(CEvent &event0)
		{
			set(event0);
		}

		/**
		 * initialize with 2 events
		 */
		inline CEvents(CEvent &event0, CEvent &event1)
		{
			set(event0, event1);
		}

		/**
		 * initialize with 1 event
		 */
		inline void set(CEvent &event0)
		{
			events_count = 1;
			events[0] = event0.event;
		}

		/**
		 * initialize with 2 events
		 */
		inline void set(CEvent &event0, CEvent &event1)
		{
			events_count = 2;
			events[0] = event0.event;
			events[1] = event1.event;
		}
	};
#endif

	/**
	 * \brief OpenCL command queue handler
	 */
	class CCommandQueue
	{
	public:
		CError error;		///< error handler

		//cl_command_queue command_queue;	///< OpenCL command queue handler
		CUstream cuda_stream;	///< CUDA stream handler

// Ignore code not required for LBM code execution to get a simple example
// working with no compiler or linking bugs.		
#if 0
		/**
		 * increment OpenCL reference counter to command queue
		 */
		inline void retain()
		{
			CL_CHECK_ERROR(clRetainCommandQueue(command_queue));
		}
#endif
		/**
		 * decrement OpenCL reference counter to command queue
		 */
		inline void release()
		{
			if (cuda_stream != 0)	// check if this is correct
			{
				//CL_CHECK_ERROR(clReleaseCommandQueue(command_queue));
				CudaCallCheckError(cuStreamDestroy(cuda_stream));
			}
		}

// Ignore code not required for LBM code execution to get a simple example
// working with no compile or linking bugs.		
#if 0
		/**
		 * initialize existing command queue
		 */
		inline void init(	const CContext &cContext,	///< existing context handler
							const CDevice &cDevice		///< existing device handler
		)
		{
			cl_command_queue_properties properties = 0;
#if PROFILE
			properties |= CL_QUEUE_PROFILING_ENABLE;
			CCL::CDeviceInfo cDeviceInfo(cDevice);
			printf("device timer resolution: %zu nanoseconds\n", cDeviceInfo.profiling_timer_resolution);
#endif
			cl_int errcode_ret;
			command_queue = clCreateCommandQueue(
									cContext.context,
									cDevice.device_id,
									properties,
									&errcode_ret
							);

			CL_CHECK_ERROR(errcode_ret);
			retain();
		}

		/**
		 * initialize from existing command queue and increment reference counter
		 */
		void initWithExistingCommandQueue(const CCommandQueue &cCommandQueue)
		{
			command_queue = cCommandQueue.command_queue;
			retain();
		}

		/**
		 * initialize command queue from existing context and device
		 */
		inline CCommandQueue(	const CContext &cContext,	///< existing context handler
								const CDevice &cDevice		///< existing device handler
		)
		{
			init(cContext, cDevice);
		}
#endif
		inline CCommandQueue()
		{
			//command_queue = 0;
			cuda_stream = 0;
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
			CUresult errcode_ret;	// cuda error type
			errcode_ret = cuStreamSynchronize(cuda_stream);
			if (errcode_ret != CUDA_SUCCESS)
				error << getCudaDrvErrorString(errcode_ret) << std::endl;	// CUDA Driver API erro handling

		}

		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr		///< host memory pointer
		)
		{
			if (block_write)	// sync write
			{
				CudaCallCheckError( cuMemcpyDtoH(buffer_ptr, cMem.memobj, buffer_size) );
			}
			else 	// async write
			{
				CudaCallCheckError( cuMemcpyDtoHAsync(buffer_ptr, cMem.memobj, buffer_size, cuda_stream) );
			}
		}

// Ignore code not required for LBM code (BASIC) execution to get a simple example
// working with no compile or linking bugs.		
#if 0
		/**
		 * copy data from buffer in host memory to CL device
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,
						cl_bool block_write,
						size_t offset,
						size_t buffer_size,
						const void *buffer_ptr
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
		}

		/**
		 * enqueue writing a buffer
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										const void *buffer_ptr,	///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list		///< list of events
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								NULL
					));
		}

		/**
		 * enqueue writing a buffer
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										const void *buffer_ptr,	///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list,	///< list of events
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								&event.event
					));
		}

		/**
		 * enqueue writing a buffer
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										const void *buffer_ptr,	///< host memory pointer
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								&event.event
					));
		}

		/**
		 * copy data from buffer in host memory to CL device
		 */
		inline void enqueueWriteBufferRect(	CMem &cMem,
						cl_bool block_write,
						const size_t buffer_origin[3],
						const size_t host_origin[3],
						const size_t region[3],
						const void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueWriteBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								0,
								0,
								0,
								0,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif
		}
		inline void enqueueWriteBufferRect(	CMem &cMem,
						cl_bool block_write,
						const size_t buffer_origin[3],
						const size_t host_origin[3],
						const size_t region[3],
						const size_t buffer_row_pitch,
						const size_t buffer_slice_pitch,
						const size_t host_row_pitch,
						const size_t host_slice_pitch,
						const void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueWriteBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								buffer_row_pitch,
								buffer_slice_pitch,
								host_row_pitch,
								host_slice_pitch,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif
		}
		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr		///< host memory pointer
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
		}


		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr,		///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list		///< list of events
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								NULL
					));
		}


		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr,		///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list,	///< list of events
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								&event.event
					));
		}


		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr,		///< host memory pointer
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								&event.event
					));
		}

		inline void enqueueReadBufferRect(	CMem &cMem,
				cl_bool block_write,
				const size_t buffer_origin[3],
				const size_t host_origin[3],
				const size_t region[3],
				void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueReadBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								0,
								0,
								0,
								0,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif

		}

		inline void enqueueReadBufferRect(	CMem &cMem,
				cl_bool block_write,
				const size_t buffer_origin[3],
				const size_t host_origin[3],
				const size_t region[3],
				const size_t buffer_row_pitch,
				const size_t buffer_slice_pitch,
				const size_t host_row_pitch,
				const size_t host_slice_pitch,
				void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueReadBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								buffer_row_pitch,
								buffer_slice_pitch,
								host_row_pitch0,
								host_slice_pitch,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif

		}

		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBuffer(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										size_t src_offset,	///< offset in source memory
										size_t dst_offset,	///< offset in destination memory
										size_t cb			///< number of bytes to be copied
		)
		{
			CL_CHECK_ERROR(	clEnqueueCopyBuffer(
								command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_offset,
								dst_offset,
								cb,
								0,
								NULL,
								NULL
					));
		}

		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBuffer(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										size_t src_offset,	///< offset in source memory
										size_t dst_offset,	///< offset in destination memory
										size_t cb,			///< number of bytes to be copied
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list,	///< list of events
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueCopyBuffer(
								command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_offset,
								dst_offset,
								cb,
								num_events_in_wait_list,
								event_wait_list,
								&event.event
					));
		}

		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBufferRect(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										const size_t src_origin[3],
										const size_t dst_origin[3],
										const size_t region[3]
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueCopyBufferRect(	command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_origin,
								dst_origin,
								region,
								0,
								0,
								0,
								0,
								0,
								NULL,
								NULL
					));
#endif
		}

		inline void enqueueCopyBufferRect(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										const size_t src_origin[3],
										const size_t dst_origin[3],
										const size_t region[3],
										const size_t src_row_pitch,
										const size_t src_slice_pitch,
										const size_t dst_row_pitch,
										const size_t dst_slice_pitch
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueCopyBufferRect(	command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_origin,
								dst_origin,
								region,
								src_row_pitch,
								src_slice_pitch,
								dst_row_pitch,
								dst_slice_pitch,
								0,
								NULL,
								NULL
					));
#endif
		}
		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBufferToImage(	CMem &cSrcMem,		///< source memory object
												CMem &cDstMem,		///< destination image object
												size_t src_offset,	///< offset in source memory
												const size_t dst_origin[3],	///< coordinates in source image
												const size_t region[3]		///< area to be copied
		)
		{
			CL_CHECK_ERROR(	clEnqueueCopyBufferToImage(
								command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_offset,
								dst_origin,
								region,
								0,
								NULL,
								NULL
					));
		}
#endif
		/**
		 * set the number of grid size and threads per block for GPU
		 */
		void setGridAndBlockSize(	dim3 &grid, dim3 &block, unsigned int work_dim,
									const int grid_size_x, const size_t total_elems, 
									size_t *local_work_size, size_t *global_work_size
		)
		{
			size_t threads_per_block;

			// if no work-group size specified set it to default defined in lbm_defaults.h
			// if (local_work_size == NULL)
			// {
			// 	if (work_dim == 1)
			// 	{
			// 		block = dim3(LOCAL_WORK_GROUP_SIZE, 1 ,1);
			// 		global_work_size[0] = total_elems;
					
			// 		grid = dim3((total_elems + block.x - 1)/block.x, 1, 1);
			// 	}
			// 	else if(work_dim == 2)
			// 	{
			// 		block = dim3(LOCAL_WORK_GROUP_SIZE/2, 2, 1);

			// 		global_work_size[0] = grid_size_x * LOCAL_WORK_GROUP_SIZE;
			// 		global_work_size[1] = total_elems/global_work_size[0];

			// 		grid = dim3(grid_size_x, 
			// 					(total_elems + global_work_size[0] - 1)/global_work_size[0], 1);
			// 	}
			// 	else
			// 	{
			// 		block = dim3(LOCAL_WORK_GROUP_SIZE/4, 2 ,2);
					
			// 		global_work_size[0] = grid_size_x * LOCAL_WORK_GROUP_SIZE;
			// 		global_work_size[1] = total_elems/global_work_size[0];
			// 		global_work_size[2] = 1;

			// 		grid = dim3(grid_size_x, 
			// 					(total_elems + global_work_size[0] - 1)/global_work_size[0], 1);
			// 	}
			// }
			// else

			block = dim3(local_work_size[0], local_work_size[1], local_work_size[2]);
			threads_per_block = block.x * block.y * block.z;

			if (work_dim == 1)
			{
				global_work_size[0] = total_elems;
				grid = dim3((total_elems + threads_per_block - 1)/threads_per_block, 1, 1);
			}
			else if(work_dim == 2)
			{
				global_work_size[0] = grid_size_x * threads_per_block;
				global_work_size[1] = total_elems/global_work_size[0];

				grid = dim3(grid_size_x, 
							(total_elems + global_work_size[0] - 1)/global_work_size[0], 1);
			}
			else
			{
				global_work_size[0] = grid_size_x * threads_per_block;
				global_work_size[1] = total_elems/global_work_size[0];
				global_work_size[2] = 1;

				grid = dim3(grid_size_x, 
							(total_elems + global_work_size[0] - 1)/global_work_size[0], 1);
			}
		}

		/**
		 * enqueue nd range kernel
		 */
/*		inline void enqueueNDRangeKernel( 	CKernel &cKernel, ///< enqueue a OpenCL kernel
											cl_uint work_dim, ///< number of work dimensions (0, 1 or 2)
											const size_t *global_work_offset, ///< global work offset
											const size_t *global_work_size, ///< global work size
											const size_t *local_work_size ///< local work size
											)	*/		 
		inline void enqueueNDRangeKernel(	CKernel &cKernel,					///< enqueue a OpenCL kernel
											unsigned int work_dim,				///< number of work dimensions (0, 1 or 2)
											const unsigned int grid_size_x,		///< number of blocks in x-direction (e.g., total_elems/threads_per_block[0])
											const size_t total_elements,		///< total number of elements to process
											const size_t *global_work_offset,	///< global work offset
											size_t *global_work_size,			///< global work size
											size_t *local_work_size,			///< local work size
											std::vector<void *>& kernelParams	///< kernel input arguments
		)
		{
		  //cl_event* event = NULL;
			//CUevent *event = NULL;	// cuda event type
#if PROFILE
		  event = new cl_event();
#endif
/*			CL_CHECK_ERROR(	clEnqueueNDRangeKernel(	command_queue,
								cKernel.kernel,
								work_dim,
								global_work_offset,
								global_work_size,
								local_work_size,
								0,
								NULL,
								event
					));		*/

			dim3 block; //	threads per grid dimension = dim3(32, 1, 1);
			dim3 grid;	// 	number of blocks to launch = dim3((size + block.x - 1) / block.x, 1, 1);

			setGridAndBlockSize(grid, block, work_dim, grid_size_x, total_elements, local_work_size, global_work_size);

			CudaCallCheckError( cuLaunchKernel(cKernel.kernel,
												grid.x,
												grid.y,
												grid.z,
												block.x,
												block.y,
												block.z,
												0,		// sharedMemBytes
												cuda_stream,
												&kernelParams[0],
												0 		// extra	 
												) );
// TODO: Check if size_t is a valid argument for cuLaunchKernel() for grid and block variables.

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

// Ignore code not required for LBM code execution to get a simple example
// working with no compile or linking bugs.		
#if 0
		/**
		 * enqueue nd range kernel
		 */
		inline void enqueueNDRangeKernel(	CKernel &cKernel,		///< enqueue a OpenCL kernel
											cl_uint work_dim,		///< number of work dimensions (0, 1 or 2)
											const size_t *global_work_offset,	///< global work offset
											const size_t *global_work_size,		///< global work size
											const size_t *local_work_size,		///< local work size
											CEvents &events			///< events to wait for
		)
		{
			CL_CHECK_ERROR(	clEnqueueNDRangeKernel(	command_queue,
								cKernel.kernel,
								work_dim,
								global_work_offset,
								global_work_size,
								local_work_size,
								events.events_count,
								events.events,
								NULL
					));
		}

		/**
		 * enqueue nd range kernel
		 */
		inline void enqueueNDRangeKernel(	CKernel &cKernel,		///< enqueue a OpenCL kernel
											cl_uint work_dim,		///< number of work dimensions (0, 1 or 2)
											const size_t *global_work_offset,	///< global work offset
											const size_t *global_work_size,		///< global work size
											const size_t *local_work_size,		///< local work size
											CEvent &event			///< event to wait for
		)
		{
			CL_CHECK_ERROR(	clEnqueueNDRangeKernel(	command_queue,
								cKernel.kernel,
								work_dim,
								global_work_offset,
								global_work_size,
								local_work_size,
								0,
								NULL,
								&event.event
					));
		}
#endif
		/**
		 * wait until all enqueued object from command queue are finished
		 */
		inline void finish()
		{
			//CL_CHECK_ERROR(	clFinish(command_queue));
			CudaCallCheckError( cuStreamSynchronize(cuda_stream) );
		}

#ifdef C_GL_TEXTURE_HPP
		/**
		 * acquire memory object from OpenGL context
		 */
		inline void enqueueAcquireGLObject(CMem &cMem)
		{
			CL_CHECK_ERROR(	clEnqueueAcquireGLObjects(command_queue, 1, &(cMem.memobj), 0, NULL, NULL));
		}

		/**
		 * release memory object from OpenGL context
		 */
		inline void enqueueReleaseGLObject(CMem &cMem)
		{
			CL_CHECK_ERROR(	clEnqueueReleaseGLObjects(command_queue, 1, &(cMem.memobj), 0, NULL, NULL));
		}
#endif
	};

// No CUDA platform concept. Just change it to display device info
#if 0
	/**
	 * output informations about all profiles and devices
	 */
	inline void printPlatformInfo(cl_device_type device_type = CL_DEVICE_TYPE_ALL)
	{
		CPlatforms platforms;
		std::cout << "platforms_ids: " << platforms.platform_ids_count << std::endl;

		for (cl_uint i = 0; i < platforms.platform_ids_count; i++)
		{
			CPlatform platform(platforms.platform_ids[i]);
			platform.loadPlatformInfo();

			std::cout << "  [" << i << "] PLATFORM_PROFILE: " << platform.profile << std::endl;
			std::cout << "  [" << i << "] PLATFORM_VERSION: " << platform.version << std::endl;

			CDevices cDevices(platform, device_type);

			for (cl_uint d = 0; d < cDevices.device_ids_count; d++)
			{
				CDeviceInfo deviceInfo(cDevices[d]);

				std::ostringstream sPrefix;
			       	sPrefix << "    [" << d << "] ";

				deviceInfo.printDeviceInfo(sPrefix.str());
			}
		}
	}
#endif

	/**
	 * \brief load information about a specific device
	 */
	class CDeviceInfo : public CDevice
	{
	public:
		//cl_device_type device_type;		///< OpenCL device type
		//cl_uint vendor_id;				///< OpenCL vendor id
		int max_compute_units;		///< maximum compute units available
		int max_work_item_dimensions;	///< maximum number of dimensions

		size_t *max_work_item_sizes;	///< maximum amount of work items
		int max_work_group_size;		///< maximum group size
		int max_clock_frequency;			///< maximum clock frequency
		int address_bits;					///< address bits for device
		int global_mem_cache_size;				///< size of global memory cache
		size_t global_mem_size;					///< size of global memory

		int max_constant_buffer_size;			///< maximum bytes for constant buffer
		int local_mem_size;						///< size of local memory
		int available;							///< true, if device available
		int execution_capabilities;				///< kernel execution capabilities
		int queue_properties;					///< queue properties
		
		char *name;				///< name of device
		int driver_version;		///< driver version of device

#if 0
		cl_uint preferred_vector_width_char;	///< preferred vector width for type char
		cl_uint preferred_vector_width_short;	///< preferred vector width for type short
		cl_uint preferred_vector_width_int;		///< preferred vector width for type int
		cl_uint preferred_vector_width_long;	///< preferred vector width for type long
		cl_uint preferred_vector_width_float;	///< preferred vector width for type float
		cl_uint preferred_vector_width_double;	///< preferred vector width for type float
		cl_ulong max_mem_alloc_size;			///< maximum number of allocatable bytes

		cl_bool image_support;					///< image support available
		cl_uint max_read_image_args;			///< maximum number of images as read kernel arguments
		cl_uint max_write_image_args;			///< maximum number of images as write kernel arguments

		size_t image2d_max_width;				///< maximum 2d image width
		size_t image2d_max_height;				///< maximum 2d image height

		size_t image3d_max_width;				///< maximum 3d image width
		size_t image3d_max_height;				///< maximum 3d image height
		size_t image3d_max_depth;				///< maximum 3d image depth

		cl_uint max_samplers;					///< maximum number of samplers
		size_t max_parameter_size;				///< maximum number of kernel parameters
		cl_uint mem_base_addr_align;			///< alignment of device base memory
		cl_uint min_data_type_align_size;		///< minimum alignment needed for device memory

		cl_device_fp_config single_fp_config;	///< single precision floating point capabilities
		cl_device_mem_cache_type global_mem_cache_type;	///< cache type of global memory
		cl_uint global_mem_cacheline_size;		///< size of a line of global memory cache

		cl_uint max_constant_args;				///< maximum number of constant arguments
		cl_device_local_mem_type local_mem_type;	///< type of local memory
		cl_bool error_correction;				///< error correction available
		size_t profiling_timer_resolution;		///< resolution of profiling timer
		cl_bool endian_little;					///< little endian device
		cl_bool compiler_available;				///< true, if compiler for device is available
		
		char *vendor;			///< vendor of device
		char *profile;			///< profile of device
		char *extensions;		///< extensions available for device
		char *version;			///< version of device
#endif	

		/**
		 * initialize device information with NULL data
		 */
		inline void initCDeviceInfo()
		{
			max_work_item_sizes = NULL;

			name = NULL;
//			vendor = NULL;
			driver_version = '\0';	//NULL;
//			profile = NULL;
//			version = NULL;
//			extensions = NULL;
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
//			delete[] vendor;
//			delete[] driver_version;
//			delete[] profile;
//			delete[] version;
//			delete[] extensions;

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
				case CL_DEVICE_TYPE_CPU:	return "CPU";
				case CL_DEVICE_TYPE_GPU:	return "GPU";
				case CL_DEVICE_TYPE_ACCELERATOR:	return "ACCELERATOR";
				case CL_DEVICE_TYPE_DEFAULT:	return "DEFAULT";
				case CL_DEVICE_TYPE_ALL:	return "ALL";
				default:			return "unknown";
			}
#endif
		}

		/**
		 * load device information given by device_id
		 */
		inline void loadDeviceInfo(const CDevice &cDevice)
		{
			set(cDevice.device_id);	// set device id
	
			getCudaDrvErrorString( cuDeviceGetCount(&max_compute_units) );
			max_work_item_dimensions = 3;

			delete max_work_item_sizes;
			max_work_item_sizes = new size_t[max_work_item_dimensions];
			
			getCudaDrvErrorString( cuDeviceGetAttribute(&max_work_group_size, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, device_id) );
			getCudaDrvErrorString( cuDeviceGetAttribute(&max_clock_frequency, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, device_id) );
			max_clock_frequency = max_clock_frequency * 1e-3f;	// convert to MHz from kHz
			
			getCudaDrvErrorString( cuDeviceGetAttribute(&address_bits, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, device_id) );
			getCudaDrvErrorString( cuDeviceGetAttribute(&global_mem_cache_size, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, device_id) );
			getCudaDrvErrorString( cuDeviceTotalMem(&global_mem_size, device_id) );
			getCudaDrvErrorString( cuDeviceGetAttribute(&max_constant_buffer_size, CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, device_id) );
			getCudaDrvErrorString( cuDeviceGetAttribute(&local_mem_size, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, device_id) );
			getCudaDrvErrorString( cuDeviceGetAttribute(&available, CU_DEVICE_ATTRIBUTE_COMPUTE_MODE, device_id) );
			getCudaDrvErrorString( cuDeviceComputeCapability(&execution_capabilities, NULL, device_id) );
			getCudaDrvErrorString( cuDeviceGetAttribute(&queue_properties, CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS, device_id) );

			getCudaDrvErrorString( cuDeviceGetName(name, 256, device_id) ); 
			getCudaDrvErrorString( cuDriverGetVersion(&driver_version) );
#if 0			
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_TYPE,			sizeof(cl_device_type),	&device_type, NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_VENDOR_ID,			sizeof(cl_uint),	&vendor_id,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_COMPUTE_UNITS,		sizeof(cl_uint),	&max_compute_units,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,	sizeof(cl_uint),	&max_work_item_dimensions,	NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group_size),	&max_work_group_size,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group_size),	&max_work_group_size,	NULL));

			// the role of this function is not clear? TODO: check
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*max_work_item_dimensions,	max_work_item_sizes,	NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,	sizeof(size_t),	&preferred_vector_width_char,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,	sizeof(size_t),	&preferred_vector_width_short,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,	sizeof(size_t),	&preferred_vector_width_int,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,	sizeof(size_t),	&preferred_vector_width_long,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,	sizeof(size_t),	&preferred_vector_width_float,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,	sizeof(size_t),	&preferred_vector_width_double,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY,		sizeof(size_t),	&max_clock_frequency,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_ADDRESS_BITS,		sizeof(size_t),		&address_bits,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE,	sizeof(cl_ulong),	&max_mem_alloc_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE_SUPPORT,		sizeof(cl_bool),	&image_support,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_READ_IMAGE_ARGS,	sizeof(cl_uint),	&max_read_image_args,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WRITE_IMAGE_ARGS,	sizeof(cl_uint),	&max_write_image_args,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE2D_MAX_WIDTH,		sizeof(size_t),		&image2d_max_width,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE2D_MAX_HEIGHT,	sizeof(size_t),		&image2d_max_height,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE3D_MAX_WIDTH,		sizeof(size_t),		&image3d_max_width,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE3D_MAX_HEIGHT,	sizeof(size_t),		&image3d_max_height,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE3D_MAX_DEPTH,		sizeof(size_t),		&image3d_max_depth,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_SAMPLERS,		sizeof(cl_uint),	&max_samplers,			NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_PARAMETER_SIZE,	sizeof(size_t),		&max_parameter_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MEM_BASE_ADDR_ALIGN,	sizeof(cl_uint),	&mem_base_addr_align,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE,	sizeof(cl_uint),	&min_data_type_align_size,	NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_SINGLE_FP_CONFIG,		sizeof(cl_device_fp_config),	&single_fp_config,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_CACHE_TYPE,	sizeof(cl_device_mem_cache_type),	&global_mem_cache_type,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,	sizeof(cl_uint),		&global_mem_cacheline_size,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,	sizeof(cl_ulong),		&global_mem_cache_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_SIZE,		sizeof(cl_ulong),		&global_mem_size,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,	sizeof(cl_ulong),	&max_constant_buffer_size,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_CONSTANT_ARGS,		sizeof(cl_uint),	&max_constant_args,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_LOCAL_MEM_TYPE,		sizeof(cl_device_local_mem_type),	&local_mem_type,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_LOCAL_MEM_SIZE,		sizeof(cl_ulong),	&local_mem_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_ERROR_CORRECTION_SUPPORT,	sizeof(cl_bool),	&error_correction,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PROFILING_TIMER_RESOLUTION,	sizeof(size_t),		&profiling_timer_resolution,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_ENDIAN_LITTLE,		sizeof(cl_bool),	&endian_little,			NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_AVAILABLE,			sizeof(cl_bool),	&available,			NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_COMPILER_AVAILABLE,		sizeof(cl_bool),	&compiler_available,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_EXECUTION_CAPABILITIES,	sizeof(cl_device_exec_capabilities),	&execution_capabilities,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_QUEUE_PROPERTIES,		sizeof(cl_command_queue_properties),	&queue_properties,		NULL));

			loadDeviceInfoString(device_id, CL_DEVICE_NAME, &name);
			loadDeviceInfoString(device_id, CL_DEVICE_VENDOR, &vendor);
			loadDeviceInfoString(device_id, CL_DRIVER_VERSION, &driver_version);
			loadDeviceInfoString(device_id, CL_DEVICE_PROFILE, &profile);
			loadDeviceInfoString(device_id, CL_DEVICE_VERSION, &version);
			loadDeviceInfoString(device_id, CL_DEVICE_EXTENSIONS, &extensions);
#endif
			
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


			std::cout << sLinePrefix << "LOCAL_MEM_SIZE: " << local_mem_size << std::endl;
			std::cout << sLinePrefix << "AVAILABLE: " << available << std::endl;
			std::cout << sLinePrefix << "COMPUTE_CAPABILITY: " << execution_capabilities << std::endl;
			std::cout << sLinePrefix << "CONCURRENT_KERNELS: " << queue_properties << std::endl;
			std::cout << sLinePrefix << "NAME: " << name << std::endl;
			std::cout << sLinePrefix << "DRIVER_VERSION: " << driver_version << std::endl;

#if 0
			std::cout << sLinePrefix << "TYPE: ";

			if (device_type & CL_DEVICE_TYPE_CPU)	std::cout << "CPU ";
			if (device_type & CL_DEVICE_TYPE_GPU)	std::cout << "GPU ";
			if (device_type & CL_DEVICE_TYPE_ACCELERATOR)	std::cout << "ACCELERATOR ";
			if (device_type & CL_DEVICE_TYPE_DEFAULT)	std::cout << "DEFAULT ";
			std::cout << std::endl;

			std::cout << sLinePrefix << "VENDOR_ID: " << vendor_id << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_CHAR: " << preferred_vector_width_char << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_SHORT: " << preferred_vector_width_short << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_INT: " << preferred_vector_width_int << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_LONG: " << preferred_vector_width_long << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_FLOAT: " << preferred_vector_width_float << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_DOUBLE: " << preferred_vector_width_double << std::endl;	
			std::cout << sLinePrefix << std::endl;	
			std::cout << sLinePrefix << "MAX_MEM_ALLOC_SIZE: " << max_mem_alloc_size << std::endl;
			std::cout << sLinePrefix << "IMAGE_SUPPORT: " << image_support << std::endl;
			std::cout << sLinePrefix << "MAX_READ_IMAGE_ARGS: " << max_read_image_args << std::endl;
			std::cout << sLinePrefix << "MAX_WRITE_IMAGE_ARGS: " << max_write_image_args << std::endl;
			std::cout << sLinePrefix << "IMAGE2D_MAX_WIDTH: " << image2d_max_width << std::endl;
			std::cout << sLinePrefix << "IMAGE2D_MAX_HEIGHT: " << image2d_max_height << std::endl;
			std::cout << sLinePrefix << "IMAGE3D_MAX_WIDTH: " << image3d_max_width << std::endl;
			std::cout << sLinePrefix << "IMAGE3D_MAX_HEIGHT: " << image3d_max_height << std::endl;
			std::cout << sLinePrefix << "IMAGE3D_MAX_DEPTH: " << image3d_max_depth << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "MAX_SAMPLERS: " << max_samplers << std::endl;
			std::cout << sLinePrefix << "MAX_PARAMETER_SIZE: " << max_parameter_size << std::endl;
			std::cout << sLinePrefix << "MEM_BASE_ADDR_ALIGN: " << mem_base_addr_align << std::endl;
			std::cout << sLinePrefix << "MIN_DATA_TYPE_ALIGN_SIZE: " << min_data_type_align_size << std::endl;
			std::cout << sLinePrefix << "SINGLE_FP_CONFIG: ";
			if (single_fp_config & CL_FP_DENORM)	std::cout << "FP_DENORM ";
			if (single_fp_config & CL_FP_INF_NAN)	std::cout << "FP_INF_NAN ";
			if (single_fp_config & CL_FP_ROUND_TO_NEAREST)	std::cout << "FP_ROUND_TO_NEAREST ";
			if (single_fp_config & CL_FP_ROUND_TO_ZERO)	std::cout << "FP_ROUND_TO_ZERO ";
			if (single_fp_config & CL_FP_ROUND_TO_INF)	std::cout << "FP_ROUND_TO_INF ";
			if (single_fp_config & CL_FP_FMA)	std::cout << "FP_FMA ";
			std::cout << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "GLOBAL_MEM_CACHE_TYPE: ";
			switch(global_mem_cache_type)
			{
				case CL_NONE:			std::cout << "NONE";	break;
				case CL_READ_ONLY_CACHE:	std::cout << "CL_READ_ONLY_CACHE";	break;
				case CL_READ_WRITE_CACHE:	std::cout << "CL_READ_WRITE_CACHE";	break;
			}
			std::cout << std::endl;
			std::cout << sLinePrefix << "GLOBAL_MEM_CACHELINE_SIZE: " << global_mem_cacheline_size << std::endl;
			std::cout << sLinePrefix << "MAX_CONSTANT_ARGS: " << max_constant_args << std::endl;
			std::cout << sLinePrefix << "LOCAL_MEM_TYPE: ";
			switch(local_mem_type)
			{
				case CL_LOCAL:	std::cout << "LOCAL";	break;
				case CL_GLOBAL:	std::cout << "GLOBAL";	break;
				default:	std::cout << "UNKNOWN";	break;
			}
			std::cout << std::endl;
			std::cout << sLinePrefix << "ERROR_CORRECTION_SUPPORT: " << error_correction << std::endl;
			std::cout << sLinePrefix << "PROFILING_TIMER_RESOLUTION: " << profiling_timer_resolution << std::endl;
			std::cout << sLinePrefix << "ENDIAN_LITTLE: " << endian_little << std::endl;
			std::cout << sLinePrefix << "COMPILER_AVAILABLE: " << compiler_available << std::endl;
			std::cout << sLinePrefix << "EXECUTION_CAPABILITIES: ";
			if (execution_capabilities & CL_EXEC_KERNEL)		std::cout << "EXEC_KERNEL ";
			if (execution_capabilities & CL_EXEC_NATIVE_KERNEL)	std::cout << "EXEC_NATIVE_KERNEL ";
			std::cout << std::endl;
			std::cout << sLinePrefix << "QUEUE_PROPERTIES: ";
			if (queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)		std::cout << "QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE ";
			if (queue_properties & CL_QUEUE_PROFILING_ENABLE)	std::cout << "QUEUE_PROFILING_ENABLE";
			std::cout << std::endl;
			std::cout << sLinePrefix << "VENDOR: " << vendor << std::endl;
			std::cout << sLinePrefix << "PROFILE: " << profile << std::endl;
			std::cout << sLinePrefix << "VERSION: " << version << std::endl;
			std::cout << sLinePrefix << "EXTENSIONS: " << extensions << std::endl;
#endif
		}
#if 0
private:
		inline void loadDeviceInfoString(	cl_device_id device_id,
						cl_device_info device_info,
						char **param_value)
		{
			size_t retval_size;
			CL_CHECK_ERROR(clGetDeviceInfo(device_id, device_info, 0, NULL, &retval_size));
			delete *param_value;
			*param_value = new char[retval_size];
			CL_CHECK_ERROR(clGetDeviceInfo(device_id, device_info, retval_size, *param_value, NULL));
		}
#endif
	};
};


#endif
