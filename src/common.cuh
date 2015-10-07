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

#ifndef COMMON_CUH
#define COMMON_CUH

/*
 * Workaround taken from 0_Simple/simpleTemplates/sharedmem.cuh:
 * Because dynamically sized shared memory arrays are declared "extern", we
 * can't templatize them directly.  To get around this, we declare a simple
 * wrapper struct that will declare the extern array with a different name
 * depending on the type.  This avoids compiler errors about duplicate
 * definitions.
 *
 * To use dynamically allocated shared memory in a templatized __global__ or
 * __device__ function, just replace code like this:
 *
 * template<class T>
 * __global__ void
 * foo( T* g_idata, T* g_odata)
 * {
 *     // Shared mem size is determined by the host app at run time
 *     extern __shared__  T sdata[];
 *     ...
 *     doStuff(sdata);
 *     ...
 * }
 *
 * With this
 * template<class T>
 * __global__ void
 * foo( T* g_idata, T* g_odata)
 * {
 *     // Shared mem size is determined by the host app at run time
 *     SharedMemory<T> smem;
 *     T* sdata = smem.getPointer();
 *     ...
 *     doStuff(sdata);
 *     ...
 * }
 */

/*
 * This is the un-specialized struct.  Note that we prevent instantiation of
 * this struct by putting an undefined symbol in the function body so it won't
 * compile.
 */
template <typename T>
struct SharedMemory
{
    // Ensure that we won't compile any un-specialized types
    __device__ T *getPointer()
    {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

/*
 * Following are the specializations for float and double.
 */
template <>
struct SharedMemory<float>
{
    __device__ float *getPointer()
    {
        extern __shared__ float s_float[];
        return s_float;
    }
};

template <>
struct SharedMemory<double>
{
    __device__ double *getPointer()
    {
        extern __shared__ double s_double[];
        return s_double;
    }
};

#endif
