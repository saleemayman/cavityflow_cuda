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

#ifndef COPY_BUFFER_RECT_CUH
#define COPY_BUFFER_RECT_CUH

#include "../common.h"

/*
 * COPY KERNEL
 *
 * Copy a block of  data from src to dst.
 * Based on the assumption that src and dst values are saved in x, y, z order.
 * Since each row of block in x direction is copied with one work item,
 * the number of work items to use with this kernel should be (block_size_y*block_size_z)
 */

/*
 * Comments CUDA porting
 *
 *  Check if the variable address space qualifiers are correct.
 *  For CUDA address space is defined outside the function scope.
 */

template<typename T>
__global__ void copy_buffer_rect(
        T* src,
        const int src_offset,
        const int src_origin_x,
        const int src_origin_y,
        const int src_origin_z,
        const int src_size_x,
        const int src_size_y,
        const int src_size_z,
        T* dst,
        const int dst_offset,
        const int dst_origin_x,
        const int dst_origin_y,
        const int dst_origin_z,
        const int dst_size_x,
        const int dst_size_y,
        const int dst_size_z,
        const int block_size_x,
        const int block_size_y,
        const int block_size_z);

#endif
