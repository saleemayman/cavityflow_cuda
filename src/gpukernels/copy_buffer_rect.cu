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

#include "copy_buffer_rect.cuh"

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
        const int block_size_z)
{
    // get the thread ID for each direction
    const size_t gid_j = threadIdx.x + blockDim.x * blockIdx.x;     
    const size_t gid_k = threadIdx.y + blockDim.y * blockIdx.y;

    if (gid_j >= block_size_y || gid_k >= block_size_z)
        return;

    int src_slice_z_size = src_size_x * src_size_y;
    int dst_slice_z_size = dst_size_x * dst_size_y;

    // index of the origin of the block in the src array
    int src_origin_idx = src_offset + src_origin_x +
                        (src_origin_y + gid_j) * src_size_x +
                        (src_origin_z + gid_k) * src_slice_z_size;

    // index of the origin of the block in the dst array
    int dst_origin_idx = dst_offset + dst_origin_x +
                        (dst_origin_y + gid_j) * dst_size_x +
                        (dst_origin_z + gid_k) * dst_slice_z_size;
    int x_idx;

    for (x_idx = 0; x_idx < block_size_x; x_idx++)
    {
        dst[ dst_origin_idx + x_idx] = src[src_origin_idx + x_idx];
    }
}

template __global__ void copy_buffer_rect<float>(
        float* src,
        const int src_offset,
        const int src_origin_x,
        const int src_origin_y,
        const int src_origin_z,
        const int src_size_x,
        const int src_size_y,
        const int src_size_z,
        float* dst,
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
template __global__ void copy_buffer_rect<double>(
        double* src,
        const int src_offset,
        const int src_origin_x,
        const int src_origin_y,
        const int src_origin_z,
        const int src_size_x,
        const int src_size_y,
        const int src_size_z,
        double* dst,
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