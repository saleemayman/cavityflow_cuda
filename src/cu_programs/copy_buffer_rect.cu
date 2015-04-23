#include "lbm_defaults.h"

/**
 * COPY KERNEL
 *
 * Copy a block of  data from src to dst.
 * Based on the assumption that src and dst values are saved in x, y, z order.
 * Since each row of block in x direction is copied with one work item,
 * the number of work items to use with this kernel should be (block_size_y*block_size_z)
 *
 */
/*

/* 	Comments CUDA porting
 *
 *	Check if the variable address space qualifiers are correct.
 *	For CUDA address space is defined outside the function scope.
*/
extern "C" __global__ void copy_buffer_rect(
				T*	src,
		const 	int src_offset,
		const 	int src_origin_x,
		const 	int src_origin_y,
		const 	int src_origin_z,
		const 	int src_size_x,
		const 	int src_size_y,
		const 	int src_size_z,
				T*	dst,
		const 	int dst_offset,
		const 	int dst_origin_x,
		const 	int dst_origin_y,
		const 	int dst_origin_z,
		const 	int dst_size_x,
		const 	int dst_size_y,
		const 	int dst_size_z,
		const 	int block_size_x
		)
{
	/*
	int gid_j = get_global_id(0);
	int gid_k = get_global_id(1);
	*/
	// get the thread ID for each direction
	// int gid_j = threadIdx.x + blockDim.x * blockIdx.x;
	// int gid_k = threadIdx.y + blockDim.y * blockIdx.y;
	const size_t idx_x = threadIdx.x + blockDim.x * blockIdx.x;		
	const size_t idx_y = threadIdx.y + blockDim.y * blockIdx.y;
	const size_t idx_xy = idx_y * (blockDim.x * gridDim.x) + idx_x;
	int gid_j = idx_x;
	int gid_k = idx_y;

	// if (gid_j >= DOMAIN_CELLS_Y || gid_k >= DOMAIN_CELLS_Z)
	// 	return;

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

	for	(x_idx = 0; x_idx < block_size_x; x_idx++)
	{
		dst[ dst_origin_idx + x_idx] = src[src_origin_idx + x_idx];
	}

}
/*

const size_t idx_x = threadIdx.x + blockDim.x * blockIdx.x;
const size_t idx_y = threadIdx.y + blockDim.y * blockIdx.y;
const size_t idx_z = threadIdx.z + blockDim.z * blockIdx.z;

const size_t idx_xy = idx_y * (blockDim.x * gridDim.x) + idx_x;

src_offset: 0
src_origin_x[0]: 0
src_origin_y[1]: 0
src_origin_z[2]: 0
src_size_x[0]: 14
src_size_y[1]: 1
src_size_z[2]: 2
dst_offset: 0
dst_origin_x[0]: 1
dst_origin_y[1]: 14
dst_origin_z[2]: 1
dst_size_x[0]: 16
dst_size_y[1]: 16
dst_size_z[2]: 4
block_size[0]: 14
lGlobalSize[0]: 1
lGlobalSize[1]: 2
threadsPerBlock: [16, 16, 1] 
numBlocks: [3, 3, 1] 
*/