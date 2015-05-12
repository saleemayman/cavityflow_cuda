/*
 * common_cpu_gpu.h
 *
 *  Created on: May 7, 2015
 *      Author: Ayman Saleem
 */


#ifndef COMMON_ALL_H_
#define COMMON_ALL_H_

#define FLAG_OBSTACLE	(1 << 0)
#define FLAG_FLUID	(1 << 1)
#define FLAG_VELOCITY_INJECTION	(1 << 2)
#define FLAG_GHOST_LAYER (1 << 3)
#define FLAG_GHOST_LAYER_BETA (FLAG_GHOST_LAYER | (1 << 4))

#define STORE_VELOCITY	1
#define STORE_DENSITY	1

// simulation data type
#define TYPE_FLOAT 1
typedef float	T;

// #define TYPE_DOUBLE 1
// typedef double	T;


#endif /* COMMON_ALL_H_ */
