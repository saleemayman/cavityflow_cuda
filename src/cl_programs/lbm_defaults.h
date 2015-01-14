#include "vector_types.h"

#ifndef DOMAIN_CELLS_X
	#define DOMAIN_CELLS_X	(64)
	#define DOMAIN_CELLS_Y	(64)
	#define DOMAIN_CELLS_Z	(64)

	#define LOCAL_WORK_GROUP_SIZE	(128)

	typedef float	T;
	typedef float2	T2;
	typedef float4	T4;	
//	typedef float8	T8;	

	#ifndef float8
		struct float8 {
			float a, b, c, d, e, f, g, h;
		};

		typedef struct float8 T8;
		//typedef float8 T8;
	#endif

	#define FLAG_OBSTACLE	(1 << 0)
	#define FLAG_FLUID	(1 << 1)
	#define FLAG_VELOCITY_INJECTION	(1 << 2)
	#define FLAG_GHOST_LAYER (1 << 3)

	#define SIZE_DD_HOST_BYTES	(19*sizeof(T))

	#define STORE_VELOCITY	1
	#define STORE_DENSITY	1

#endif