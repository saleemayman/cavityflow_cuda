#include "vector_types.h"

#if PROBLEM_DEFAULTS
	#define DOMAIN_CELLS_X	(1)
	#define DOMAIN_CELLS_Y	(1)
	#define DOMAIN_CELLS_Z	(1)

	#define LOCAL_WORK_GROUP_SIZE	(16)


	#define FLAG_OBSTACLE	(1 << 0)
	#define FLAG_FLUID	(1 << 1)
	#define FLAG_VELOCITY_INJECTION	(1 << 2)
	#define FLAG_GHOST_LAYER (1 << 3)

	#define SIZE_DD_HOST_BYTES	(19*sizeof(T))

	#define STORE_VELOCITY	1
	#define STORE_DENSITY	1
#endif	//PROBLEM_DEFAULTS

#ifndef SIMULATION_DATA_TYPE
#define SIMULATION_DATA_TYPE
	#if TYPE_FLOAT
		typedef float	T;
		typedef float2	T2;
		typedef float4	T4;
	//	typedef float8	T8;	

		#ifndef float8
			struct float8 {
				float a, b, c, d, e, f, g, h;
			};

			typedef struct float8 T8;
		#endif

		#ifndef float16
			struct float16 {
				float a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p;
			};

			typedef struct float16 T16;
		#endif
	#elif TYPE_DOUBLE
		typedef double	T;
		typedef double2	T2;
		typedef double4	T4;
	//	typedef double8	T8;	

		#ifndef double8
			struct double8 {
				double a, b, c, d, e, f, g, h;
			};

			typedef struct double8 T8;
		#endif

		#ifndef double16
			struct double16 {
				double a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p;
			};

			typedef struct double16 T16;
		#endif
	#endif
#endif	//SIMULATION_DATA_TYPE