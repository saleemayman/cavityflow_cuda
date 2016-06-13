################################################################################
# paths, directories and folders
################################################################################

CUDAINSTALLPATH		:=	$(CUDA_HOME)

CCLIBDIR			:=	-L$(NETCDF_LIBDIR)
CXXLIBDIR			:=	-L$(NETCDF_LIBDIR)
CUDALIBDIR			:=	

CCINCLUDES			:=	$(NETCDF_INC)
ifeq ($(USE_MPI), 1)
CCINCLUDES			+=	$(MPI_INC)
endif
CXXINCLUDES			:=	$(NETCDF_INC)
ifeq ($(USE_MPI), 1)
CXXINCLUDES			+=	$(MPI_INC)
endif
CUDAINCLUDES		:=	

CCLIB				:=	-lgomp \
						-lnetcdf
CXXLIB				:=	-lgomp \
						-lnetcdf
CUDALIB				:=	

COMPUTE_CAPABILITY	:=	20

################################################################################
# compilers and linkers
################################################################################

CC					:=	
CXX					:=	
LINKER				:=	

ifeq ($(INSTRUMENT), scalasca)
CC					+=	scalasca -instrument
CXX					+=	scalasca -instrument
LINKER				+=	scalasca -instrument
endif
ifeq ($(INSTRUMENT), scorep)
CC					+=	scorep
CXX					+=	scorep
LINKER				+=	scorep
endif

ifeq ($(USE_MPI), 1)
CC					+=	mpicc
CXX					+=	mpicxx
LINKER				+=	mpicxx
else
CC					:=	gcc
CXX					:=	g++
LINKER				:=	g++
endif

NVCC				:=	$(CUDAINSTALLPATH)/bin/nvcc
NVCCLINKER			:=	$(CUDAINSTALLPATH)/bin/nvcc

################################################################################
# compiler arguments and flags
################################################################################

CCFLAGS				:=	-O3 \
						-fopenmp \
#						-xHost \
#						-std=c11
CXXFLAGS			:=	-O3 \
						-fopenmp \
#						-xHost \
#						-std=c11

# arch: specifies the compatibility from source code to PTX stage. Can be a
#       virtual (compute_*) or real (sm_*) compatibility.
# code: specifies the compatibility from PTX stage to binary code. Can only be
#       real (sm_*). Code has to be >= arch.
# -rdc: -rdc is short for --relocatable-device-code which generates relocatable
#       device code. This is necessary to generate multiple CUDA object files
#       which can then be linked together.
NVCCFLAGS			:=	-O3 \
						-gencode arch=compute_$(COMPUTE_CAPABILITY),code=sm_$(COMPUTE_CAPABILITY) \
#						-maxrregcount=56 \
#						--ptxas-options -v \
#						-Xcompiler "-std=c++0x"

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS			:=	

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS		:=	-arch=sm_$(COMPUTE_CAPABILITY)

include common.make

