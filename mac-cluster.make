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

CCLIB				:=	
CXXLIB				:=	
CUDALIB				:=	

COMPUTE_CAPABILITY	:=	20

################################################################################
# compilers and linkers
################################################################################

ifeq ($(USE_MPI), 1)
CC					:=	mpicc
else
CC					:=	gcc
endif

ifeq ($(USE_MPI), 1)
CXX					:=	mpicxx
else
CXX					:=	g++
endif
NVCC				:=	$(CUDAINSTALLPATH)/bin/nvcc

ifeq ($(USE_MPI), 1)
LINKER				:=	mpicxx
else
LINKER				:=	g++
endif
NVCCLINKER			:=	$(CUDAINSTALLPATH)/bin/nvcc

################################################################################
# compiler arguments and flags
################################################################################

CCFLAGS				:=	-O3 \
#						-std=c0x
CXXFLAGS			:=	-O3 \
#						-std=c++0x

# arch: specifies the compatibility from source code to PTX stage. Can be a
#       virtual (compute_*) or real (sm_*) compatibility.
# code: specifies the compatibility from PTX stage to binary code. Can only be
#       real (sm_*). Code has to be >= arch.
# -rdc: -rdc is short for --relocatable-device-code which generates relocatable
#       device code. This is necessary to generate multiple CUDA object files
#       which can then be linked together.
NVCCFLAGS			:=	-O3 \
						-gencode arch=compute_$(COMPUTE_CAPABILITY),code=sm_$(COMPUTE_CAPABILITY) \
#						--ptxas-options -v
#						-Xcompiler "-std=c++0x"

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS			:=	

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS		:=	-arch=sm_$(COMPUTE_CAPABILITY)

include common.make

