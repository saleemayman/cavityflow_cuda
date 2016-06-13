################################################################################
# paths, directories and folders
################################################################################

CUDAINSTALLPATH		:=	$(CUDA_HOME)

CCLIBDIR			:=	-L$(NETCDF_HOME)/lib
CXXLIBDIR			:=	-L$(NETCDF_HOME)/lib
CUDALIBDIR			:=	

CCINCLUDES			:=	-I$(NETCDF_HOME)/include
CXXINCLUDES			:=	-I$(NETCDF_HOME)/include
CUDAINCLUDES		:=	

CCLIB				:=	-lnetcdf
CXXLIB				:=	-lnetcdf
CUDALIB				:=	

COMPUTE_CAPABILITY	:=	35

################################################################################
# compilers and linkers
################################################################################

ifeq ($(USE_MPI), 1)
CC					:=	mpiicc
CXX					:=	mpiicpc
LINKER				:=	mpiicpc
else
CC					:=	icc
CXX					:=	icpc
LINKER				:=	icpc
endif

NVCC				:=	$(CUDAINSTALLPATH)/bin/nvcc
NVCCLINKER			:=	$(CUDAINSTALLPATH)/bin/nvcc

################################################################################
# compiler arguments and flags
################################################################################

CCFLAGS				:=	-O3 \
						-qopenmp \
						-parallel \
						-xHost \
						-DMPICH_IGNORE_CXX_SEEK \
#						-std=c11
CXXFLAGS			:=	-O3 \
						-qopenmp \
						-parallel \
						-xHost \
						-DMPICH_IGNORE_CXX_SEEK \
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
#						-maxrregcount=136 \
#						--ptxas-options -v \
#						-Xcompiler "-std=c++0x" \

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS			:=	-parallel

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS		:=	-arch=sm_$(COMPUTE_CAPABILITY)

include common.make

