################################################################################
# paths, directories and folders
################################################################################

CUDAINSTALLPATH		:=	$(CUDA_HOME)

CCLIBDIR			:= -L$(NETCDF_HOME)/lib
CXXLIBDIR			:= -L$(NETCDF_HOME)/lib
CUDALIBDIR			:=	

CCINCLUDES			:=	-I$(NETCDF_HOME)/include
CXXINCLUDES			:=	-I$(NETCDF_HOME)/include
CUDAINCLUDES		:=	

CCLIB				:=	
CXXLIB				:=	
CUDALIB				:=	

COMPUTE_CAPABILITY	:=	35

################################################################################
# compilers and linkers
################################################################################

CC					:=	mpicc
CXX					:=	mpicxx
NVCC				:=	$(CUDAINSTALLPATH)/bin/nvcc
LINKER				:=	mpicxx
NVCCLINKER			:=	$(CUDAINSTALLPATH)/bin/nvcc

################################################################################
# compiler arguments and flags
################################################################################

CCFLAGS				:=	-O3 \
#						-D DEBUG \
#						-std=c0x
CXXFLAGS			:=	-O3 \
#						-D DEBUG \
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
#						-D DEBUG \
#						-Xcompiler "-std=c++0x"

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS			:=	

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS		:=	-arch=sm_$(COMPUTE_CAPABILITY)

include common.make
