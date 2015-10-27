################################################################################
# paths, directories and folders
################################################################################

CUDAINSTALLPATH		:=	/usr/local/cuda

CCLIBDIR			:=	-L/opt/netcdf/4.3.3.1/lib
CXXLIBDIR			:=	-L/opt/netcdf/4.3.3.1/lib
CUDALIBDIR			:=	

CCINCLUDES			:=	-I/usr/lib/openmpi/include \
						-I/opt/netcdf/4.3.3.1/include
CXXINCLUDES			:=	-I/usr/lib/openmpi/include \
						-I/opt/netcdf/4.3.3.1/include
CUDAINCLUDES		:=	

CCLIB				:=	
CXXLIB				:=	
CUDALIB				:=	

COMPUTE_CAPABILITY	:=	50

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
#						-D PAR_NETCDF \
#						-std=c11
CXXFLAGS			:=	-O3 \
#						-D PAR_NETCDF \
#						-std=c++11

# arch: specifies the compatibility from source code to PTX stage. Can be a
#       virtual (compute_*) or real (sm_*) compatibility.
# code: specifies the compatibility from PTX stage to binary code. Can only be
#       real (sm_*). Code has to be >= arch.
# -rdc: -rdc is short for --relocatable-device-code which generates relocatable
#       device code. This is necessary to generate multiple CUDA object files
#       which can then be linked together.
NVCCFLAGS			:=	-O3 \
						-gencode arch=compute_$(COMPUTE_CAPABILITY),code=sm_$(COMPUTE_CAPABILITY) \
						-D PAR_NETCDF \
#						-Xcompiler "-std=c++11"

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS			:=	

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS		:=	-arch=sm_$(COMPUTE_CAPABILITY)

include common.make

