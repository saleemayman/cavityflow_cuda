################################################################################
# paths, directories and folders
################################################################################

CUDAINSTALLPATH		:=	/usr/local/cuda

CCLIBDIR			:= 
CXXLIBDIR			:= 
CUDALIBDIR			:=	

CCINCLUDES			:=	-I/usr/include \
						-I/usr/lib/openmpi/include
CXXINCLUDES			:=	-I/usr/include \
						-I/usr/lib/openmpi/include
CUDAINCLUDES		:=	-I/usr/include \
						-I/usr/lib/openmpi/include

CCLIB				:=	
CXXLIB				:=	
CUDALIB				:=	

COMPUTE_CAPABILITY	:=	20

################################################################################
# compilers and linkers
################################################################################

CC					:=	gcc
CXX					:=	g++
NVCC				:=	$(CUDAINSTALLPATH)/bin/nvcc
LINKER				:=	g++
NVCCLINKER			:=	$(CUDAINSTALLPATH)/bin/nvcc

################################################################################
# compiler arguments and flags
################################################################################

CCFLAGS				:=	-O3 \
#						-std=c11
CXXFLAGS			:=	-O3 \
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
#						-Xcompiler "-std=c++11"

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS			:=	

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS		:=	-arch=sm_$(COMPUTE_CAPABILITY)

include common.make

