################################################################################
# paths, directories and folders
################################################################################

BINDIR			:=	bin
OBJDIR			:=	obj

CCLIBDIR		+= 
CXXLIBDIR		+= 
CUDABINDIR		:=	$(CUDAINSTALLPATH)/bin
CUDALIBDIR		+=	-L$(CUDAINSTALLPATH)/lib64

CCINCLUDES		+=	-I$(CUDAINSTALLPATH)/include
CXXINCLUDES		+=	-I$(CUDAINSTALLPATH)/include
CUDAINCLUDES	+=	-I$(CUDAINSTALLPATH)/include

CCLIB			+=	-lnetcdf \
					-lrt
CXXLIB			+=	-lnetcdf \
					-lrt
# libcuda:      required for access to driver api
# libcudart:    required for execution of a cuda program
# libcudadevrt: required for dynamic parallelism which is again required for
#               valid linking of multiple cuda files
CUDALIB			+=	-lcudart \
					-lcudadevrt

EXECUTABLE		:=	lbm

################################################################################
# source files
################################################################################

# c/c++ source files (compiled with $(CC))
CCFILES			+=	

# c/c++ source files (compiled with $(CXX))
CXXFILES		+=	external/tinyxml2/tinyxml2.cpp \
					src/cpukernels/CLbmInitCPU.cpp \
					src/cpukernels/CLbmAlphaCPU.cpp \
					src/cpukernels/CLbmBetaCPU.cpp \
					src/libmath/CMath.cpp \
					src/libmath/CVector2.cpp \
					src/libmath/CVector3.cpp \
					src/libmath/CVector4.cpp \
					src/libvis/CLbmVisualizationVTK.cpp \
					src/CConfiguration.cpp \
					src/CComm.cpp \
					src/CController.cpp \
					src/CDomain.cpp \
					src/CManager.cpp \
					src/CLbmSolver.cpp \
					src/CLbmSolverCPU.cpp \
					src/main.cpp

# cuda source files (compiled with $(NVCC))
CUDAFILES		+=	src/gpukernels/lbm_alpha.cu \
					src/gpukernels/lbm_beta.cu \
					src/gpukernels/lbm_init.cu \
					src/CLbmSolverGPU.cu

################################################################################
# compiler arguments and flags
################################################################################

CCFLAGS			+=	-Wall
CXXFLAGS		+=	-Wall

# arch: specifies the compatibility from source code to PTX stage. Can be a
#       virtual (compute_*) or real (sm_*) compatibility.
# code: specifies the compatibility from PTX stage to binary code. Can only be
#       real (sm_*). Code has to be >= arch.
# -rdc: -rdc is short for --relocatable-device-code which generates relocatable
#       device code. This is necessary to generate multiple CUDA object files
#       which can then be linked together.
NVCCFLAGS		+=	-lineinfo \
					-rdc=true \
					-use_fast_math \
					--compiler-options -Wall \
#					--ptxas-options -v

################################################################################
# linker arguments and flags
################################################################################

LINKERFLAGS		+=	

# -dlink: Necessary linker option to link multiple CUDA object files together.
NVCCLINKERFLAGS	+=	-dlink

################################################################################
# set up virtual path to enable subfolders for source files
################################################################################

VPATH 			:=	external/tinyxml2/ \
					src/cpukernels/ \
					src/gpukernels/ \
					src/libmath/ \
					src/libtools/ \
					src/libvis/ \
					src/

################################################################################
# set up object files
#
# semantics patsubst(a, b, c): replace b by a in c.
################################################################################

CCOBJS			:=	$(patsubst %.c,   $(OBJDIR)/%.c.o,   $(notdir $(CCFILES)))
CXXOBJS			:=	$(patsubst %.cpp, $(OBJDIR)/%.cpp.o, $(notdir $(CXXFILES)))
CUDAOBJS		:=	$(patsubst %.cu,  $(OBJDIR)/%.cu.o,  $(notdir $(CUDAFILES)))

CCXXOBJS		:=	$(CCOBJS)
CCXXOBJS		+=	$(CXXOBJS)

OBJS			:=  $(CCOBJS)
OBJS			+=  $(CXXOBJS)
OBJS			+=  $(CUDAOBJS)

################################################################################
# set up link process
################################################################################

LINKLINE		:=	$(LINKER) $(LINKERFLAGS) $(OBJS) $(OBJDIR)/cuda.cu.o -o $(BINDIR)/$(EXECUTABLE) $(CCLIBDIR) $(CXXLIBDIR) $(CUDALIBDIR) $(CCLIB) $(CXXLIB) $(CUDALIB)
NVCCLINKLINE	:=	$(NVCCLINKER) $(NVCCLINKERFLAGS) $(CUDAOBJS) -o $(OBJDIR)/cuda.cu.o

################################################################################
# targets
################################################################################

# target to compile c files
$(OBJDIR)/%.c.o: %.c
	$(CC) $(CCINCLUDES) $(CCFLAGS) -c $< -o $@

# target to compile c++ files
$(OBJDIR)/%.cpp.o: %.cpp
	$(CXX) $(CXXINCLUDES) $(CXXFLAGS) -c $< -o $@

# target to compile cuda files
$(OBJDIR)/%.cu.o: %.cu
	$(NVCC) $(CUDAINCLUDES) $(NVCCFLAGS) -c $< -o $@
	
# misc targets
makedirectories:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

clean:
	rm -f $(OBJDIR)/*
	rm -f $(BINDIR)/$(EXECUTABLE)
	rmdir $(OBJDIR)
	rmdir $(BINDIR)
	
# link targets (results are executables)
link: makedirectories $(CCXXOBJS) cudaobject
	@echo '-- Invoking C/C++ linker: Link C/C++ objects, CUDA objects and single CUDA object --'
	$(LINKLINE)
	@echo '-- End invoking C/C++ linker --'
	
# compile targets (results are object files)
cudaobject: $(CUDAOBJS)
	@echo '-- Invoking CUDA linker: Linking all CUDA objects to one single object --'
	$(NVCCLINKLINE)
	@echo '-- End invoking CUDA linker --'

# frontend targets (sould be called as make option)
all: link
	@echo '-- Everything went fine --'

