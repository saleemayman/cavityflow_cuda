all:
	nvcc -x cu -keep -keep-dir src/cu_programs                             \
	-I.  src/cu_programs/lbm_init.cu -o src/cu_programs/lbm_init.o

	nvcc -x cu -keep -keep-dir src/cu_programs                             \
	-I.  src/cu_programs/lbm_alpha.cu -o src/cu_programs/lbm_alpha.o

	nvcc -x cu -keep -keep-dir src/cu_programs                            \
	-I.  src/cu_programs/lbm_beta.cu -o src/cu_programs/lbm_beta.o

	nvcc -x cu -keep -keep-dir src/cu_programs                             
	-I.  src/cu_programs/copy_buffer_rect.cu -o src/cu_programs/copy_buffer_rect.o

ifeq ($(strip $(PROF)),)
#	scons --compiler=$(COMP) 
else
	#scons --compiler=$(strip $(COMP)) --profiler=$(PROF)
#	scons --compiler=mpiCC --profiler=$(PROF)
endif


doc:
	doxygen ./docs/Doxyfile

cleandoc:
	rm -rf ./docs/html ./docs/latex

cleanvtk:
	rm -f ./output/vtk/*

cleanbenchmark:
	rm -f ./output/benchmark/*

cleanprofile:
	rm -f ./output/profile/*

cleanlog: 
	rm -f ./output/log/*

cleanbuild:
	rm -rf build/*

cleanresultsoutput:
	rm -rf output

clean: cleanbuild cleanlog cleanprofile cleanbenchmark cleanvtk
	rm src/cu_programs/*.cuda* src/cu_programs/*.cpp* src/cu_programs/*.fatbin* src/cu_programs/*.hash src/cu_programs/*.module_id src/cu_programs/*.sm_*
