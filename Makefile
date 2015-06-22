all:
#	nvcc -x cu -keep -keep-dir src/cu_programs                             \
#	-I.  src/cu_programs/lbm_init.cu -o src/cu_programs/lbm_init.o
#
#	nvcc -x cu -keep -keep-dir src/cu_programs                             \
#	-I.  src/cu_programs/lbm_alpha.cu -o src/cu_programs/lbm_alpha.o
#
#	nvcc -x cu -keep -keep-dir src/cu_programs                            \
#	-I.  src/cu_programs/lbm_beta.cu -o src/cu_programs/lbm_beta.o
#
#	nvcc -x cu -keep -keep-dir src/cu_programs                             
#	-I.  src/cu_programs/copy_buffer_rect.cu -o src/cu_programs/copy_buffer_rect.o
	nvcc -ptx  src/cu_programs/lbm_init.cu -o src/cu_programs/lbm_init.ptx
                                                                 
	nvcc -ptx  src/cu_programs/lbm_alpha.cu -o src/cu_programs/lbm_alpha.ptx
                                                                   
	nvcc -ptx  src/cu_programs/lbm_beta.cu -o src/cu_programs/lbm_beta.ptx

	nvcc -ptx  src/cu_programs/copy_buffer_rect.cu -o src/cu_programs/copy_buffer_rect.ptx

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
	rm src/cu_programs/*.ptx
