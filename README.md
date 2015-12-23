#+TITLE: README for lbm
#+AUTHOR: Christoph Riesinger
#+EMAIL: riesinge@in.tum.de

Features
========

Installation
============
lbm is build via make. For every target system, an individual makefile with system specific parameters have to be provided. In addition, three additional options can be passed to make influencing the resulting executable. There are two build targets: clean and all.

All system specific individual makefiles have to include common.make at the end of the file. They enable to set system specific parameters like compilers, compute capabilities, compiler options, etc.

There are three options which can be passed to make:
- USE_MPI: If set to 1, a MPI capable parallel binary is generated. Having a non-MPI version can be useful, e.g. to profile GPU kernels.
- PAR_NETCDF: If set to 1, parallel netCDF is enabled. This is especially useful if it comes to output of simulation data on many ranks at the same time. Not all systems support parallel but only serial netCDF.
- INSTRUMENT: To profile lbm, the code has to be instrumented during compile time. Two different instrument options are available which lead to the usage of the corresponding tools: scalasca and scorep.

A valid make call to compile lbm could look as follows:
make -f mac-cluster.make -j32 all USE_MPI=1 PAR_NETCDF=1 INSTRUMENT=scorep

Running
=======

