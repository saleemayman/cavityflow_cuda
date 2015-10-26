/*
 * Copyright
 * 2010 Martin Schreiber
 * 2013 Arash Bakhtiari
 * 2016 Christoph Riesinger, Ayman Saleem
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "CLbmVisualizationNetCDF.hpp"

#include <cstdio>
#include <typeinfo>

#include <mpi.h>

template <class T>
CLbmVisualizationNetCDF<T>::CLbmVisualizationNetCDF(
		int id,
		int visualizationRate,
		CLbmSolver<T>* solver,
		std::string filePath) :
		CLbmVisualization<T>(id, visualizationRate, solver),
		filePath(filePath)
{
}

template <class T>
void CLbmVisualizationNetCDF<T>::openFile(int iteration)
{
	if (true)
	{
		int error;
	    std::stringstream fileName;
#ifdef PAR_NETCDF
	    fileName << filePath << "/visualization_" << iteration << ".nc";
	    error = nc_create_par(fileName.str().c_str(), NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &fileId);
#else
	    fileName << filePath << "/visualization_" << id << "_" << iteration << ".nc";
	    error = nc_create(fileName.str().c_str(), NC_CLOBBER, &fileId);
#endif

	    if (error)
	    {
	        std::cerr << "----- CLbmVisualizationNetCDF<T>::openFile() -----" << std::endl;
	        std::cerr << "An error occurred while opening parallel netCDF file \"" << fileName.str() << "\"" << std::endl;std::cerr << nc_strerror(error) << std::endl;
	        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
	        std::cerr << "--------------------------------------------------" << std::endl;

	        exit (EXIT_FAILURE);
	    }
	} else {
        std::cerr << "----- CLbmVisualizationNetCDF<T>::openFile() -----" << std::endl;
#ifdef PAR_NETCDF
        std::cerr << "netCDF file \"visualization_" << iteration << ".nc is already open" << std::endl;
#else
        std::cerr << "netCDF file \"visualization_" << id << "_" << iteration << ".nc is already open" << std::endl;
#endif
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "--------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationNetCDF<T>::closeFile()
{
	if (true)
	{
		int error;
		error = nc_close(fileId);

	    if (error)
	    {
	        std::cerr << "----- CLbmVisualizationNetCDF<T>::closeFile() -----" << std::endl;
	        std::cerr << "An error occurred while closing parallel netCDF file" << std::endl;
	        std::cerr << nc_strerror(error) << std::endl;
	        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
	        std::cerr << "---------------------------------------------------" << std::endl;

	        exit (EXIT_FAILURE);
	    }
	} else {
        std::cerr << "----- CLbmVisualizationNetCDF<T>::closeFile() -----" << std::endl;
        std::cerr << "There is no open netCDF file to close" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "---------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationNetCDF<T>::defineData()
{
	if (true)
	{
		const std::string DIM_UNIT("meters");
		const std::string FLAG_UNIT("flag");
		const std::string DENSITY_UNIT("kg/m^3");
		const std::string VELOCITY_UNIT("m/s");
		int dimIds[3];

		/*
		 * Due to the strange fact that netCDF takes the first coordinate as
		 * the fastest running, dimensions have to be inverted. On the other
		 * hand, since that's exactly the way how our data is stored, no final
		 * matrix inversion is necessary.
		 */
		nc_def_dim(fileId, "z", solver->getDomain()->getSize()[2], &dimIds[0]);
		nc_def_dim(fileId, "y", solver->getDomain()->getSize()[1], &dimIds[1]);
		nc_def_dim(fileId, "x", solver->getDomain()->getSize()[0], &dimIds[2]);

		nc_def_var(fileId, "z", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 1, &dimIds[0], &dimVarIds[0]);
		nc_def_var(fileId, "y", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 1, &dimIds[1], &dimVarIds[1]);
		nc_def_var(fileId, "x", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 1, &dimIds[2], &dimVarIds[2]);
		nc_def_var(fileId, "flags", NC_INT, 3, dimIds, &flagsVarId);
		nc_def_var(fileId, "densities", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 3, dimIds, &densitiesVarId);

		nc_put_att_text(fileId, dimVarIds[0], "units", DIM_UNIT.length(), DIM_UNIT.c_str());
		nc_put_att_text(fileId, dimVarIds[1], "units", DIM_UNIT.length(), DIM_UNIT.c_str());
		nc_put_att_text(fileId, dimVarIds[2], "units", DIM_UNIT.length(), DIM_UNIT.c_str());
		nc_put_att_text(fileId, flagsVarId, "units", FLAG_UNIT.length(), FLAG_UNIT.c_str());
		nc_put_att_text(fileId, densitiesVarId, "units", DENSITY_UNIT.length(), DENSITY_UNIT.c_str());
		nc_put_att_text(fileId, velocitiesVarId, "units", VELOCITY_UNIT.length(), VELOCITY_UNIT.c_str());

		nc_enddef(fileId);
	} else {
        std::cerr << "----- CLbmVisualizationNetCDF<T>::defineData() -----" << std::endl;
        std::cerr << "There is no open NetCDF file to write dataset" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "----------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationNetCDF<T>::writeData()
{
	if (true)
	{
		T x[solver->getDomain()->getSize()[0]];
		T y[solver->getDomain()->getSize()[1]];
		T z[solver->getDomain()->getSize()[2]];

		for (int i = 0; i < solver->getDomain()->getSize()[0]; i++)
		{
			x[i] = ((T)(solver->getDomain()->getOrigin()[0] + i) + (T)0.5) * solver->getDomain()->getLength()[0] / (T)solver->getDomain()->getSize()[0];
		}
		for (int i = 0; i < solver->getDomain()->getSize()[1]; i++)
		{
			y[i] = ((T)(solver->getDomain()->getOrigin()[1] + i) + (T)0.5) * solver->getDomain()->getLength()[1] / (T)solver->getDomain()->getSize()[1];
		}
		for (int i = 0; i < solver->getDomain()->getSize()[2]; i++)
		{
			z[i] = ((T)(solver->getDomain()->getOrigin()[2] + i) + (T)0.5) * solver->getDomain()->getLength()[2] / (T)solver->getDomain()->getSize()[2];
		}

		solver->getFlags(flags);
		solver->getDensities(densities);
		solver->getVelocities(velocities);

		if (typeid(T) == typeid(double))
		{
			nc_put_var_double(fileId, dimVarIds[2], (double*)x);
			nc_put_var_double(fileId, dimVarIds[1], (double*)y);
			nc_put_var_double(fileId, dimVarIds[0], (double*)z);
			nc_put_var_double(fileId, densitiesVarId, (double*)densities);
		} else if(typeid(T) == typeid(float)) {
			nc_put_var_float(fileId, dimVarIds[2], (float*)x);
			nc_put_var_float(fileId, dimVarIds[1], (float*)y);
			nc_put_var_float(fileId, dimVarIds[0], (float*)z);
			nc_put_var_float(fileId, densitiesVarId, (float*)densities);
		} else {
	        std::cerr << "----- CLbmVisualizationNetCDF<T>::writeData() -----" << std::endl;
	        std::cerr << "A datatype was specified (neither float nor double) which is not supported" << std::endl;
	        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
	        std::cerr << "----------------------------------------------------" << std::endl;

	        exit (EXIT_FAILURE);
		}

#ifdef PAR_NETCDF
		nc_var_par_access(fileId, flagsVarId, NC_COLLECTIVE);
#else
		nc_put_var_int(fileId, flagsVarId, (int*)flags);
#endif

	} else {
        std::cerr << "----- CLbmVisualizationNetCDF<T>::writeData() -----" << std::endl;
        std::cerr << "There is no open NetCDF file to write dataset" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "----------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationNetCDF<T>::render(int iteration)
{
	if(iteration % visualizationRate == 0)
	{
		openFile(iteration);
		defineData();
		writeData();
		closeFile();
	}
}

template class CLbmVisualizationNetCDF<double>;
template class CLbmVisualizationNetCDF<float>;
