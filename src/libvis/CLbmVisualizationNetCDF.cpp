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
#ifdef PAR_NETCDF
        CVector<3, int> numOfSubdomains,
#endif
        std::string filePath) :
        CLbmVisualization<T>(id, visualizationRate, solver),
#ifdef PAR_NETCDF
        numOfSubdomains(numOfSubdomains),
#endif
        filePath(filePath)
{
}

template <class T>
void CLbmVisualizationNetCDF<T>::openFile(int iteration)
{
    int error;
    std::stringstream fileName;
#ifdef PAR_NETCDF
    fileName << filePath << "/visualization_" << iteration << ".nc";
    error = nc_create_par(fileName.str().c_str(), NC_CLOBBER | NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &fileId);
#else
    fileName << filePath << "/visualization_" << id << "_" << iteration << ".nc";
    error = nc_create(fileName.str().c_str(), NC_CLOBBER, &fileId);
#endif

    if (error)
    {
        std::cerr << "----- CLbmVisualizationNetCDF<T>::openFile() -----" << std::endl;
#ifdef PAR_NETCDF
        std::cerr << "An error occurred while opening parallel netCDF file \"" << fileName.str() << "\"" << std::endl;std::cerr << nc_strerror(error) << std::endl;
#else
        std::cerr << "An error occurred while opening serial netCDF file \"" << fileName.str() << "\"" << std::endl;std::cerr << nc_strerror(error) << std::endl;
#endif
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "--------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
    }
}

template <class T>
void CLbmVisualizationNetCDF<T>::closeFile()
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
}

template <class T>
void CLbmVisualizationNetCDF<T>::defineData()
{
    const std::string DIM_UNIT("meters");
    const std::string FLAG_UNIT("flag");
    const std::string DENSITY_UNIT("kg/m^3");
    const std::string VELOCITY_UNIT("m/s");
    int dimIds[3];
    int oldFill;

    /*
     * Due to the strange fact that netCDF takes the first coordinate as
     * the fastest running, dimensions have to be inverted. On the other
     * hand, since that's exactly the way how our data is stored, no final
     * matrix inversion is necessary.
     */
#ifdef PAR_NETCDF
    nc_def_dim(fileId, "x", numOfSubdomains[0] * solver->getDomain()->getSize()[0], &dimIds[2]);
    nc_def_dim(fileId, "y", numOfSubdomains[1] * solver->getDomain()->getSize()[1], &dimIds[1]);
    nc_def_dim(fileId, "z", numOfSubdomains[2] * solver->getDomain()->getSize()[2], &dimIds[0]);
#else
    nc_def_dim(fileId, "x", solver->getDomain()->getSize()[0], &dimIds[2]);
    nc_def_dim(fileId, "y", solver->getDomain()->getSize()[1], &dimIds[1]);
    nc_def_dim(fileId, "z", solver->getDomain()->getSize()[2], &dimIds[0]);
#endif

    nc_def_var(fileId, "x", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 1, &dimIds[2], &dimVarIds[2]);
    nc_def_var(fileId, "y", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 1, &dimIds[1], &dimVarIds[1]);
    nc_def_var(fileId, "z", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 1, &dimIds[0], &dimVarIds[0]);
    nc_def_var(fileId, "flags", NC_INT, 3, dimIds, &flagsVarId);
    nc_def_var(fileId, "densities", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 3, dimIds, &densitiesVarId);
    nc_def_var(fileId, "velocities x", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 3, dimIds, &velocitiesVarId[0]);
    nc_def_var(fileId, "velocities y", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 3, dimIds, &velocitiesVarId[1]);
    nc_def_var(fileId, "velocities z", ((typeid(T) == typeid(float)) ? NC_FLOAT : NC_DOUBLE), 3, dimIds, &velocitiesVarId[2]);

    nc_put_att_text(fileId, dimVarIds[0], "units", DIM_UNIT.length(), DIM_UNIT.c_str());
    nc_put_att_text(fileId, dimVarIds[1], "units", DIM_UNIT.length(), DIM_UNIT.c_str());
    nc_put_att_text(fileId, dimVarIds[2], "units", DIM_UNIT.length(), DIM_UNIT.c_str());
    nc_put_att_text(fileId, flagsVarId, "units", FLAG_UNIT.length(), FLAG_UNIT.c_str());
    nc_put_att_text(fileId, densitiesVarId, "units", DENSITY_UNIT.length(), DENSITY_UNIT.c_str());
    nc_put_att_text(fileId, velocitiesVarId[0], "units", VELOCITY_UNIT.length(), VELOCITY_UNIT.c_str());
    nc_put_att_text(fileId, velocitiesVarId[1], "units", VELOCITY_UNIT.length(), VELOCITY_UNIT.c_str());
    nc_put_att_text(fileId, velocitiesVarId[2], "units", VELOCITY_UNIT.length(), VELOCITY_UNIT.c_str());

    nc_set_fill(fileId, NC_NOFILL, &oldFill);

    nc_enddef(fileId);

#ifdef PAR_NETCDF
    nc_var_par_access(fileId, flagsVarId, NC_INDEPENDENT);
    nc_var_par_access(fileId, densitiesVarId, NC_INDEPENDENT);
    nc_var_par_access(fileId, velocitiesVarId[0], NC_INDEPENDENT);
    nc_var_par_access(fileId, velocitiesVarId[1], NC_INDEPENDENT);
    nc_var_par_access(fileId, velocitiesVarId[2], NC_INDEPENDENT);
#endif
}

template <class T>
void CLbmVisualizationNetCDF<T>::writeData()
{
    if (id == 0)
    {
        T x[solver->getDomain()->getSize()[0]];
        T y[solver->getDomain()->getSize()[1]];
        T z[solver->getDomain()->getSize()[2]];

#ifdef PAR_NETCDF
        for (int i = 0; i < numOfSubdomains[0] * solver->getDomain()->getSize()[0]; i++)
#else
        for (int i = 0; i < solver->getDomain()->getSize()[0]; i++)
#endif
        {
            x[i] = ((T)i + (T)0.5) * solver->getDomain()->getLength()[0] / (T)solver->getDomain()->getSize()[0];
        }
#ifdef PAR_NETCDF
        for (int i = 0; i < numOfSubdomains[1] * solver->getDomain()->getSize()[1]; i++)
#else
        for (int i = 0; i < solver->getDomain()->getSize()[1]; i++)
#endif
        {
            y[i] = ((T)i + (T)0.5) * solver->getDomain()->getLength()[1] / (T)solver->getDomain()->getSize()[1];
        }
#ifdef PAR_NETCDF
        for (int i = 0; i < numOfSubdomains[2] * solver->getDomain()->getSize()[2]; i++)
#else
        for (int i = 0; i < solver->getDomain()->getSize()[2]; i++)
#endif
        {
            z[i] = ((T)i + (T)0.5) * solver->getDomain()->getLength()[2] / (T)solver->getDomain()->getSize()[2];
        }

        if (typeid(T) == typeid(double))
        {
            nc_put_var_double(fileId, dimVarIds[2], (double*)x);
            nc_put_var_double(fileId, dimVarIds[1], (double*)y);
            nc_put_var_double(fileId, dimVarIds[0], (double*)z);
        } else if(typeid(T) == typeid(float)) {
            nc_put_var_float(fileId, dimVarIds[2], (float*)x);
            nc_put_var_float(fileId, dimVarIds[1], (float*)y);
            nc_put_var_float(fileId, dimVarIds[0], (float*)z);
        } else {
            std::cerr << "----- CLbmVisualizationNetCDF<T>::writeData() -----" << std::endl;
            std::cerr << "A datatype was specified (neither float nor double) which is not supported" << std::endl;
            std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
            std::cerr << "----------------------------------------------------" << std::endl;

            exit (EXIT_FAILURE);
        }
    }

    solver->getFlags(flags);
    solver->getDensities(densities);
    solver->getVelocities(velocities);

    size_t start[3] = {solver->getDomain()->getOrigin()[0], solver->getDomain()->getOrigin()[1], solver->getDomain()->getOrigin()[2]};
    size_t count[3] = {solver->getDomain()->getSize()[2], solver->getDomain()->getSize()[1], solver->getDomain()->getSize()[0]};

    nc_put_vara_int(fileId, flagsVarId, start, count, (int*)flags);

    if (typeid(T) == typeid(double))
    {
        nc_put_vara_double(fileId, densitiesVarId, start, count, (double*)densities);
        nc_put_vara_double(fileId, velocitiesVarId[0], start, count, (double*)&(velocities[0]));
        nc_put_vara_double(fileId, velocitiesVarId[1], start, count, (double*)&(velocities[solver->getDomain()->getNumOfCells()]));
        nc_put_vara_double(fileId, velocitiesVarId[2], start, count, (double*)&(velocities[2 * solver->getDomain()->getNumOfCells()]));
    } else if(typeid(T) == typeid(float)) {
        nc_put_vara_float(fileId, densitiesVarId, start, count, (float*)densities);
        nc_put_vara_float(fileId, velocitiesVarId[0], start, count, (float*)&(velocities[0]));
        nc_put_vara_float(fileId, velocitiesVarId[1], start, count, (float*)&(velocities[solver->getDomain()->getNumOfCells()]));
        nc_put_vara_float(fileId, velocitiesVarId[2], start, count, (float*)&(velocities[2 * solver->getDomain()->getNumOfCells()]));
    } else {
        std::cerr << "----- CLbmVisualizationNetCDF<T>::writeData() -----" << std::endl;
        std::cerr << "A datatype was specified (neither float nor double) which is not supported" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "----------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
    }
}

template <class T>
void CLbmVisualizationNetCDF<T>::render(int iteration)
{
    if (iteration % visualizationRate == 0)
    {
        openFile(iteration);
        defineData();
        writeData();
        closeFile();
    }
}

template class CLbmVisualizationNetCDF<double>;
template class CLbmVisualizationNetCDF<float>;
