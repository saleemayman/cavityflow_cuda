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

#ifndef CLBMVISUALIZATIONNETCDF_HPP
#define CLBMVISUALIZATIONNETCDF_HPP

#include "CLbmVisualization.hpp"

#include <fstream>
#include <sstream>

#ifdef PAR_NETCDF
#include <netcdf_par.h>
#endif
#include <netcdf.h>

template <typename T>
class CLbmVisualizationNetCDF : public virtual CLbmVisualization<T>
{
private:
    using CLbmVisualization<T>::id;
    using CLbmVisualization<T>::visualizationRate;
    using CLbmVisualization<T>::flags;
    using CLbmVisualization<T>::densities;
    using CLbmVisualization<T>::velocities;
    using CLbmVisualization<T>::solver;

#ifdef PAR_NETCDF
    CVector<3, int> numOfSubdomains;
#endif
	std::string filePath;

	int fileId;
	int dimVarIds[3];
	int flagsVarId, densitiesVarId;
	int velocitiesVarId[3];

	void openFile(int iteration);
    void closeFile();
	void defineData();
	void writeData();

public:
#ifdef PAR_NETCDF
	CLbmVisualizationNetCDF(int id, int visualizationRate, CLbmSolver<T>* solver, CVector<3, int> numOfSubdomains, std::string filePath);
#else
	CLbmVisualizationNetCDF(int id, int visualizationRate, CLbmSolver<T>* solver, std::string filePath);
#endif

	void render(int iteration);
};

#endif
