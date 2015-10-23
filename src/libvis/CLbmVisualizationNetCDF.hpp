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

#include <netcdf.h>
// #include <netcdf_par.h>

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

	std::string filePath;
	int fileId;
	int dimVarIds[3];
	int flagsVarId, densitiesVarId, velocitiesVarId;

	void openFile(int iteration);
    void closeFile();
	void defineData();
	void writeData();

public:
	CLbmVisualizationNetCDF(int id, int visualizationRate, CLbmSolver<T>* solver, std::string filePath);

	void render(int iteration);
};

#endif
