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

#include "CLbmVisualizationVTK.hpp"

template <class T>
CLbmVisualizationVTK<T>::CLbmVisualizationVTK(
		int id,
		CLbmSolver<T>* solver,
		std::string filePath) :
		CLbmVisualization<T>(id, solver),
		filePath(filePath)
{
}

template <class T>
void CLbmVisualizationVTK<T>::openFile(int iteration)
{
	if (!file.is_open())
	{
	    std::stringstream fileName;
	    fileName << filePath << "/visualization_" << id << "_" << iteration << ".vtk";
	    std::ofstream benchmarkFile(fileName.str().c_str(), std::ios::out);
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::closeFile()
{
	if (file.is_open())
	{
	    file.close();
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeHeader()
{
	if (file.is_open())
	{
		file << "# vtk DataFile Version 3.0\n";
		file << "Heterogenous (MPI/OpenMPI) and hybrid (CPU/GPU) LBM simulation on rank " << id << "\n";
		file << "ASCII" << std::endl;
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeDataset()
{
	if (file.is_open())
	{
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeFlags()
{
	if (file.is_open())
	{
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeDensities()
{
	if (file.is_open())
	{
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeVelocities()
{
	if (file.is_open())
	{
	} else {
		/*
		 * TODO
		 */
	}
}

template <class T>
void CLbmVisualizationVTK<T>::render(int iteration)
{
	openFile(iteration);
	writeHeader();
	writeDataset();
	writeFlags();
	writeDensities();
	writeVelocities();
	closeFile();
}

template class CLbmVisualizationVTK<double>;
template class CLbmVisualizationVTK<float>;
