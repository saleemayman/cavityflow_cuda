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

#include <typeinfo>

template <class T>
CLbmVisualizationVTK<T>::CLbmVisualizationVTK(
		int id,
		int visualizationRate,
		CLbmSolver<T>* solver,
		std::string filePath) :
		CLbmVisualization<T>(id, visualizationRate, solver),
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
	    file.open(fileName.str().c_str(), std::ios::out);
	} else {
        std::cerr << "----- CLbmVisualizationVTK<T>::openFile() -----" << std::endl;
        std::cerr << "VTK file \"visualization_" << id << "_" << iteration << ".vtk is already open." << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "-----------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationVTK<T>::closeFile()
{
	if (file.is_open())
	{
	    file.close();
	} else {
        std::cerr << "----- CLbmVisualizationVTK<T>::closeFile() -----" << std::endl;
        std::cerr << "There is no open VTK file to close." << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
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
        std::cerr << "----- CLbmVisualizationVTK<T>::writeHeader() -----" << std::endl;
        std::cerr << "There is no open VTK file to write header" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "--------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeDataset()
{
	if (file.is_open())
	{
		file << "DATASET STRUCTURED_GRID\n";
		file << "DIMENSIONS " << (solver->getDomain()->getSize()[0] + 1) << " " << (solver->getDomain()->getSize()[1] + 1) << " " << (solver->getDomain()->getSize()[2] + 1) << "\n";
		file << "POINTS " << ((solver->getDomain()->getSize()[0] + 1) * (solver->getDomain()->getSize()[1] + 1) * (solver->getDomain()->getSize()[2] + 1)) << " " << ((typeid(T) == typeid(float)) ? "float" : "double") << "\n";

		for (int k = 0; k < solver->getDomain()->getSize()[2] + 1; k++) {
			for (int j = 0; j < solver->getDomain()->getSize()[1] + 1; j++) {
				for (int i = 0; i < solver->getDomain()->getSize()[0] + 1; i++) {
					file << ((T)(solver->getDomain()->getOrigin()[0] + i) * solver->getDomain()->getLength()[0] / (T)solver->getDomain()->getSize()[0]) << " " <<
							((T)(solver->getDomain()->getOrigin()[1] + j) * solver->getDomain()->getLength()[1] / (T)solver->getDomain()->getSize()[1]) << " " <<
							((T)(solver->getDomain()->getOrigin()[2] + k) * solver->getDomain()->getLength()[2] / (T)solver->getDomain()->getSize()[2]) <<"\n";
				}
			}
		}
	} else {
        std::cerr << "----- CLbmVisualizationVTK<T>::writeDataset() -----" << std::endl;
        std::cerr << "There is no open VTK file to write dataset" << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "---------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeFlags()
{
	if (file.is_open())
	{
		solver->getFlags(flags);

		file << "SCALARS flags INT 1\n";
		file << "LOOKUP_TABLE default\n";

		for (int i = 0; i < solver->getDomain()->getNumOfCells(); i++)
		{
			file << flags[i] << "\n";
		}
	} else {
        std::cerr << "----- CLbmVisualizationVTK<T>::writeFlags() -----" << std::endl;
        std::cerr << "There is no open VTK file to write flags." << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "-----------------------------------------------  " << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeDensities()
{
	if (file.is_open())
	{
		solver->getDensities(densities);

		file << "SCALARS densities " << ((typeid(T) == typeid(float)) ? "float" : "double") << " 1\n";
		file << "LOOKUP_TABLE default\n";

		for (int i = 0; i < solver->getDomain()->getNumOfCells(); i++)
		{
			file << densities[i] << "\n";
		}
	} else {
        std::cerr << "----- CLbmVisualizationVTK<T>::writeDensities() -----" << std::endl;
        std::cerr << "There is no open VTK file to write densities." << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "-----------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationVTK<T>::writeVelocities()
{
	if (file.is_open())
	{
		solver->getVelocities(velocities);

		T* velocitiesX = velocities;
		T* velocitiesY = velocities + solver->getDomain()->getNumOfCells();
		T* velocitiesZ = velocities + 2 * solver->getDomain()->getNumOfCells();

		file << "VECTORS velocities " << ((typeid(T) == typeid(float)) ? "float" : "double") << "\n";

		for (int i = 0; i < solver->getDomain()->getNumOfCells(); i++)
		{
			file << velocitiesX[i] << " " << velocitiesY[i] << " " << velocitiesZ[i] << "\n";
		}
	} else {
        std::cerr << "----- CLbmVisualizationVTK<T>::writeVelocities() -----" << std::endl;
        std::cerr << "There is no open VTK file to write velocities." << std::endl;
        std::cerr << "EXECUTION WILL BE TERMINATED IMMEDIATELY" << std::endl;
        std::cerr << "------------------------------------------------------" << std::endl;

        exit (EXIT_FAILURE);
	}
}

template <class T>
void CLbmVisualizationVTK<T>::render(int iteration)
{
	if(iteration % visualizationRate == 0)
	{
		openFile(iteration);
		writeHeader();
		writeDataset();
		file << "CELL_DATA " << solver->getDomain()->getNumOfCells() << "\n";
		writeFlags();
		writeDensities();
		writeVelocities();
		closeFile();

	}
}

template class CLbmVisualizationVTK<double>;
template class CLbmVisualizationVTK<float>;
