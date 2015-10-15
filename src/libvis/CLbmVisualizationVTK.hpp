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

#ifndef CLBMVISUALIZATIONVTK_HPP
#define CLBMVISUALIZATIONVTK_HPP

#include <fstream>
#include <sstream>

#include "CLbmVisualization.hpp"

template <typename T>
class CLbmVisualizationVTK : public virtual CLbmVisualization<T>
{
private:
    using CLbmVisualization<T>::id;
    using CLbmVisualization<T>::flags;
    using CLbmVisualization<T>::densities;
    using CLbmVisualization<T>::velocities;
    using CLbmVisualization<T>::solver;

    std::string filePath;
    std::ofstream file;

    void openFile(int iteration);
    void closeFile();
	void writeHeader();
	void writeDataset();
	void writeFlags();
	void writeDensities();
	void writeVelocities();

public:
	CLbmVisualizationVTK(int id, CLbmSolver<T>* solver, std::string filePath);

	void render(int iteration);
};

#endif
