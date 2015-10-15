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

#ifndef CLBMVISUALIZATION_HPP
#define CLBMVISUALIZATION_HPP

#include "../CLbmSolver.hpp"

template <typename T>
class CLbmVisualization
{
protected:
	int id;
	Flag* flags;
	T* densities;
	T* velocities;
	CLbmSolver<T>* solver;

public:
	CLbmVisualization(int id, CLbmSolver<T>* solver) :
		id(id), solver(solver)
	{
		flags = new Flag[this->solver->getDomain()->getNumOfCells()];
		densities = new T[this->solver->getDomain()->getNumOfCells()];
		velocities = new T[3 * this->solver->getDomain()->getNumOfCells()];
	}

    virtual ~CLbmVisualization()
    {
		delete[] velocities;
		delete[] densities;
		delete[] flags;
    };

	virtual void render(int iteration = -1) = 0;
};

#endif
