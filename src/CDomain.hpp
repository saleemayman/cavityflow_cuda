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

#ifndef CDOMAIN_HPP
#define CDOMAIN_HPP

#include "libmath/CVector.hpp"

/*
 * Class CDomain contains meta data related to a subdomain of simulation grid.
 * It does not contain the domain data itself like physical values (densities,
 * velocities, etc.) or values which are specific for a certain cell of the
 * domain.
 */
template<typename T>
class CDomain
{
private:
    int id;                // Unique id of the sub/domain
    CVector<3,T> length;   // Length of the domain in each direction
    CVector<3,int> origin; // Origin of the domain points if it is part of a bigger domain
    CVector<3,int> size;   // Size of the domain in each direction

public:
    CDomain(int id, CVector<3,int> size, CVector<3,int> originCell, CVector<3,T> length);
    CDomain(int id, CVector<3,int> size);
    ~CDomain();

    int getId() const;
    CVector<3,T> getLength() const;
    CVector<3,T> getLengthWithHalo() const;
    int getNumOfCells() const;
    int getNumOfCellsWithHalo() const;
    int getNumOfXFaceCells() const;
    int getNumOfXFaceCellsWithHalo() const;
    int getNumOfYFaceCells() const;
    int getNumOfYFaceCellsWithHalo() const;
    int getNumOfZFaceCells() const;
    int getNumOfZFaceCellsWithHalo() const;
    CVector<3,int> getOrigin() const;
    CVector<3,int> getSize() const;
    CVector<3,int> getSizeWithHalo() const;
};

#endif
