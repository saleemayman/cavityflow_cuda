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

#include "CDomain.hpp"

template <class T>
CDomain<T>::CDomain(int id, CVector<3, int> size) :
        id(id), size(size)
{
    origin = CVector<3, int>(0, 0, 0);
    length = CVector<3, T>((T)0.05, (T)0.05, (T)0.05);
}

template <class T>
CDomain<T>::CDomain(int id, CVector<3, int> size, CVector<3, int> originCell, CVector<3, T> length) :
        id(id), size(size), origin(origin), length(length)
{
}

template <class T>
CDomain<T>::~CDomain()
{
}

template <class T>
int CDomain<T>::getId() const
{
    return id;
}

template <class T>
CVector<3, T> CDomain<T>::getLength() const
{
    return length;
}

template <class T>
int CDomain<T>::getNumOfCells() const
{
    return (size.data[0] * size.data[1] * size.data[2]);
}

template <class T>
int CDomain<T>::getNumOfXFaceCells() const
{
    return (size.data[1] * size.data[2]);
}

template <class T>
int CDomain<T>::getNumOfYFaceCells() const
{
    return (size.data[0] * size.data[2]);
}

template <class T>
int CDomain<T>::getNumOfZFaceCells() const
{
    return (size.data[0] * size.data[1]);
}

template <class T>
CVector<3, int> CDomain<T>::getOrigin() const
{
    return origin;
}

template <class T>
CVector<3, int> CDomain<T>::getSize() const
{
    return size;
}

template class CDomain<double>;
template class CDomain<float>;
