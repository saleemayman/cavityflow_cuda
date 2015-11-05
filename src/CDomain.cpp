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
    length = CVector<3, T>((T)1, (T)1, (T)1);
}

template <class T>
CDomain<T>::CDomain(int id, CVector<3, int> size, CVector<3, int> origin, CVector<3, T> length) :
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
CVector<3, T> CDomain<T>::getLengthWithHalo() const
{
	CVector<3, T> lengthWithHalo(
			length.data[0] * ((T)(size.data[0] + 2) / (T)size.data[0]),
			length.data[1] * ((T)(size.data[1] + 2) / (T)size.data[1]),
			length.data[2] * ((T)(size.data[2] + 2) / (T)size.data[2]));

    return lengthWithHalo;
}

template <class T>
int CDomain<T>::getNumOfCells() const
{
    return (size.data[0] * size.data[1] * size.data[2]);
}

template <class T>
int CDomain<T>::getNumOfCellsWithHalo() const
{
    return ((size.data[0] + 2) * (size.data[1] + 2) * (size.data[2] + 2));
}

template <class T>
CVector<3, int> CDomain<T>::getOrigin() const
{
    return origin;
}

template <class T>
CVector<3, int> CDomain<T>::getOriginWithHalo() const
{
	CVector<3, int> originWithHalo(
			origin.data[0] - 1,
			origin.data[1] - 1,
			origin.data[2] - 1);

	return originWithHalo;
}

template <class T>
CVector<3, int> CDomain<T>::getSize() const
{
    return size;
}

template <class T>
CVector<3, int> CDomain<T>::getSizeWithHalo() const
{
	CVector<3, int> sizeWithHalo(
			size.data[0] + 2,
			size.data[1] + 2,
			size.data[2] + 2);

	return sizeWithHalo;
}

template class CDomain<double>;
template class CDomain<float>;
