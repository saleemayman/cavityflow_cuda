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

#include "CVector.hpp"

#include <cassert>

#include "CMath.hpp"

template <class T>
CVector<4, T>::CVector()
{
    setZero();
}

template <class T>
CVector<4, T>::CVector(const T x0, const T x1, const T x2, const T x3)
{
    data[0] = x0;
    data[1] = x1;
    data[2] = x2;
    data[3] = x3;
}

template <class T>
CVector<4, T>::CVector(const T x)
{
    data[0] = x;
    data[1] = x;
    data[2] = x;
    data[3] = x;
}

template <class T>
CVector<4, T>::CVector(const CVector<4, T>& v)
{
    data[0] = v.data[0];
    data[1] = v.data[1];
    data[2] = v.data[2];
    data[3] = v.data[3];
}

template <class T>
CVector<4, T>::CVector(const T v[4])
{
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    data[3] = v[3];
}


template <class T>
void CVector<4, T>::setZero()
{
    data[0] = T(0);
    data[1] = T(0);
    data[2] = T(0);
    data[3] = T(0);
}

template <class T>
int CVector<4, T>::getGlobalIdx(CVector<4, T>& size)
{
    assert(data[0] < size[0]);
    assert(data[1] < size[1]);
    assert(data[2] < size[2]);
    assert(data[3] < size[3]);

    return (data[3] * size[1] * size[2] * size[3] + data[2] * size[1] * size[2] + data[1] * size[0] + data[0]);
}

template <class T>
T CVector<4, T>::length()
{
    return CMath<T>::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3]);
}

template <class T>
T CVector<4, T>::length2()
{
    return data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3];
}

template <class T>
CVector<4, T>& CVector<4, T>::operator=(const T a[4])
{
    data[0] = a[0];
    data[1] = a[1];
    data[2] = a[2];
    data[3] = a[3];
    return *this;
}

template <class T>
CVector<4, T>& CVector<4, T>::operator=(CVector<4, T> const & a)
{
    data[0] = a.data[0];
    data[1] = a.data[1];
    data[2] = a.data[2];
    data[3] = a.data[3];
    return *this;
}

template <class T>
CVector<4, T>& CVector<4, T>::operator=(CVector<3, T> const & a)
{
    data[0] = a.data[0];
    data[1] = a.data[1];
    data[2] = a.data[2];
    data[3] = 1;
    return *this;
};

template <class T>
CVector<4, T> CVector<4, T>::operator+(const T a)
{
    return CVector<4, T>(data[0]+a, data[1]+a, data[2]+a, data[3]+a);
}

template <class T>
CVector<4, T> CVector<4, T>::operator-(const T a)
{
    return CVector<4, T>(data[0]-a, data[1]-a, data[2]-a, data[3]-a);
}

template <class T>
CVector<4, T> CVector<4, T>::operator*(const T a)
{
    return CVector<4, T>(data[0]*a, data[1]*a, data[2]*a, data[3]*a);
}

template <class T>
CVector<4, T> CVector<4, T>::operator/(const T a)
{
    return CVector<4, T>(data[0]/a, data[1]/a, data[2]/a, data[3]/a);
}

template <class T>
CVector<4, T>&  CVector<4, T>::operator+=(const T a)
{
    data[0] += a;
    data[1] += a;
    data[2] += a;
    data[3] += a;
    return *this;
}

template <class T>
CVector<4, T>& CVector<4, T>::operator-=(const T a)
{
    data[0] -= a;
    data[1] -= a;
    data[2] -= a;
    data[3] -= a;
    return *this;
}

template <class T>
CVector<4, T>& CVector<4, T>::operator*=(const T a)
{
    data[0] *= a;
    data[1] *= a;
    data[2] *= a;
    data[3] *= a;
    return *this;
}

template <class T>
CVector<4, T>& CVector<4, T>::operator/=(const T a)
{
    data[0] /= a;
    data[1] /= a;
    data[2] /= a;
    data[3] /= a;
    return *this;
}

template <class T>
CVector<4, T> CVector<4, T>::operator+(const CVector<4, T>& v)
{
    return CVector<4, T>(data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2], data[3]+v.data[3]);
}

template <class T>
CVector<4, T> CVector<4, T>::operator-(const CVector<4, T>& v)
{
    return CVector<4, T>(data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2], data[3]-v.data[3]);
}

template <class T>
CVector<4, T> CVector<4, T>::operator*(const CVector<4, T>& v)
{
    return CVector<4, T>(data[0]*v.data[0], data[1]*v.data[1], data[2]*v.data[2], data[3]*v.data[3]);
}

template <class T>
CVector<4, T> CVector<4, T>::operator/(const CVector<4, T>& v)
{
    return CVector<4, T>(data[0]/v.data[0], data[1]/v.data[1], data[2]/v.data[2], data[3]/v.data[3]);
}

template <class T>
CVector<4, T>& CVector<4, T>::operator+=(const CVector<4, T>& v)
{
    data[0] += v.data[0];
    data[1] += v.data[1];
    data[2] += v.data[2];
    data[3] += v.data[3];
    return *this;
}

template <class T>
CVector<4, T>& CVector<4, T>::operator-=(const CVector<4, T>& v)
{
    data[0] -= v.data[0];
    data[1] -= v.data[1];
    data[2] -= v.data[2];
    data[3] -= v.data[3];
    return *this;
}

template <class T>
bool CVector<4, T>::operator==(const CVector<4, T>& v)
{
    return bool(data[0] == v.data[0] && data[1] == v.data[1] && data[2] == v.data[2] && data[3] == v.data[3]);
}

template <class T>
bool CVector<4, T>::operator!=(const CVector<4, T>& v)
{
    return bool(data[0] != v.data[0] || data[1] != v.data[1] || data[2] != v.data[2] || data[3] != v.data[3]);
}

template <class T>
T& CVector<4, T>::operator[](const int i)
{
    return data[i];
}

template <class T>
std::ostream& operator<<(std::ostream& co, const CVector<4, T>& v)
{
    return co << "[" << v.data[0] << ", " << v.data[1] << ", " << v.data[2] << ", " << v.data[3] << "]";
}

template class CVector<4, double>;
template class CVector<4, float>;
template class CVector<4, int>;

template std::ostream& operator<<(std::ostream& co, const CVector<4, double>& v);
template std::ostream& operator<<(std::ostream& co, const CVector<4, float>& v);
template std::ostream& operator<<(std::ostream& co, const CVector<4, int>& v);
