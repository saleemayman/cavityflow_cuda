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
CVector<2, T>::CVector()
{
    setZero();
}

template <class T>
CVector<2, T>::CVector(const T x0, const T x1)
{
    data[0] = x0;
    data[1] = x1;
}

template <class T>
CVector<2, T>:: CVector(const T x)
{
    data[0] = x;
    data[1] = x;
}

template <class T>
CVector<2, T>::CVector(const T v[2])
{
    data[0] = v[0];
    data[1] = v[1];
}

template <class T>
void CVector<2, T>::setZero()
{
    data[0] = T(0);
    data[1] = T(0);
}

template <class T>
int CVector<2, T>::getGlobalIdx(CVector<2, T>& size)
{
    assert(data[0] < size[0]);
    assert(data[1] < size[1]);

    return (data[1] * size[0] + data[0]);
}

template <class T>
CVector<2, T> CVector<2, T>::getNormal()
{
    CVector<2, T> v = *this;
    T inv_length = (T) 1/CMath<T>::sqrt(v[0]*v[0] + v[1]*v[1]);
    v[0] *= inv_length;
    v[1] *= inv_length;
    return v;
}

template <class T>
T CVector<2, T>::dotProd(const CVector<2, T>& v)
{
    return v.data[0]*data[0] + v.data[1]*data[1];
}

template <class T>
T CVector<2, T>::crossProd(CVector<2, T> &a)
{
    return data[0]*a.data[1] - data[1]*a.data[0];
}


template <class T>
T CVector<2, T>::elements()
{
    return (data[0]*data[1]);
}

template <class T>
T CVector<2, T>::length()
{
    return CMath<T>::sqrt(data[0]*data[0] + data[1]*data[1]);
}

template <class T>
T CVector<2, T>::length2()
{
    return data[0]*data[0] + data[1]*data[1];
}

template <class T>
T CVector<2, T>::dist2(const CVector<2, T>& v)
{
    CVector<2, T> d = CVector<2, T>(v.data[0] - data[0], v.data[1] - data[1]);
    return (d[0]*d[0] + d[1]*d[1]);
}

template <class T>
T CVector<2, T>::dist(const CVector<2, T>& v)
{
    CVector<2, T> d = CVector<2, T>(v.data[0] - data[0], v.data[1] - data[1]);
    return CMath<T>::sqrt(d[0]*d[0] + d[1]*d[1]);
}

template <class T>
void CVector<2, T>::normalize()
{
    T il = 1/length();
    data[0] *= il;
    data[1] *= il;
}

template <class T>
void CVector<2, T>::clamp1_1()
{
    data[0] = (data[0] < -1 ? -1 : (data[0] > 1 ? 1 : data[0]));
    data[1] = (data[1] < -1 ? -1 : (data[1] > 1 ? 1 : data[1]));
}

template <class T>
CVector<2, T>& CVector<2, T>::operator=(const T a[2])
{
    data[0] = a[0];
    data[1] = a[1];
    return *this;
}

template <class T>
CVector<2, T>& CVector<2, T>::operator=(CVector<2, T> const& a)
{
    data[0] = a.data[0];
    data[1] = a.data[1];
    return *this;
}

template <class T>
CVector<2, T> CVector<2, T>::operator-()
{
    return CVector<2, T>(-data[0], -data[1]);
}

template <class T>
CVector<2, T> CVector<2, T>::operator+(const T a)
{
    return CVector<2, T>(data[0]+a, data[1]+a);
}

template <class T>
CVector<2, T> CVector<2, T>::operator-(const T a)
{
    return CVector<2, T>(data[0]-a, data[1]-a);
}

template <class T>
CVector<2, T> CVector<2, T>::operator*(const T a)
{
    return CVector<2, T>(data[0]*a, data[1]*a);
}

template <class T>
CVector<2, T> CVector<2, T>::operator/(const T a)
{
    return CVector<2, T>(data[0]/a, data[1]/a);
}

template <class T>
CVector<2, T>& CVector<2, T>::operator+=(const T a)
{
    data[0] += a;
    data[1] += a;
    return *this;
}

template <class T>
CVector<2, T>& CVector<2, T>::operator-=(const T a)
{
    data[0] -= a;
    data[1] -= a;
    return *this;
}

template <class T>
CVector<2, T>& CVector<2, T>::operator*=(const T a)
{
    data[0] *= a; data[1] *= a;
    return *this;
}

template <class T>
CVector<2, T>& CVector<2, T>::operator/=(const T a)
{
    data[0] /= a;
    data[1] /= a;
    return *this;
}

template <class T>
CVector<2, T> CVector<2, T>::operator+(const CVector<2, T>& v)
{
    return CVector<2, T>(data[0]+v.data[0], data[1]+v.data[1]);
}

template <class T>
CVector<2, T> CVector<2, T>::operator-(const CVector<2, T>& v)
{
    return CVector<2, T>(data[0]-v.data[0], data[1]-v.data[1]);
}

template <class T>
CVector<2, T> CVector<2, T>::operator*(const CVector<2, T>& v)
{
    return CVector<2, T>(data[0]*v.data[0], data[1]*v.data[1]);
}

template <class T>
CVector<2, T> CVector<2, T>::operator/(const CVector<2, T>& v)
{
    return CVector<2, T>(data[0]/v.data[0], data[1]/v.data[1]);
}

template <class T>
CVector<2, T>& CVector<2, T>::operator+=(const CVector<2, T>& v)
{
    data[0] += v.data[0];
    data[1] += v.data[1];
    return *this;
}

template <class T>
CVector<2, T>& CVector<2, T>::operator-=(const CVector<2, T>& v)
{
    data[0] -= v.data[0];
    data[1] -= v.data[1];
    return *this;
}

template <class T>
bool CVector<2, T>::operator==(const CVector<2, T>& v)
{
    return bool(data[0] == v.data[0] && data[1] == v.data[1]);
}

template <class T>
bool CVector<2, T>::operator!=(const CVector<2, T>& v)
{
    return bool(data[0] != v.data[0] || data[1] != v.data[1]);
}

template <class T>
T& CVector<2, T>::operator[](const int i)
{
    return data[i];
}

template <class T>
std::ostream& operator<<(std::ostream& co, const CVector<2, T>& v)
{
    return co << "[" << v.data[0] << ", " << v.data[1] << "]";
}

template class CVector<2, double>;
template class CVector<2, float>;
template class CVector<2, int>;

template std::ostream& operator<<(std::ostream& co, const CVector<2, double>& v);
template std::ostream& operator<<(std::ostream& co, const CVector<2, float>& v);
template std::ostream& operator<<(std::ostream& co, const CVector<2, int>& v);
