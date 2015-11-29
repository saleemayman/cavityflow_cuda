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
CVector<3, T>::CVector()
{
    setZero();
}

template <class T>
CVector<3, T>::CVector(const T x0, const T x1, const T x2)
{
    data[0] = x0;
    data[1] = x1;
    data[2] = x2;
}

template <class T>
CVector<3, T>::CVector(const T x)
{
    data[0] = x;
    data[1] = x;
    data[2] = x;
}

template <class T>
CVector<3, T>::CVector(const CVector<3, T>& v)
{
    data[0] = v.data[0];
    data[1] = v.data[1];
    data[2] = v.data[2];
}

template <class T>
CVector<3, T>::CVector(const T v[3])
{
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
}

template <class T>
void CVector<3, T>::set(const T x, const T y, const T z)
{
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

template <class T>
void CVector<3, T>::setZero()
{
    data[0] = T(0);
    data[1] = T(0);
    data[2] = T(0);
}

template <class T>
int CVector<3, T>::getGlobalIdx(CVector<3, T>& size)
{
    assert(data[0] < size[0]);
    assert(data[1] < size[1]);
    assert(data[2] < size[2]);

    return (data[2] * size[1] * size[2] + data[1] * size[0] + data[0]);
}

template <class T>
CVector<3, T> CVector<3, T>::getNormal()
{
    CVector<3, T> v = *this;
    T inv_length = (T) 1/CMath<T>::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] *= inv_length;
    v[1] *= inv_length;
    v[2] *= inv_length;
    return v;
}

template <class T>
T CVector<3, T>::dotProd(const CVector<3, T>& v)
{
    return v.data[0]*data[0] + v.data[1]*data[1] + v.data[2]*data[2];
}

template <class T>
CVector<3, T> CVector<3, T>::crossProd(CVector<3, T> &a)
{
    return CVector<3, T>(data[1]*a.data[2] - data[2]*a.data[1], data[2]*a.data[0] - data[0]*a.data[2], data[0]*a.data[1] - data[1]*a.data[0]);
}

template <class T>
T CVector<3, T>::dist2(const CVector<3, T>& v)
{
    CVector<3, T> d = CVector<3, T>(v.data[0] - data[0], v.data[1] - data[1], v.data[2] - data[2]);
    return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

template <class T>
T CVector<3, T>::dist(const CVector<3, T>& v)
{
    CVector<3, T> d = CVector<3, T>(v.data[0] - data[0], v.data[1] - data[1], v.data[2] - data[2]);
    return CMath<T>::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

template <class T>
T CVector<3, T>::elements()
{
    return (data[0]*data[1]*data[2]);
}

template <class T>
T CVector<3, T>::length()
{
    return CMath<T>::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
}

template <class T>
T CVector<3, T>::length2()
{
    return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
}

template <class T>
T CVector<3, T>::min()
{
    T min = (data[0] < data[1] ? data[0] : data[1]);
    return (data[2] < min ? data[2] : min);
}

template <class T>
T CVector<3, T>::max()
{
    T max = (data[0] > data[1] ? data[0] : data[1]);
    return (data[2] > max ? data[2] : max);
}

template <class T>
void CVector<3, T>::normalize()
{
    T il = 1/length();
    data[0] *= il;
    data[1] *= il;
    data[2] *= il;
}

template <class T>
void CVector<3, T>::clamp1_1()
{
    data[0] = (data[0] < -1 ? -1 : (data[0] > 1 ? 1 : data[0]));
    data[1] = (data[1] < -1 ? -1 : (data[1] > 1 ? 1 : data[1]));
    data[2] = (data[2] < -1 ? -1 : (data[2] > 1 ? 1 : data[2]));
}

template <class T>
CVector<3, T>& CVector<3, T>::operator=(const T a[3])
{
    data[0] = a[0];
    data[1] = a[1];
    data[2] = a[2];
    return *this;
}

template <class T>
CVector<3, T>& CVector<3, T>::operator=(CVector<3, T> const & a)
{
    data[0] = a.data[0];
    data[1] = a.data[1];
    data[2] = a.data[2];
    return *this;
}

template <class T>
CVector<3, T> CVector<3, T>::operator+(const T a)
{
    return CVector<3, T>(data[0]+a, data[1]+a, data[2]+a);
}

template <class T>
CVector<3, T> CVector<3, T>::operator-(const T a)
{
    return CVector<3, T>(data[0]-a, data[1]-a, data[2]-a);
}

template <class T>
CVector<3, T> CVector<3, T>::operator*(const T a)
{
    return CVector<3, T>(data[0]*a, data[1]*a, data[2]*a);
}

template <class T>
CVector<3, T> CVector<3, T>::operator/(const T a)
{
    return CVector<3, T>(data[0]/a, data[1]/a, data[2]/a);
}

template <class T>
CVector<3, T>& CVector<3, T>::operator+=(const T a)
{
    data[0] += a;
    data[1] += a;
    data[2] += a;
    return *this;
}

template <class T>
CVector<3, T>& CVector<3, T>::operator-=(const T a)
{
    data[0] -= a;
    data[1] -= a;
    data[2] -= a;
    return *this;
}

template <class T>
CVector<3, T>& CVector<3, T>::operator*=(const T a)
{
    data[0] *= a;
    data[1] *= a;
    data[2] *= a;
    return *this;
}

template <class T>
CVector<3, T>& CVector<3, T>::operator/=(const T a)
{
    data[0] /= a;
    data[1] /= a;
    data[2] /= a;
    return *this;
}

template <class T>
CVector<3, T> CVector<3, T>::operator+(const CVector<3, T>& v)
{
    return CVector<3, T>(data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2]);
}

template <class T>
CVector<3, T> CVector<3, T>::operator-(const CVector<3, T>& v)
{
    return CVector<3, T>(data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2]);
}

template <class T>
CVector<3, T> CVector<3, T>::operator*(const CVector<3, T>& v)
{
    return CVector<3, T>(data[0]*v.data[0], data[1]*v.data[1], data[2]*v.data[2]);
}

template <class T>
CVector<3, T> CVector<3, T>::operator/(const CVector<3, T>& v)
{
    return CVector<3, T>(data[0]/v.data[0], data[1]/v.data[1], data[2]/v.data[2]);
}

template <class T>
CVector<3, T>& CVector<3, T>::operator+=(const CVector<3, T>& v)
{
    data[0] += v.data[0];
    data[1] += v.data[1];
    data[2] += v.data[2];
    return *this;
}

template <class T>
CVector<3, T>& CVector<3, T>::operator-=(const CVector<3, T>& v)
{
    data[0] -= v.data[0];
    data[1] -= v.data[1];
    data[2] -= v.data[2];
    return *this;
}

template <class T>
bool CVector<3, T>::operator==(const CVector<3, T>& v)
{
    return bool(data[0] == v.data[0] && data[1] == v.data[1] && data[2] == v.data[2]);
}

template <class T>
bool CVector<3, T>::operator!=(const CVector<3, T>& v) {
    return bool(data[0] != v.data[0] || data[1] != v.data[1] || data[2] != v.data[2]);
}

template <class T>
T& CVector<3, T>::operator[](const int i)
{
    return data[i];
}

template <class T>
std::ostream& operator<<(std::ostream& co, const CVector<3, T>& v)
{
    return co << "[" << v.data[0] << ", " << v.data[1] << ", " << v.data[2] << "]";
}

template class CVector<3, double>;
template class CVector<3, float>;
template class CVector<3, int>;

template std::ostream& operator<<(std::ostream& co, const CVector<3, double>& v);
template std::ostream& operator<<(std::ostream& co, const CVector<3, float>& v);
template std::ostream& operator<<(std::ostream& co, const CVector<3, int>& v);
