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

#ifndef __CVECTOR_HH
    #error "dont include CVector4.hpp directly!"
#endif

#ifndef __CVECTOR4_HH
#define __CVECTOR4_HH

#include <iostream>

/**
 * \brief    4D Vector handler
 */
template <typename T>
class CVector<4,T>
{
public:
    T data[4];        ///< vector data

    CVector();
    CVector(const T x0, const T x1, const T x2, const T x3);
    CVector(const T x);
    CVector(const T v[4]);

    void setZero();
    T length();
    T length2();

    CVector<4,T>& operator=(const T a[4]);
    CVector<4,T>& operator=(CVector<4,T> const & a);
    CVector<4,T>& operator=(CVector<3,T> const & a);
    CVector<4,T> operator+(const T a);
    CVector<4,T> operator-(const T a);
    CVector<4,T> operator*(const T a);
    CVector<4,T> operator/(const T a);
    CVector<4,T>& operator+=(const T a);
    CVector<4,T>& operator-=(const T a);
    CVector<4,T>& operator*=(const T a);
    CVector<4,T>& operator/=(const T a);
    CVector<4,T> operator+(const CVector<4,T> &v);
    CVector<4,T> operator-(const CVector<4,T> &v);
    CVector<4,T> operator*(const CVector<4,T> &v);
    CVector<4,T> operator/(const CVector<4,T> &v);
    CVector<4,T>& operator+=(const CVector<4,T> &v);
    CVector<4,T>& operator-=(const CVector<4,T> &v);
    bool operator==(const CVector<4,T> &v);
    bool operator!=(const CVector<4,T> &v);
    T& operator[](const int i);
};

template <class T>
std::ostream& operator<<(std::ostream &co, CVector<4,T> &v);

#endif
