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
    #error "dont include CVector3.hpp directly!"
#endif

#ifndef __CVECTOR3_HH
#define __CVECTOR3_HH

#include <iostream>

/**
 * \brief    3D Vector handler
 */
template <typename T>
class CVector<3,T>
{
public:
    T data[3];    ///< vector data

    CVector();
    CVector(const T x0, const T x1, const T x2);
    CVector(const T x);
    CVector(const CVector<3,T> &v);
    CVector(const T v[3]);

    void setZero();
    CVector<3,T> normal();
    T dotProd(const CVector<3,T> &v);
    CVector<3,T> crossProd(CVector<3,T> &a);
    T dist2(const CVector<3,T> &v);
    T dist(const CVector<3,T> &v);
    T elements();
    T length();
    T length2();
    T min();
    T max();
    void normalize();
    void clamp1_1();

    CVector<3,T>& operator=(const T a[3]);
    CVector<3,T>& operator=(CVector<3,T> const & a);
    CVector<3,T> operator+(const T a);
    CVector<3,T> operator-(const T a);
    CVector<3,T> operator*(const T a);
    CVector<3,T> operator/(const T a);
    CVector<3,T>& operator+=(const T a);
    CVector<3,T>& operator-=(const T a);
    CVector<3,T>& operator*=(const T a);
    CVector<3,T>& operator/=(const T a);
    CVector<3,T> operator+(const CVector<3,T> &v);
    CVector<3,T> operator-(const CVector<3,T> &v);
    CVector<3,T> operator*(const CVector<3,T> &v);
    CVector<3,T> operator/(const CVector<3,T> &v);
    CVector<3,T>& operator+=(const CVector<3,T> &v);
    CVector<3,T>& operator-=(const CVector<3,T> &v);
    bool operator==(const CVector<3,T> &v);
    bool operator!=(const CVector<3,T> &v);
    T& operator[](const int i);
};

template <class T>
std::ostream& operator<<(std::ostream &co, CVector<3,T> &v);

#endif
