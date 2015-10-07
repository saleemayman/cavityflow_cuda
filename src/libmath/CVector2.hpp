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

#ifndef CVECTOR_HPP
	#error "dont include CVector2.hpp directly!"
#endif

#ifndef CVECTOR2_HPP
#define CVECTOR2_HPP

#include <iostream>

/**
 * \brief	2D Vector handler
 */
template <typename T>
class CVector<2, T>
{
public:
	T data[2];

	CVector();
	CVector(const T x0, const T x1);
	CVector(const T x);
	CVector(const T v[2]);

	void setZero();
	CVector<2, T> getNormal();
	T dotProd(const CVector<2, T>& v);
	T crossProd(CVector<2, T> &a);
	T elements();
	T length();
	T length2();
	T dist2(const CVector<2, T>& v);
	T dist(const CVector<2, T>& v);
	void normalize();
	void clamp1_1();

	CVector<2, T>& operator=(const T a[2]);
	CVector<2, T>& operator=(CVector<2, T> const& a);
    CVector<2, T> operator-();
    CVector<2, T> operator+(const T a);
	CVector<2, T> operator-(const T a);
	CVector<2, T> operator*(const T a);
	CVector<2, T> operator/(const T a);
	CVector<2, T>& operator+=(const T a);
	CVector<2, T>& operator-=(const T a);
	CVector<2, T>& operator*=(const T a);
	CVector<2, T>& operator/=(const T a);
	CVector<2, T> operator+(const CVector<2, T>& v);
	CVector<2, T> operator-(const CVector<2, T>& v);
	CVector<2, T> operator*(const CVector<2, T>& v);
	CVector<2, T> operator/(const CVector<2, T>& v);
	CVector<2, T>& operator+=(const CVector<2, T>& v);
	CVector<2, T>& operator-=(const CVector<2, T>& v);
	bool operator==(const CVector<2, T>& v);
	bool operator!=(const CVector<2, T>& v);
	T& operator[](const int i);
};

template <class T>
std::ostream& operator<<(std::ostream& co, const CVector<2, T>& v);

#endif
