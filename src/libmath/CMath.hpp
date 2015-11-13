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

#ifndef CMATH_HPP
#define CMATH_HPP

#include <iostream>

/**
 * \brief    math handler to use same function names for different types
 */
template <typename T>
class CMath
{
public:
    T value;    ///< the value itself!

    CMath();
    CMath(T a);

    /// return PI
    static T PI();
    static T MAX();
    static T MIN();
    static T INF();

    static T abs(T a);
    static T pow(T base, T exp);
    static T floor(T a);
    static T ceil(T a);
    static T ceil2(T a);
    static T sqrt(T a);
    static T round(T a);
    static T digits2(T a);
    static T exp(T a);
    static T sin(T a);
    static T cos(T a);
    static T tan(T a);
    static T max(T a, T b);
    static T min(T a, T b);
    static T aton(const char *s);
    static T gcd(T a, T b);

    CMath operator+(const CMath &a);
    CMath operator-(const CMath &a);
    CMath operator*(const CMath a);
    CMath operator/(const CMath &a);
    CMath& operator+=(const CMath &a);
    CMath& operator-=(const CMath &a);
    CMath& operator*=(const CMath &a);
    CMath& operator/=(const CMath &a);
    CMath operator+(const T a);
    CMath operator-(const T a);
    CMath operator*(const T a);
    CMath operator/(const T a);
    CMath& operator+=(const T a);
    CMath& operator-=(const T a);
    CMath& operator*=(const T a);
    CMath& operator/=(const T a);
};

template <typename T>
std::ostream& operator<<(std::ostream& os, CMath<T> &v);

#endif
