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

#include "CMath.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>

template <class T>
CMath<T>::CMath()
{
}

template <class T>
CMath<T>::CMath(T a)
{
    value = a;
}

template <class T>
T CMath<T>::PI()
{
    return M_PI;
}

template <class T>
T CMath<T>::MAX()
{
    return std::numeric_limits<T>::max();
}

template <class T>
T CMath<T>::MIN()
{
    return std::numeric_limits<T>::min();
}

template <class T>
T CMath<T>::INF()
{
    return std::numeric_limits<T>::infinity();
}

template <class T>
T CMath<T>::max(T a, T b)
{
    return (a < b ? b : a);
}

template <class T>
T CMath<T>::min(T a, T b)
{
    return (a > b ? b : a);
}

template <class T>
CMath<T> CMath<T>::operator+(const CMath &a)
{
    return CMath(value + a.value);
}

template <class T>
CMath<T> CMath<T>::operator-(const CMath &a)
{
    return CMath(value - a.value);
}

template <class T>
CMath<T> CMath<T>::operator*(const CMath a)
{
    return CMath(value * a.value);
}

template <class T>
CMath<T> CMath<T>::operator/(const CMath &a)
{
    return CMath(value / a.value);
}

template <class T>
CMath<T>& CMath<T>::operator+=(const CMath &a)
{
    value += a.value;
    return *this;
}

template <class T>
CMath<T>& CMath<T>::operator-=(const CMath &a)
{
    value -= a.value;
    return *this;
}

template <class T>
CMath<T>& CMath<T>::operator*=(const CMath &a)
{
    value *= a.value;
    return *this;
}

template <class T>
CMath<T>& CMath<T>::operator/=(const CMath &a)
{
    value /= a.value;
    return *this;
}

template <class T>
CMath<T> CMath<T>::operator+(const T a)
{
    return CMath(value + a);
}

template <class T>
CMath<T> CMath<T>::operator-(const T a)
{
    return CMath(value - a);
}

template <class T>
CMath<T> CMath<T>::operator*(const T a)
{
    return CMath(value * a);
}

template <class T>
CMath<T> CMath<T>::operator/(const T a)
{
    return CMath(value / a);
}

template <class T>
CMath<T>& CMath<T>::operator+=(const T a)
{
    value += a;
    return *this;
}

template <class T>
CMath<T>& CMath<T>::operator-=(const T a)
{
    value -= a;
    return *this;
}

template <class T>
CMath<T>& CMath<T>::operator*=(const T a)
{
    value *= a;
    return *this;
}

template <class T>
CMath<T>& CMath<T>::operator/=(const T a)
{
    value /= a;
    return *this;
}

template <>
float CMath<float>::abs(float a)
{
    return fabsf(a);
};
template <>
double CMath<double>::abs(double a)
{
    return fabs(a);
}
template <>
int CMath<int>::abs(int a)
{
    return std::abs(a);
}

template <>
float CMath<float>::pow(float base, float exp)
{
    return powf(base, exp);
}
template <>
double CMath<double>::pow(double base, double exp)
{
    return ::pow(base, exp);
}
template <>
int CMath<int>::pow(int base, int exp)
{
    return powl(base, exp);
}

template <>
float CMath<float>::floor(float a)
{
    return floorf(a);
}
template <>
double CMath<double>::floor(double a)
{
    return ::floor(a);
}
template <>
int CMath<int>::floor(int a)
{
    return a;
}

template <>
float CMath<float>::ceil(float a)
{
    return ceilf(a);
}
template <>
double CMath<double>::ceil(double a)
{
    return ::ceil(a);
}
template <>
int CMath<int>::ceil(int a)
{
    return a;
}

/*
template int CMath<int>::digits2(int a)
{
    if (a > 0x40000000)
        return 0;

    if (a == 0)
        return 0;

    int r = 1;
    int c = 1;

    while (r < a)
    {
        r <<= 1;
        c++;
    }
    return c;
}
*/

template <>
float CMath<float>::sqrt(float a)
{
    return sqrtf(a);
}
template <>
double CMath<double>::sqrt(double a)
{
    return ::sqrt(a);
}
template <>
int CMath<int>::sqrt(int a)
{
    return sqrtl((unsigned long) a);
}

template <>
float CMath<float>::round(float a)
{
    return ::round(a);
}
template <>
double CMath<double>::round(double a)
{
    return roundf(a);
}
template <>
int CMath<int>::round(int a)
{
    return a;
}

template <>
float CMath<float>::aton(const char *s)
{
    return atof(s);
}
template <>
double CMath<double>::aton(const char *s)
{
    return atof(s);
}
template <>
int CMath<int>::aton(const char *s)
{
    return atoi(s);
}

template <>
float CMath<float>::sin(float x)
{
    return sinf(x);
}
template <>
double CMath<double>::sin(double x)
{
    return ::sin(x);
}

template <>
float CMath<float>::cos(float x)
{
    return cosf(x);
}
template <>
double CMath<double>::cos(double x)
{
    return ::cos(x);
}

template <>
float CMath<float>::tan(float x)
{
    return tanf(x);
}
template <>
double CMath<double>::tan(double x)
{
    return ::tan(x);
}

template <>
float CMath<float>::exp(float x)
{
    return ::expf(x);
}
template <>
double CMath<double>::exp(double x)
{
    return ::exp(x);
}

template <class T>
std::ostream& operator<<(std::ostream& os, CMath<T> &v)
{
    return os << v.value;
}

template class CMath<double>;
template class CMath<float>;
template class CMath<int>;

template std::ostream& operator<<(std::ostream &os, CMath<double> &v);
template std::ostream& operator<<(std::ostream &os, CMath<float> &v);
template std::ostream& operator<<(std::ostream &os, CMath<int> &v);
