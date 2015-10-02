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

#include "CSingleton.hpp"

// #include <cassert>

#include "CConfiguration.hpp"

template <class T>
T* CSingleton<T>::getInstance()
{
    if(!instance)
        instance = new T;
    // assert(instance != NULL);

    return instance;
}

template<class T>
T* CSingleton<T>::instance = NULL;

template class CSingleton<CConfiguration<float> >;
template class CSingleton<CConfiguration<double> >;