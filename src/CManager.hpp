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

#ifndef CMANAGER_HPP
#define CMANAGER_HPP

#include "libmath/CVector.hpp"
#include "CDomain.hpp"
#include "CController.hpp"

/*
 * Class CManager is responsible for dividing and assigning the subdomains to different processors.
 */
template<typename T>
class CManager
{
private:
    CDomain<T> _domain; ///< The simulation domain.
    CVector<3,int> _subdomain_size; ///< Each subdomain have the same size which is specified with this class member.
    CVector<3,T> _subdomain_length; ///< Each subdomain have the same lengthes which is specified with this class member.
    CVector<3,int> _subdomain_nums; ///< number of subdomains in each direction.
    CController<T>* _lbm_controller;

public:
    CManager(CDomain<T> domain, CVector<3, int> subdomainNums);
    ~CManager();

    CDomain<T> getDomain() const;
    void setDomain(CDomain<T> grid);
    CVector<3, int> getSubdomainNums() const;
    void setSubdomainNums(CVector<3, int> subdomainNums);
    void initSimulation(int my_rank);
    void startSimulation();
    CController<T>* getController() const;
    void setController(CController<T>* lbmController);
    CVector<3, int> getSubdomainSize() const;
};

#endif
