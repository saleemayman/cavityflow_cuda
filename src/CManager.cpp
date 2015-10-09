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

#include "CManager.hpp"
#include "CSingleton.hpp"

template <class T>
CManager<T>::CManager(CDomain<T> domain, CVector<3, int> subdomainNums) :
        _domain(domain), _lbm_controller(NULL)
{
    this->setSubdomainNums(subdomainNums);
}

template <class T>
CManager<T>::~CManager()
{
    if (_lbm_controller)
        delete _lbm_controller;
}

template <class T>
CDomain<T> CManager<T>::getDomain() const
{
    return _domain;
}

template <class T>
void CManager<T>::setDomain(CDomain<T> grid)
{
    _domain = grid;
}

template <class T>
CVector<3,int> CManager<T>::getSubdomainNums() const
{
    return _subdomain_nums;
}

template <class T>
void CManager<T>::setSubdomainNums(CVector<3,int> subdomainNums)
{
    CVector<3, int> do_size = _domain.getSize();
    if ((do_size[0] % subdomainNums[0] != 0)
            || (do_size[1] % subdomainNums[1] != 0)
            || (do_size[2] % subdomainNums[2] != 0)) {
        throw "Number of subdomains does not match with the grid size!";
    }

    CVector<3, int> tmpSD_size;
    tmpSD_size[0] = do_size[0] / subdomainNums[0];
    tmpSD_size[1] = do_size[1] / subdomainNums[1];
    tmpSD_size[2] = do_size[2] / subdomainNums[2];

    _subdomain_size = tmpSD_size;
    _subdomain_nums = subdomainNums;

    // subdomain lengths
    CVector<3, T> domain_length = _domain.getLength();
    _subdomain_length[0] = domain_length[0] / _subdomain_nums[0];
    _subdomain_length[1] = domain_length[1] / _subdomain_nums[1];
    _subdomain_length[2] = domain_length[2] / _subdomain_nums[2];

#if DEBUG
    std::cout << "NUMBER OF SUBDOMAINS: " << _subdomain_nums << std::endl;
    std::cout << "SUBDOMAIN_SIZE: " << _subdomain_size << std::endl;
    std::cout << "SUBDOMAIN_LENGTHS: " << _subdomain_length << std::endl;
#endif
}

template <class T>
void CManager<T>::initSimulation(int my_rank)
{
    // initialize the boundary condition
	std::vector<Flag> boundaryConditions(6, GHOST_LAYER);
    int id = my_rank;
    if (id < 0)
        id = 0;

    int tmpid = id;
    int nx, ny, nz;
    nx = tmpid % _subdomain_nums[0];
    tmpid /= _subdomain_nums[0];
    ny = tmpid % _subdomain_nums[1];
    tmpid /= _subdomain_nums[1];
    nz = tmpid;
#if DEBUG
    std::cout << "ID: "<< id << " NX: " << nx << " NY: " << ny << " NZ: " << nz << std::endl;
#endif

    // create the subdomains instances for the whole domain
    CVector<3, int> origin(nx * _subdomain_size[0], ny * _subdomain_size[1],
            nz * _subdomain_size[2]);
    CDomain<T> *subdomain = new CDomain<T>(id, _subdomain_size, origin,
            _subdomain_length);

    // Setting the boundary conditions for the current Controller
    if (nx == 0)
    	boundaryConditions[0] = OBSTACLE;
    if (nx == (_subdomain_nums[0] - 1))
    	boundaryConditions[1] = OBSTACLE;
    if (ny == 0)
    	boundaryConditions[2] = OBSTACLE;
    if (ny == (_subdomain_nums[1] - 1))
    	boundaryConditions[3] = OBSTACLE;
    if (nz == 0)
    	boundaryConditions[4] = OBSTACLE;
    if (nz == (_subdomain_nums[2] - 1))
    	boundaryConditions[5] = OBSTACLE;

    _lbm_controller = new CController<T>(id, *subdomain, boundaryConditions, CSingleton<CConfiguration<T> >::getInstance());

    // Initializing the Controller's communication classes based on the already computed boundary conditions
    if (boundaryConditions[0] == GHOST_LAYER) {
        int comm_destination = id - 1;
        CVector<3, int> send_size(1, _subdomain_size[1], _subdomain_size[2]);
        CVector<3, int> recv_size(1, _subdomain_size[1], _subdomain_size[2]);
        CVector<3, int> send_origin(1, 0, 0);
        CVector<3, int> recv_origin(0, 0, 0);
        CVector<3, int> comm_direction(1, 0, 0);
        _lbm_controller->addCommunication(
                new CComm<T>(comm_destination, send_size, recv_size,
                        send_origin, recv_origin, comm_direction));
    }
    if (boundaryConditions[1] == GHOST_LAYER) {
        int comm_destination = id + 1;
        CVector<3, int> send_size(1, _subdomain_size[1], _subdomain_size[2]);
        CVector<3, int> recv_size(1, _subdomain_size[1], _subdomain_size[2]);
        CVector<3, int> send_origin(_subdomain_size[0] - 2, 0, 0);
        CVector<3, int> recv_origin(_subdomain_size[0] - 1, 0, 0);
        CVector<3, int> comm_direction(-1, 0, 0);
        _lbm_controller->addCommunication(
                new CComm<T>(comm_destination, send_size, recv_size,
                        send_origin, recv_origin, comm_direction));
    }
    if (boundaryConditions[2] == GHOST_LAYER) {
        int comm_destination = id - _subdomain_nums[0];
        CVector<3, int> send_size(_subdomain_size[0], 1, _subdomain_size[2]);
        CVector<3, int> recv_size(_subdomain_size[0], 1, _subdomain_size[2]);
        CVector<3, int> send_origin(0, 1, 0);
        CVector<3, int> recv_origin(0, 0, 0);
        CVector<3, int> comm_direction(0, 1, 0);
        _lbm_controller->addCommunication(
                new CComm<T>(comm_destination, send_size, recv_size,
                        send_origin, recv_origin, comm_direction));
    }
    if (boundaryConditions[3] == GHOST_LAYER) {
        int comm_destination = id + _subdomain_nums[0];
        CVector<3, int> send_size(_subdomain_size[0], 1, _subdomain_size[2]);
        CVector<3, int> recv_size(_subdomain_size[0], 1,  _subdomain_size[2]);
        CVector<3, int> send_origin(0, _subdomain_size[1] - 2, 0);
        CVector<3, int> recv_origin(0, _subdomain_size[1] - 1, 0);
        CVector<3, int> comm_direction(0, -1, 0);
        _lbm_controller->addCommunication(
                new CComm<T>(comm_destination, send_size, recv_size,
                        send_origin, recv_origin, comm_direction));
    }
    if (boundaryConditions[4] == GHOST_LAYER) {
        int comm_destination = id - _subdomain_nums[0] * _subdomain_nums[1];
        CVector<3, int> send_size(_subdomain_size[0], _subdomain_size[1], 1);
        CVector<3, int> recv_size(_subdomain_size[0], _subdomain_size[1], 1);
        CVector<3, int> send_origin(0, 0, 1);
        CVector<3, int> recv_origin(0, 0, 0);
        CVector<3, int> comm_direction(0, 0, 1);
        _lbm_controller->addCommunication(
                new CComm<T>(comm_destination, send_size, recv_size,
                        send_origin, recv_origin, comm_direction));
    }
    if (boundaryConditions[5] == GHOST_LAYER) {
        int comm_destination = id + _subdomain_nums[0] * _subdomain_nums[1];
        CVector<3, int> send_size(_subdomain_size[0], _subdomain_size[1], 1);
        CVector<3, int> recv_size(_subdomain_size[0], _subdomain_size[1], 1);
        CVector<3, int> send_origin(0, 0, _subdomain_size[2] - 2);
        CVector<3, int> recv_origin(0, 0, _subdomain_size[2] - 1);
        CVector<3, int> comm_direction(0, 0, -1);
        _lbm_controller->addCommunication(
                new CComm<T>(comm_destination, send_size, recv_size,
                        send_origin, recv_origin, comm_direction));
    }
    if (ny == _subdomain_nums[1] - 1) {
        _lbm_controller->setGeometry();
    }

    /*
     * TODO
     * This piece of code should be obsolete since the ghost layer sizes are now
     * set implicitly by CLbmSolver depending on the domain size.
     */
    //_lbm_controller->addCommToSolver();

}

template <class T>
void CManager<T>::startSimulation()
{
    if (!_lbm_controller)
        throw "CManager: Initialize the simulation before starting it!";
    _lbm_controller->run();

}

template <class T>
CController<T>* CManager<T>::getController() const
{
    return _lbm_controller;
}

template <class T>
void CManager<T>::setController(CController<T>* lbmController)
{
    _lbm_controller = lbmController;
}

template <class T>
CVector<3,int> CManager<T>::getSubdomainSize() const
{
    return _subdomain_size;
}

template class CManager<double>;
template class CManager<float>;
