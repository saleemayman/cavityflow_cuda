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

#include "CComm.hpp"

template <class T>
CComm<T>::CComm(int dstID, CVector<3, int> send_size, CVector<3, int> recv_size,
        CVector<3, int> send_origin, CVector<3, int> recv_origin,
        CVector<3, int> comm_direction) :
        _dstID(dstID), _send_size(send_size), _recv_size(recv_size), _send_origin(
                send_origin), _recv_origin(recv_origin), _comm_direction(
                comm_direction)
{
}

template <class T>
CComm<T>::~CComm()
{
}

template <class T>
CVector<3,int> CComm<T>::getCommDirection() const
{
    return _comm_direction;
}

template <class T>
void CComm<T>::setCommDirection(CVector<3, int> normal)
{
    _comm_direction = normal;
}

template <class T>
CVector<3, int> CComm<T>::getRecvOrigin() const
{
    return _recv_origin;
}

template <class T>
void CComm<T>::setRecvOrigin(CVector<3, int> recvOrigin)
{
    _recv_origin = recvOrigin;
}

template <class T>
CVector<3, int> CComm<T>::getRecvSize() const
{
    return _recv_size;
}

template <class T>
void CComm<T>::setRecvSize(CVector<3, int> recvSize)
{
    _recv_size = recvSize;
}

template <class T>
CVector<3, int> CComm<T>::getSendOrigin() const
{
    return _send_origin;
}

template <class T>
void CComm<T>::setSendOrigin(CVector<3, int> sendOrigin)
{
    _send_origin = sendOrigin;
}

template <class T>
CVector<3, int> CComm<T>::getSendSize() const
{
    return _send_size;
}

template <class T>
void CComm<T>::setSendSize(CVector<3, int> sendSize)
{
    _send_size = sendSize;
}

template <class T>
int CComm<T>::getDstId() const
{
    return _dstID;
}

template <class T>
void CComm<T>::setDstId(int dstId)
{
    _dstID = dstId;
}

template class CComm<double>;
template class CComm<float>;
