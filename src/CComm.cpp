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
CComm<T>::CComm(int dstId,
        CVector<3, int> sendSize, CVector<3, int> recvSize,
        CVector<3, int> sendOrigin, CVector<3, int> recvOrigin,
        Direction direction/*,
        Direction sendDirection, Direction recvDirection*/) :
        dstId(dstId),
        sendSize(sendSize), recvSize(recvSize),
        sendOrigin(sendOrigin), recvOrigin(recvOrigin),
        direction(direction)/*,
        sendDirection(sendDirection), recvDirection(recvDirection) */
{
}

template <class T>
CComm<T>::~CComm()
{
}

template <class T>
int CComm<T>::getDstId()
{
    return dstId;
}

template <class T>
void CComm<T>::setDstId(int dstId)
{
    this->dstId = dstId;
}

template <class T>
CVector<3, int> CComm<T>::getSendSize()
{
    return sendSize;
}

template <class T>
void CComm<T>::setSendSize(CVector<3, int> sendSize)
{
    this->sendSize = sendSize;
}

template <class T>
CVector<3, int> CComm<T>::getRecvSize()
{
    return recvSize;
}

template <class T>
void CComm<T>::setRecvSize(CVector<3, int> recvSize)
{
    this->recvSize = recvSize;
}

template <class T>
CVector<3, int> CComm<T>::getSendOrigin()
{
    return sendOrigin;
}

template <class T>
void CComm<T>::setSendOrigin(CVector<3, int> sendOrigin)
{
    this->sendOrigin = sendOrigin;
}

template <class T>
CVector<3, int> CComm<T>::getRecvOrigin()
{
    return recvOrigin;
}

template <class T>
void CComm<T>::setRecvOrigin(CVector<3, int> recvOrigin)
{
    this->recvOrigin = recvOrigin;
}

template <class T>
Direction CComm<T>::getDirection()
{
    return direction;
}

template <class T>
void CComm<T>::setDirection(Direction direction)
{
    this->direction = direction;
}

/*
template <class T>
Direction CComm<T>::getSendDirection()
{
    return sendDirection;
}

template <class T>
void CComm<T>::setSendDirection(Direction sendDirection)
{
    this->sendDirection = sendDirection;
}

template <class T>
Direction CComm<T>::getRecvDirection()
{
    return recvDirection;
}

template <class T>
void CComm<T>::setRecvDirection(Direction recvDirection)
{
    this->recvDirection = recvDirection;
}
*/

template class CComm<double>;
template class CComm<float>;
