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

#ifndef CCOMM_HPP
#define CCOMM_HPP

#include "libmath/CVector.hpp"

/*
 * Class CComm provides necessary information for communication of data between two subdomains.
 */
template<typename T>
class CComm
{
private:
    int _dstID;
    CVector<3, int> _send_size;
    CVector<3, int> _recv_size;
    CVector<3, int> _send_origin;
    CVector<3, int> _recv_origin;
    CVector<3, int> _comm_direction;

public:
    CComm(int dstID, CVector<3, int> send_size, CVector<3, int> recv_size,
            CVector<3, int> send_origin, CVector<3, int> recv_origin,
            CVector<3, int> comm_direction);
    ~CComm();

    CVector<3, int> getCommDirection() const;
    void setCommDirection(CVector<3, int> normal);
    CVector<3, int> getRecvOrigin() const;
    void setRecvOrigin(CVector<3, int> recvOrigin);
    CVector<3, int> getRecvSize() const;
    void setRecvSize(CVector<3, int> recvSize);
    CVector<3, int> getSendOrigin() const;
    void setSendOrigin(CVector<3, int> sendOrigin);
    CVector<3, int> getSendSize() const;
    void setSendSize(CVector<3, int> sendSize);
    int getDstId() const;
    void setDstId(int dstId);
};

#endif
