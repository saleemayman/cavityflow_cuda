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
#include "common.h"

/*
 * Class CComm provides necessary information for communication of data between two subdomains.
 */
template<typename T>
class CComm
{
private:
    int dstId;
    CVector<3, int> sendSize;
    CVector<3, int> recvSize;
    CVector<3, int> sendOrigin;
    CVector<3, int> recvOrigin;
    Direction direction;
    /*
    Direction sendDirection;
    Direction recvDirection;
    */

public:
    CComm(int dstId,
            CVector<3, int> sendSize, CVector<3, int> recvSize,
            CVector<3, int> sendOrigin, CVector<3, int> recvOrigin,
            Direction direction/*,
            Direction sendDirection, Direction recvDirection*/);
    ~CComm();

    int getDstId();
    void setDstId(int dstId);
    CVector<3, int> getSendSize();
    void setSendSize(CVector<3, int> sendSize);
    CVector<3, int> getRecvSize();
    void setRecvSize(CVector<3, int> recvSize);
    CVector<3, int> getSendOrigin();
    void setSendOrigin(CVector<3, int> sendOrigin);
    CVector<3, int> getRecvOrigin();
    void setRecvOrigin(CVector<3, int> recvOrigin);
    Direction getDirection();
    void setDirection(Direction direction);
    /*
    Direction getSendDirection();
    void setSendDirection(Direction sendDirection);
    Direction getRecvDirection();
    void setRecvDirection(Direction recvDirection);
    */
};

#endif
