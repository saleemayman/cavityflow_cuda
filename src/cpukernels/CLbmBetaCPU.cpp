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

#include "CLbmBetaCPU.hpp"

#define GRAVITATION         0

template<class T>
CLbmBetaCPU<T>::CLbmBetaCPU(
                    int domainCellsCPUWithHalo,
                    CVector<3, int> domainSizeWithHalo,
                    //std::vector<int> *localToGlobalIndexMap,
                    CLbmInitCPU<T> *initLbmCPU,
                    CVector<3, T> gravitation) :
                        domainCellsCPUWithHalo(domainCellsCPUWithHalo),
                        domainSizeWithHalo(domainSizeWithHalo),
                        //localToGlobalIndexMap(localToGlobalIndexMap),
                        initLbmCPU(initLbmCPU),
                        gravitation(gravitation)
{
    domainCells = domainSizeWithHalo[0] * domainSizeWithHalo[1] * domainSizeWithHalo[2];
    domainCellsInXYPlane = domainSizeWithHalo[0] * domainSizeWithHalo[1];

    deltaPosX = 1;
    deltaNegX = domainCells - 1;
    deltaPosY = domainSizeWithHalo[0];
    deltaNegY = domainCells - domainSizeWithHalo[0];
    deltaPosZ = domainCellsInXYPlane;
    deltaNegZ = domainCells - domainCellsInXYPlane;
}


template<class T>
CLbmBetaCPU<T>::~CLbmBetaCPU()
{
	if (isSubRegion)
		delete localIndices;
}

template<class T>
int CLbmBetaCPU<T>::domainWrap(int A, int domainCells, bool isPowTwo)
{
    int globalWrappedIndex = (int)isPowTwo*(A & (domainCellsCPUWithHalo-1)) + (int)(!isPowTwo)*(A % domainCellsCPUWithHalo);

#if !TOP_DOWN_DECOMPOSITION
    return initLbmCPU->getLocalIndex(globalWrappedIndex);
#else
	return globalWrappedIndex;
#endif

//    /* 
//     * Map the global index to local (CPU linear id) index by 
//     * searching the map-vector for globalIndex and return the
//     * local index (position of vector element globalWrappedIndex).
//     */
//    std::vector<int>::iterator localIndex;
//    localIndex = std::lower_bound(localToGlobalIndexMap->begin(), 
//                                    localToGlobalIndexMap->end(), 
//                                    globalWrappedIndex);
//    if (localIndex != localToGlobalIndexMap->end() && *localIndex == globalWrappedIndex)
//    {
//        return (int)(localIndex - localToGlobalIndexMap->begin());
//    }
//    else
//    {
//        return -1;
//    }
}


template<class T>
void CLbmBetaCPU<T>::betaKernelCPU(
        std::vector<T> &densityDistributions,
        std::vector<Flag> &flags,
        std::vector<T> &velocities,
        std::vector<T> &densities,
        const T inv_tau,
        const T drivenCavityVelocity, 
        CVector<3, int> origin,
        CVector<3, int> size,
        const bool isDomainPowOfTwo,
        const bool storeDensities,
        const bool storeVelocities)
{
    int cellID; // cell local linear id
    int gid;    // cell global linear id
	int startIndex, endIndex;

    /*
     * Check if computation has to be done for a sub-region, if yes then create a local vector
     * containing the local linear indices of all the linear global cells for which we need to 
     * do the computation (the global cells are specified in size parameter).
     */
    if (size[0] == domainSizeWithHalo[0] && size[1] == domainSizeWithHalo[1] && size[2] == domainSizeWithHalo[2])
    {
        startIndex = 0; 
        endIndex = domainCellsCPUWithHalo;
        isSubRegion = (bool)0;
    }
    else
        isSubRegion = (bool)1;
    
    if (isSubRegion)
    {
        startIndex = 0; 
        endIndex = size[0] * size[1] * size[2]; 
        isSubRegion = (bool)1;
        localIndices = new std::vector<int>(endIndex);

        for (int i = 0; i < size[2]; i++)
        {
            for (int j = 0; j < size[1]; j++)
            {
                for (int k = 0; k < size[0]; k++)
                {
#if !TOP_DOWN_DECOMPOSITION
                    localIndices->operator[](k + j * size[0] + i * size[0] * size[1]) = initLbmCPU->getLocalIndex((origin[0] +  k) + (origin[1] + j) * domainSizeWithHalo[0] + (origin[2] + i) * (domainSizeWithHalo[1] * domainSizeWithHalo[2]));
#else
					localIndices->operator[](k + j * size[0] + i * size[0] * size[1]) = (origin[0] +  k) + (origin[1] + j) * domainSizeWithHalo[0] + (origin[2] + i) * (domainSizeWithHalo[1] * domainSizeWithHalo[2]);
#endif
                }
            }
        }
    }

    /*
     * Iterate over all CPU domain cells in the following order:
     * x-cells, y-cells, z-cells.
     */
    for (int id = startIndex; id < endIndex; id++)
    {
        if (isSubRegion)
        {
            cellID = localIndices->operator[](id);
        }
        else
        {
            cellID = id;
        }
#if !TOP_DOWN_DECOMPOSITION
        gid = initLbmCPU->getGlobalIndex(cellID);
#else
		gid = cellID;
#endif

        /*
         * dd 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
         */
        dd1 = densityDistributions[ domainWrap(gid + deltaPosX, domainCells, isDomainPowOfTwo) + 0*domainCellsCPUWithHalo ];
        dd0 = densityDistributions[ domainWrap(gid + deltaNegX, domainCells, isDomainPowOfTwo) + 1*domainCellsCPUWithHalo ];  
    
        rho = dd0;
        velocity_x = dd0;
    
        /*
         * we have to sum the densities up in a specific order.
         * otherwise it seems that we run into numerical errors for fluids with zero velocity.
         */
        rho += dd1;
        velocity_x -= dd1;
    
        dd3 = densityDistributions[ domainWrap(gid + deltaPosY, domainCells, isDomainPowOfTwo) + 2*domainCellsCPUWithHalo ];
        dd2 = densityDistributions[ domainWrap(gid + deltaNegY, domainCells, isDomainPowOfTwo) + 3*domainCellsCPUWithHalo ];      
    
        rho += dd2;
        velocity_y = dd2;
    
        rho += dd3;
        velocity_y -= dd3;
    
        /*
         * dd 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
         */
        dd5 = densityDistributions[ domainWrap(gid + (deltaPosX + deltaPosY), domainCells, isDomainPowOfTwo) + 4*domainCellsCPUWithHalo ];
        dd4 = densityDistributions[ domainWrap(gid + (deltaNegX + deltaNegY), domainCells, isDomainPowOfTwo) + 5*domainCellsCPUWithHalo ];      
    
        rho += dd4;
        velocity_x += dd4;
        velocity_y += dd4;
    
        rho += dd5;
        velocity_x -= dd5;
        velocity_y -= dd5;
    
        dd7 = densityDistributions[ domainWrap(gid + (deltaPosX + deltaNegY), domainCells, isDomainPowOfTwo) + 6*domainCellsCPUWithHalo ];      
        dd6 = densityDistributions[ domainWrap(gid + (deltaNegX + deltaPosY), domainCells, isDomainPowOfTwo) + 7*domainCellsCPUWithHalo ];
    
        rho += dd6;
        velocity_x += dd6;
        velocity_y -= dd6;
    
        rho += dd7;
        velocity_x -= dd7;
        velocity_y += dd7;
    
        /*
         * dd 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
         */
        dd9 = densityDistributions[ domainWrap(gid + (deltaPosX + deltaPosZ), domainCells, isDomainPowOfTwo) + 8*domainCellsCPUWithHalo ];
        dd8 = densityDistributions[ domainWrap(gid + (deltaNegX + deltaNegZ), domainCells, isDomainPowOfTwo) + 9*domainCellsCPUWithHalo ];
    
        rho += dd8;
        velocity_x += dd8;
        velocity_z = dd8;
    
        rho += dd9;
        velocity_x -= dd9;
        velocity_z -= dd9;
    
        dd11 = densityDistributions[ domainWrap(gid + (deltaPosX + deltaNegZ), domainCells, isDomainPowOfTwo) + 10*domainCellsCPUWithHalo ];
        dd10 = densityDistributions[ domainWrap(gid + (deltaNegX + deltaPosZ), domainCells, isDomainPowOfTwo) + 11*domainCellsCPUWithHalo ];
    
        rho += dd10;
        velocity_x += dd10;
        velocity_z -= dd10;
    
        rho += dd11;
        velocity_x -= dd11;
        velocity_z += dd11;
    
        /*
         * dd 12-15: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
         */
        dd13 = densityDistributions[ domainWrap(gid + (deltaPosY + deltaPosZ), domainCells, isDomainPowOfTwo) + 12*domainCellsCPUWithHalo ];
        dd12 = densityDistributions[ domainWrap(gid + (deltaNegY + deltaNegZ), domainCells, isDomainPowOfTwo) + 13*domainCellsCPUWithHalo ];
    
        rho += dd12;
        velocity_y += dd12;
        velocity_z += dd12;
    
        rho += dd13;
        velocity_y -= dd13;
        velocity_z -= dd13;
    
        dd15 = densityDistributions[ domainWrap(gid + (deltaPosY + deltaNegZ), domainCells, isDomainPowOfTwo) + 14*domainCellsCPUWithHalo ];
        dd14 = densityDistributions[ domainWrap(gid + (deltaNegY + deltaPosZ), domainCells, isDomainPowOfTwo) + 15*domainCellsCPUWithHalo ];
    
        rho += dd14;
        velocity_y += dd14;
        velocity_z -= dd14;
    
        rho += dd15;
        velocity_y -= dd15;
        velocity_z += dd15;
    
        /*
         * dd 16-18: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
         */
        dd17 = densityDistributions[ domainWrap(gid + deltaPosZ, domainCells, isDomainPowOfTwo) + 16*domainCellsCPUWithHalo ];
        dd16 = densityDistributions[ domainWrap(gid + deltaNegZ, domainCells, isDomainPowOfTwo) + 17*domainCellsCPUWithHalo ];
    
        rho += dd16;
        velocity_z += dd16;
    
        rho += dd17;
        velocity_z -= dd17;
    
        dd18 = densityDistributions[cellID + 18*domainCellsCPUWithHalo];
        rho += dd18;
    
        vela_velb = vel2;
        vela_velb_2 = vela2;
        dd_param = rho;
    
        switch(flags[cellID])
        {
            case FLUID:    // this is the whole collision operator
                vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
                dd_param = rho - (T)(3.0f/2.0f)*(vel2);
    
                vela2 = velocity_x*velocity_x;
                dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
                dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
    
                vela2 = velocity_y*velocity_y;
                dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
                dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
    
                /***********************
                 * DD1
                 ***********************/
                vela_velb = velocity_x+velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
    
                dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
                dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
    
                vela_velb = velocity_x-velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
    
                dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
                dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
    
                /***********************
                 * DD2
                 ***********************/
                vela_velb = velocity_x+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
                dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
    
                vela_velb = velocity_x-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
                dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
    
                /***********************
                 * DD3
                 ***********************/
                vela_velb = velocity_y+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
                dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
    
                vela_velb = velocity_y-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
                dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
    
                /***********************
                 * DD4
                 ***********************/
                vela2 = velocity_z*velocity_z;
                dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
                dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
    
                dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
                break;
    
            case OBSTACLE: // in case of an obstacle, we bounce back the values
                // set to zero velocity and no fluid density
    
                if (storeVelocities)
                {
                    velocity_x = (T)0;
                    velocity_y = (T)0;
                    velocity_z = (T)0;
                }
    
                // use simple bounce back
                vela2 = dd1;    dd1 = dd0;      dd0 = vela2;
                vela2 = dd3;    dd3 = dd2;      dd2 = vela2;
                vela2 = dd5;    dd5 = dd4;      dd4 = vela2;
                vela2 = dd7;    dd7 = dd6;      dd6 = vela2;
                vela2 = dd9;    dd9 = dd8;      dd8 = vela2;
                vela2 = dd11;   dd11 = dd10;    dd10 = vela2;
                vela2 = dd13;   dd13 = dd12;    dd12 = vela2;
                vela2 = dd15;   dd15 = dd14;    dd14 = vela2;
                vela2 = dd17;   dd17 = dd16;    dd16 = vela2;
    
                break;
    
            case VELOCITY_INJECTION:   // this flag specifies the injection area of the fluid
                rho = 1.0f;
                velocity_x = drivenCavityVelocity;
                velocity_y = 0;
                velocity_z = 0;
    
                vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
                dd_param = rho - (T)(3.0f/2.0f)*(vel2);
    
                /***********************
                 * DD0
                 ***********************/
                vela2 = velocity_x*velocity_x;
                dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
                dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
    
                vela2 = velocity_y*velocity_y;
                dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
                dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
    
                /***********************
                 * DD1
                 ***********************/
                vela_velb = velocity_x+velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
    
                dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    
                vela_velb = velocity_x-velocity_y;
                vela_velb_2 = vela_velb*vela_velb;
                dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    
                /***********************
                 * DD2
                 ***********************/
                vela_velb = velocity_x+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
    
                dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    
                vela_velb = velocity_x-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    
                /***********************
                 * DD3
                 ***********************/
                vela_velb = velocity_y+velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
    
                dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    
                vela_velb = velocity_y-velocity_z;
                vela_velb_2 = vela_velb*vela_velb;
                dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
    
                /***********************
                 * DD4
                 ***********************/
                vela2 = velocity_z*velocity_z;
                dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
                dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
    
                dd18 = eq_dd18(dd_param);
                break;
            case ( GHOST_LAYER):
                break;
        }
    
#if GRAVITATION
        T tmp = vela2;
        if (flags[cellID] != OBSTACLE)
        {
            tmp = gravitation[0]*(T)(1.0f/18.0f)*rho;
            dd0 += tmp;
            dd1 -= tmp;
            tmp = gravitation[1]*(T)(-1.0f/18.0f)*rho;
            dd2 += tmp;
            dd3 -= tmp;
    
            tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;
            dd4 += tmp;
            dd5 -= tmp;
            tmp = (gravitation[0] + gravitation[1])*(T)((T)1/(T)36)*rho;
            dd6 += tmp;
            dd7 -= tmp;
    
            tmp = (gravitation[0] + gravitation[2])*(T)((T)1/(T)36)*rho;
            dd8 += tmp;
            dd9 -= tmp;
            tmp = (gravitation[0] - gravitation[2])*(T)((T)1/(T)36)*rho;
            dd10 += tmp;
            dd11 -= tmp;
    
            tmp = (gravitation[2] - gravitation[1])*(T)((T)1/(T)36)*rho;
            dd12 += tmp;
            dd13 -= tmp;
            tmp = (gravitation[2] + gravitation[1])*(T)((T)-1/(T)36)*rho;
            dd14 += tmp;
            dd15 -= tmp;
    
            tmp = gravitation[2]*(T)(1.0f/18.0f)*rho;
            dd16 += tmp;
            dd17 -= tmp;
        }
#endif
        /* f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0) */
        densityDistributions[ domainWrap(gid + deltaPosX, domainCells, isDomainPowOfTwo) + 0*domainCellsCPUWithHalo ] = dd0;
        densityDistributions[ domainWrap(gid + deltaNegX, domainCells, isDomainPowOfTwo) + 1*domainCellsCPUWithHalo ] = dd1;
        densityDistributions[ domainWrap(gid + deltaPosY, domainCells, isDomainPowOfTwo) + 2*domainCellsCPUWithHalo ] = dd2;
        densityDistributions[ domainWrap(gid + deltaNegY, domainCells, isDomainPowOfTwo) + 3*domainCellsCPUWithHalo ] = dd3;
    
        /* f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0) */
        densityDistributions[ domainWrap(gid + deltaPosX + deltaPosY, domainCells, isDomainPowOfTwo) + 4*domainCellsCPUWithHalo ] = dd4;
        densityDistributions[ domainWrap(gid + deltaNegX + deltaNegY, domainCells, isDomainPowOfTwo) + 5*domainCellsCPUWithHalo ] = dd5;
        densityDistributions[ domainWrap(gid + deltaPosX + deltaNegY, domainCells, isDomainPowOfTwo) + 6*domainCellsCPUWithHalo ] = dd6;
        densityDistributions[ domainWrap(gid + deltaNegX + deltaPosY, domainCells, isDomainPowOfTwo) + 7*domainCellsCPUWithHalo ] = dd7;
    
        /* f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1) */
        densityDistributions[ domainWrap(gid + deltaPosX + deltaPosZ, domainCells, isDomainPowOfTwo) + 8*domainCellsCPUWithHalo ] = dd8;
        densityDistributions[ domainWrap(gid + deltaNegX + deltaNegZ, domainCells, isDomainPowOfTwo) + 9*domainCellsCPUWithHalo ] = dd9;
        densityDistributions[ domainWrap(gid + deltaPosX + deltaNegZ, domainCells, isDomainPowOfTwo) + 10*domainCellsCPUWithHalo ] = dd10;
        densityDistributions[ domainWrap(gid + deltaNegX + deltaPosZ, domainCells, isDomainPowOfTwo) + 11*domainCellsCPUWithHalo ] = dd11;
    
        /* f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1) */
        densityDistributions[ domainWrap(gid + deltaPosY + deltaPosZ, domainCells, isDomainPowOfTwo) + 12*domainCellsCPUWithHalo ] = dd12;
        densityDistributions[ domainWrap(gid + deltaNegY + deltaNegZ, domainCells, isDomainPowOfTwo) + 13*domainCellsCPUWithHalo ] = dd13;
        densityDistributions[ domainWrap(gid + deltaPosY + deltaNegZ, domainCells, isDomainPowOfTwo) + 14*domainCellsCPUWithHalo ] = dd14;
        densityDistributions[ domainWrap(gid + deltaNegY + deltaPosZ, domainCells, isDomainPowOfTwo) + 15*domainCellsCPUWithHalo ] = dd15;
    
        /* f(0,0,1), f(0,0,-1),  f(0,0,0) */
        densityDistributions[ domainWrap(gid + deltaPosZ, domainCells, isDomainPowOfTwo) + 16*domainCellsCPUWithHalo ] = dd16;
        densityDistributions[ domainWrap(gid + deltaNegZ, domainCells, isDomainPowOfTwo) + 17*domainCellsCPUWithHalo ] = dd17;
        densityDistributions[cellID + 18*domainCellsCPUWithHalo] = dd18;
    
        if ( flags[cellID] == GHOST_LAYER)
            continue;
    
        if (storeVelocities)
        {
            velocities[cellID] = velocity_x;
            velocities[cellID + 1*domainCellsCPUWithHalo] = velocity_y;
            velocities[cellID + 2*domainCellsCPUWithHalo] = velocity_z;
        }
    
        if (storeDensities)
        {
            // store density (not necessary)
            densities[cellID] = rho;
        }
    }
}

template class CLbmBetaCPU<double>;
template class CLbmBetaCPU<float>;
