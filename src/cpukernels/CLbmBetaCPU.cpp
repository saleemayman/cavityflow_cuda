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

#include "../libmath/CMath.hpp"
#include "CLbmBetaCPU.hpp"

template<class T>
CLbmBetaCPU<T>::CLbmBetaCPU(
        CVector<3, int> domainSize,
        CVector<3, T> gravitation) :
        domainSize(domainSize),
        gravitation(gravitation),
        numOfCells(domainSize.elements())
{
    deltaPosX = 1;
    deltaNegX = numOfCells - 1;
    deltaPosY = domainSize[0];
    deltaNegY = numOfCells - domainSize[0];
    deltaPosZ = domainSize[0] * domainSize[1];
    deltaNegZ = numOfCells - domainSize[0] * domainSize[1];
}

template<class T>
CLbmBetaCPU<T>::~CLbmBetaCPU()
{
}

template<class T>
int CLbmBetaCPU<T>::domainWrap(int A, int numOfCells, bool isPowTwo)
{
    int globalWrappedIndex = (int)isPowTwo*(A & (numOfCells-1)) + (int)(!isPowTwo)*(A % numOfCells);

    return globalWrappedIndex;
}

template<class T>
void CLbmBetaCPU<T>::betaKernelCPU(
        T* densityDistributions,
        Flag* flags,
        T* densities,
        T* velocities,
        const T inv_tau,
        const T drivenCavityVelocity, 
        CVector<3, int> origin,
        CVector<3, int> size,
        CVector<3, int> blockDim,
        const bool isDomainPowOfTwo,
        const bool storeDensities,
        const bool storeVelocities)
{
	for (int blockIdxZ = 0; blockIdxZ < ((size[2] - 1) / blockDim[2]) + 1; blockIdxZ++) {
		for (int blockIdxY = 0; blockIdxY < ((size[1] - 1) / blockDim[1]) + 1; blockIdxY++) {
			for (int blockIdxX = 0; blockIdxX < ((size[0] - 1) / blockDim[0]) + 1; blockIdxX++) {
				// #pragma omp task
				{
					for (int elementIdxZ = 0; elementIdxZ < CMath<int>::min(blockDim[2], size[2] - blockIdxZ * blockDim[2]); elementIdxZ++)
					{
						for (int elementIdxY = 0; elementIdxY < CMath<int>::min(blockDim[1], size[1] - blockIdxY * blockDim[1]); elementIdxY++)
						{
							for (int elementIdxX = 0; elementIdxX < CMath<int>::min(blockDim[0], size[0] - blockIdxX * blockDim[0]); elementIdxX++)
							{
								const int X = blockIdxX * blockDim[0] + elementIdxX;
								const int Y = blockIdxY * blockDim[1] + elementIdxY;
								const int Z = blockIdxZ * blockDim[2] + elementIdxZ;
								const int id = (origin[2] + Z) * (domainSize[0] * domainSize[1]) + (origin[1] + Y) * domainSize[0] + (origin[0] + X);

								T *current_dds = densityDistributions;

								/*
								 * dd 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
								 */
								dd1 = current_dds[domainWrap(id + deltaPosX, numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd0 = current_dds[domainWrap(id + deltaNegX, numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho = dd0;
								velocity_x = dd0;

								/*
								 * we have to sum the densities up in a specific order.
								 * otherwise it seems that we run into numerical errors for fluids with zero velocity.
								 */
								rho += dd1;
								velocity_x -= dd1;

								dd3 = current_dds[domainWrap(id + deltaPosY, numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd2 = current_dds[domainWrap(id + deltaNegY, numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd2;
								velocity_y = dd2;

								rho += dd3;
								velocity_y -= dd3;

								/*
								 * dd 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
								 */
								dd5 = current_dds[domainWrap(id + (deltaPosX + deltaPosY), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd4 = current_dds[domainWrap(id + (deltaNegX + deltaNegY), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd4;
								velocity_x += dd4;
								velocity_y += dd4;

								rho += dd5;
								velocity_x -= dd5;
								velocity_y -= dd5;

								dd7 = current_dds[domainWrap(id + (deltaPosX + deltaNegY), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd6 = current_dds[domainWrap(id + (deltaNegX + deltaPosY), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;


								rho += dd6;
								velocity_x += dd6;
								velocity_y -= dd6;

								rho += dd7;
								velocity_x -= dd7;
								velocity_y += dd7;

								/*
								 * dd 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
								 */
								dd9 = current_dds[domainWrap(id + (deltaPosX + deltaPosZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd8 = current_dds[domainWrap(id + (deltaNegX + deltaNegZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd8;
								velocity_x += dd8;
								velocity_z = dd8;

								rho += dd9;
								velocity_x -= dd9;
								velocity_z -= dd9;

								dd11 = current_dds[domainWrap(id + (deltaPosX + deltaNegZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd10 = current_dds[domainWrap(id + (deltaNegX + deltaPosZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd10;
								velocity_x += dd10;
								velocity_z -= dd10;

								rho += dd11;
								velocity_x -= dd11;
								velocity_z += dd11;

								/*
								 * dd 12-15: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
								 */
								dd13 = current_dds[domainWrap(id + (deltaPosY + deltaPosZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd12 = current_dds[domainWrap(id + (deltaNegY + deltaNegZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd12;
								velocity_y += dd12;
								velocity_z += dd12;

								rho += dd13;
								velocity_y -= dd13;
								velocity_z -= dd13;

								dd15 = current_dds[domainWrap(id + (deltaPosY + deltaNegZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd14 = current_dds[domainWrap(id + (deltaNegY + deltaPosZ), numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd14;
								velocity_y += dd14;
								velocity_z -= dd14;

								rho += dd15;
								velocity_y -= dd15;
								velocity_z += dd15;

								/*
								 * dd 16-18: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
								 */
								dd17 = current_dds[domainWrap(id + deltaPosZ, numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;
								dd16 = current_dds[domainWrap(id + deltaNegZ, numOfCells, isDomainPowOfTwo)];
								current_dds += numOfCells;

								rho += dd16;
								velocity_z += dd16;

								rho += dd17;
								velocity_z -= dd17;

								dd18 = current_dds[id];
								rho += dd18;

								T vel2;     // vel*vel
								T vela2;

							#define vela_velb   vel2
							#define vela_velb_2 vela2
							#define dd_param    rho

								switch(flags[id])
								{
									case (FLUID):    // this is the whole collision operator
										vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
										dd_param = rho - ((T)3/(T)2)*(vel2);

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

									case (OBSTACLE): // in case of an obstacle, we bounce back the values
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

									case (VELOCITY_INJECTION):   // this flag specifies the injection area of the fluid
										velocity_x = drivenCavityVelocity;
										velocity_y = 0;
										velocity_z = 0;

										rho = (T)1;

										vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
										dd_param = rho - ((T)3/(T)2)*vel2;

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
									case (GHOST_LAYER):
										break;
								}

								current_dds = densityDistributions;

								/* f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0) */
								current_dds[domainWrap(id + deltaPosX, numOfCells, isDomainPowOfTwo)] = dd0;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegX, numOfCells, isDomainPowOfTwo)] = dd1;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaPosY, numOfCells, isDomainPowOfTwo)] = dd2;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegY, numOfCells, isDomainPowOfTwo)] = dd3;
								current_dds += numOfCells;

								/* f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0) */
								current_dds[domainWrap(id + deltaPosX + deltaPosY, numOfCells, isDomainPowOfTwo)] = dd4;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegX + deltaNegY, numOfCells, isDomainPowOfTwo)] = dd5;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaPosX + deltaNegY, numOfCells, isDomainPowOfTwo)] = dd6;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegX + deltaPosY, numOfCells, isDomainPowOfTwo)] = dd7;
								current_dds += numOfCells;

								/* f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1) */
								current_dds[domainWrap(id + deltaPosX + deltaPosZ, numOfCells, isDomainPowOfTwo)] = dd8;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegX + deltaNegZ, numOfCells, isDomainPowOfTwo)] = dd9;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaPosX + deltaNegZ, numOfCells, isDomainPowOfTwo)] = dd10;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegX + deltaPosZ, numOfCells, isDomainPowOfTwo)] = dd11;
								current_dds += numOfCells;

								/* f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1) */
								current_dds[domainWrap(id + deltaPosY + deltaPosZ, numOfCells, isDomainPowOfTwo)] = dd12;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegY + deltaNegZ, numOfCells, isDomainPowOfTwo)] = dd13;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaPosY + deltaNegZ, numOfCells, isDomainPowOfTwo)] = dd14;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegY + deltaPosZ, numOfCells, isDomainPowOfTwo)] = dd15;
								current_dds += numOfCells;

								/* f(0,0,1), f(0,0,-1),  f(0,0,0) */
								current_dds[domainWrap(id + deltaPosZ, numOfCells, isDomainPowOfTwo)] = dd16;
								current_dds += numOfCells;
								current_dds[domainWrap(id + deltaNegZ, numOfCells, isDomainPowOfTwo)] = dd17;
								current_dds += numOfCells;
								current_dds[id] = dd18;

								if (storeVelocities)
								{
									velocities[id] = velocity_x;
									velocities[id + 1*numOfCells] = velocity_y;
									velocities[id + 2*numOfCells] = velocity_z;
								}

								if (storeDensities)
								{
									// store density (not necessary)
									densities[id] = rho;
								}
							}
						}
					}
				}
			}
		}
	}
}

template class CLbmBetaCPU<double>;
template class CLbmBetaCPU<float>;
