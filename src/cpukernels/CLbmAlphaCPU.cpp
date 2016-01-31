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

#include <omp.h>

#include "../libmath/CMath.hpp"
#include "CLbmAlphaCPU.hpp"

template<class T>
CLbmAlphaCPU<T>::CLbmAlphaCPU(
        CVector<3, int> domainSize,
        CVector<3, T> gravitation) :
        domainSize(domainSize),
        gravitation(gravitation),
        numOfCells(domainSize.elements())
{
}

template<class T>
CLbmAlphaCPU<T>::~CLbmAlphaCPU()
{
}

template<class T>
void CLbmAlphaCPU<T>::alphaKernelCPU(
        T* densityDistributions,
        Flag* flags,
        T* densities,
        T* velocities,
        const T inv_tau,
        const T drivenCavityVelocity,
        CVector<3, int> origin,
        CVector<3, int> size,
        CVector<3, int> blockDim,
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

								/*
								 * skip cell if it is a ghost cell. Note: a ghost cell also means
								 * that the cell lies at the boundary of the GPU domain and we
								 * have to skip them as well.
								 */
								if (flags[id] == GHOST_LAYER)
									continue;

								// pointer to density distributions
								//__global T *current_dds = &densityDistributions[id];
								T *current_dds = &densityDistributions[id];

								// +++++++++++
								// +++ DD0 +++
								// +++++++++++
								//
								// 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
								dd0 = *current_dds;     current_dds += numOfCells;
								rho = dd0;
								velocity_x = dd0;
								dd1 = *current_dds;     current_dds += numOfCells;
								rho += dd1;
								velocity_x -= dd1;

								dd2 = *current_dds;     current_dds += numOfCells;
								rho += dd2;
								velocity_y = dd2;
								dd3 = *current_dds;     current_dds += numOfCells;
								rho += dd3;
								velocity_y -= dd3;

								// +++++++++++
								// +++ DD1 +++
								// +++++++++++
								//
								// 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
								dd4 = *current_dds;     current_dds += numOfCells;
								rho += dd4;
								velocity_x += dd4;
								velocity_y += dd4;
								dd5 = *current_dds;     current_dds += numOfCells;
								rho += dd5;
								velocity_x -= dd5;
								velocity_y -= dd5;

								dd6 = *current_dds;     current_dds += numOfCells;
								rho += dd6;
								velocity_x += dd6;
								velocity_y -= dd6;
								dd7 = *current_dds;     current_dds += numOfCells;
								rho += dd7;
								velocity_x -= dd7;
								velocity_y += dd7;

								// +++++++++++
								// +++ DD2 +++
								// +++++++++++
								//
								// 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
								dd8 = *current_dds;     current_dds += numOfCells;
								rho += dd8;
								velocity_x += dd8;
								velocity_z = dd8;
								dd9 = *current_dds;     current_dds += numOfCells;
								rho += dd9;
								velocity_x -= dd9;
								velocity_z -= dd9;

								dd10 = *current_dds;        current_dds += numOfCells;
								rho += dd10;
								velocity_x += dd10;
								velocity_z -= dd10;
								dd11 = *current_dds;        current_dds += numOfCells;
								rho += dd11;
								velocity_x -= dd11;
								velocity_z += dd11;

								// +++++++++++
								// +++ DD3 +++
								// +++++++++++

								// dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
								dd12 = *current_dds;        current_dds += numOfCells;
								rho += dd12;
								velocity_y += dd12;
								velocity_z += dd12;
								dd13 = *current_dds;        current_dds += numOfCells;
								rho += dd13;
								velocity_y -= dd13;
								velocity_z -= dd13;

								dd14 = *current_dds;        current_dds += numOfCells;
								rho += dd14;
								velocity_y += dd14;
								velocity_z -= dd14;
								dd15 = *current_dds;        current_dds += numOfCells;
								rho += dd15;
								velocity_y -= dd15;
								velocity_z += dd15;

								// +++++++++++
								// +++ DD4 +++
								// +++++++++++
								//
								// dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
								dd16 = *current_dds;        current_dds += numOfCells;
								rho += dd16;
								velocity_z += dd16;
								dd17 = *current_dds;        current_dds += numOfCells;
								rho += dd17;
								velocity_z -= dd17;
								dd18 = *current_dds;
								rho += dd18;

							//  rho *= (float)(flag != FLAG_OBSTACLE);


								// to add something to a pointer is faster than subtracting it.
								// thus we restart here with pointer to dd0
								current_dds = &densityDistributions[id];

								/**
								 * instead of storing the density distributions after modification,
								 * we store it during the modifications to hide memory waiting stalls
								 */

								T vel2;     // vel*vel
								T vela2;
								T vela_velb;
							#define tmp rho

								T dd_param; // modified rho as temporary variable
								switch(flags[id])
								{
									case (FLUID):    // this is the whole collision operator
										vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
										dd_param = rho - ((T)3/(T)2)*vel2;

										tmp = gravitation[0]*((T)1/(T)18)*rho;
										vela2 = velocity_x*velocity_x;
										dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
										dd1 -= tmp;
										*current_dds = dd1;     current_dds += numOfCells;

										dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
										dd0 += tmp;
										*current_dds = dd0;     current_dds += numOfCells;

										tmp = gravitation[1]*((T)-1/(T)18)*rho;
										vela2 = velocity_y*velocity_y;
										dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
										dd3 -= tmp;
										*current_dds = dd3;     current_dds += numOfCells;

										dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
										dd2 += tmp;
										*current_dds = dd2;     current_dds += numOfCells;


							#define vela_velb_2 vela2
										/***********************
										 * DD1
										 ***********************/
										vela_velb = velocity_x+velocity_y;
										vela_velb_2 = vela_velb*vela_velb;

										tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;

										dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
										dd5 -= tmp;
										*current_dds = dd5;     current_dds += numOfCells;

										dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
										dd4 += tmp;
										*current_dds = dd4;     current_dds += numOfCells;

										vela_velb = velocity_x-velocity_y;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;
										dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
										dd7 -= tmp;
										*current_dds = dd7;     current_dds += numOfCells;

										dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
										dd6 += tmp;
										*current_dds = dd6;     current_dds += numOfCells;

										/***********************
										 * DD2
										 ***********************/
										vela_velb = velocity_x+velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;

										dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
										dd9 -= tmp;
										*current_dds = dd9;     current_dds += numOfCells;

										dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
										dd8 += tmp;
										*current_dds = dd8;     current_dds += numOfCells;

										tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;
										vela_velb = velocity_x-velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
										dd11 -= tmp;
										*current_dds = dd11;        current_dds += numOfCells;

										dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
										dd10 += tmp;
										*current_dds = dd10;        current_dds += numOfCells;

										/***********************
										 * DD3
										 ***********************/
										vela_velb = velocity_y+velocity_z;
										vela_velb_2 = vela_velb*vela_velb;

										tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;
										dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
										dd13 -= tmp;
										*current_dds = dd13;        current_dds += numOfCells;

										dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
										dd12 += tmp;
										*current_dds = dd12;        current_dds += numOfCells;

										vela_velb = velocity_y-velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;

										dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
										dd15 -= tmp;
										*current_dds = dd15;        current_dds += numOfCells;

										dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
										dd14 += tmp;
										*current_dds = dd14;        current_dds += numOfCells;

							#undef vela_velb_2
										/***********************
										 * DD4
										 ***********************/
										vela2 = velocity_z*velocity_z;

										tmp = gravitation[2]*((T)1/(T)18)*rho;
										dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
										dd17 -= tmp;
										*current_dds = dd17;        current_dds += numOfCells;

										dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
										dd16 += tmp;
										*current_dds = dd16;        current_dds += numOfCells;

										dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
										*current_dds = dd18;

										break;

									case (OBSTACLE): // in case of an obstacle, we bounce back the values

										/**
										 * if we are using only bounce back method, it's not necessary to write back the results.
										 * we even wouldn't need to read the density distribution values.
										 */

										if (storeVelocities)
										{
											velocity_x = (T)0;
											velocity_y = (T)0;
											velocity_z = (T)0;
										}
										break;

									case (VELOCITY_INJECTION):   // this flag specifies the injection area of the fluid
										velocity_x = drivenCavityVelocity;
										velocity_y = (T)0;
										velocity_z = (T)0;

										rho = (T)1;

										vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
										dd_param = rho - ((T)3/(T)2)*(vel2);

										/***********************
										 * DD0
										 ***********************/
										vela2 = velocity_x*velocity_x;
										tmp = gravitation[0]*((T)1/(T)18)*rho;

										dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
										dd1 -= tmp;
										*current_dds = dd1;     current_dds += numOfCells;

										dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
										dd0 += tmp;
										*current_dds = dd0;     current_dds += numOfCells;

										vela2 = velocity_y*velocity_y;
										tmp = gravitation[1]*((T)-1/(T)18)*rho;

										dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
										dd3 -= tmp;
										*current_dds = dd3;     current_dds += numOfCells;

										dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
										dd2 += tmp;
										*current_dds = dd2;     current_dds += numOfCells;


							#define vela_velb_2 vela2
										/***********************
										 * DD1
										 ***********************/
										vela_velb = velocity_x+velocity_y;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;

										dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
										dd5 -= tmp;
										*current_dds = dd5;     current_dds += numOfCells;

										dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
										dd4 += tmp;
										*current_dds = dd4;     current_dds += numOfCells;

										vela_velb = velocity_x-velocity_y;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;

										dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
										dd7 -= tmp;
										*current_dds = dd7;     current_dds += numOfCells;

										dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
										dd6 += tmp;
										*current_dds = dd6;     current_dds += numOfCells;

										/***********************
										 * DD2
										 ***********************/
										vela_velb = velocity_x+velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;

										dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
										dd9 -= tmp;
										*current_dds = dd9;     current_dds += numOfCells;

										dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
										dd8 += tmp;
										*current_dds = dd8;     current_dds += numOfCells;

										vela_velb = velocity_x-velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;

										dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
										dd11 -= tmp;
										*current_dds = dd11;        current_dds += numOfCells;

										dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
										dd10 += tmp;
										*current_dds = dd10;        current_dds += numOfCells;

										/***********************
										 * DD3
										 ***********************/
										vela_velb = velocity_y+velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;

										dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
										dd13 -= tmp;
										*current_dds = dd13;        current_dds += numOfCells;

										dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
										dd12 += tmp;
										*current_dds = dd12;        current_dds += numOfCells;

										vela_velb = velocity_y-velocity_z;
										vela_velb_2 = vela_velb*vela_velb;
										tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;

										dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
										dd15 -= tmp;
										*current_dds = dd15;        current_dds += numOfCells;

										dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
										dd14 += tmp;
										*current_dds = dd14;        current_dds += numOfCells;
							#undef vela_velb_2

										/***********************
										 * DD4
										 ***********************/
										vela2 = velocity_z*velocity_z;

										tmp = gravitation[2]*((T)1/(T)18)*rho;
										dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
										dd17 -= tmp;
										*current_dds = dd17;        current_dds += numOfCells;

										dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
										dd16 += tmp;
										*current_dds = dd16;        current_dds += numOfCells;

										dd18 = eq_dd18(dd_param);
										*current_dds = dd18;
										break;

									case (GHOST_LAYER):
											break;
								}

								if (storeVelocities)
								{
									velocities[id + 0*numOfCells] = velocity_x;
									velocities[id + 1*numOfCells] = velocity_y;
									velocities[id + 2*numOfCells] = velocity_z;
								}

								if (storeDensities)
									densities[id] = rho;
							}
						}
					}
				}
			}
		}
	}
}

template class CLbmAlphaCPU<double>;
template class CLbmAlphaCPU<float>;
