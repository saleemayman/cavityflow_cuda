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
        const bool storeDensities,
        const bool storeVelocities)
{
    int id;

    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                id = (origin[2] + k) * (domainSize[0] * domainSize[1]) + (origin[1] + j) * domainSize[0] + (origin[0] + i);

                /*
                 * skip cell if it is a ghost cell. Note: a ghost cell also means
                 * that the cell lies at the boundary of the GPU domain and we
                 * have to skip them as well.
                 */
                if (flags[id] == GHOST_LAYER)
                    continue;

                /*
                 * we have to sum the densities up in a specific order.
                 * otherwise it seems that we run into numerical errors.
                 */

                /* +++++++++++
                 * +++ DD0 +++
                 * +++++++++++
                 * 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
                 */
                dd0 = densityDistributions[id + 0*numOfCells];
                rho = dd0;
                velocity_x = dd0;

                dd1 = densityDistributions[id + 1*numOfCells];
                rho += dd1;
                velocity_x -= dd1;

                dd2 = densityDistributions[id + 2*numOfCells];
                rho += dd2;
                velocity_y = dd2;

                dd3 = densityDistributions[id + 3*numOfCells];
                rho += dd3;
                velocity_y -= dd3;

                /* +++++++++++
                 * +++ DD1 +++
                 * +++++++++++
                 * 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
                 */
                dd4 = densityDistributions[id + 4*numOfCells];
                rho += dd4;
                velocity_x += dd4;
                velocity_y += dd4;

                dd5 = densityDistributions[id + 5*numOfCells];
                rho += dd5;
                velocity_x -= dd5;
                velocity_y -= dd5;

                dd6 = densityDistributions[id + 6*numOfCells];
                rho += dd6;
                velocity_x += dd6;
                velocity_y -= dd6;

                dd7 = densityDistributions[id + 7*numOfCells];
                rho += dd7;
                velocity_x -= dd7;
                velocity_y += dd7;

                /* +++++++++++
                 * +++ DD2 +++
                 * +++++++++++
                 * 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
                 */
                dd8 = densityDistributions[id + 8*numOfCells];
                rho += dd8;
                velocity_x += dd8;
                velocity_z = dd8;

                dd9 = densityDistributions[id + 9*numOfCells];
                rho += dd9;
                velocity_x -= dd9;
                velocity_z -= dd9;

                dd10 = densityDistributions[id + 10*numOfCells];
                rho += dd10;
                velocity_x += dd10;
                velocity_z -= dd10;

                dd11 = densityDistributions[id + 11*numOfCells];
                rho += dd11;
                velocity_x -= dd11;
                velocity_z += dd11;

                /* +++++++++++
                 * +++ DD3 +++
                 * +++++++++++
                 * dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
                 */
                dd12 = densityDistributions[id + 12*numOfCells];
                rho += dd12;
                velocity_y += dd12;
                velocity_z += dd12;

                dd13 = densityDistributions[id + 13*numOfCells];
                rho += dd13;
                velocity_y -= dd13;
                velocity_z -= dd13;

                dd14 = densityDistributions[id + 14*numOfCells];
                rho += dd14;
                velocity_y += dd14;
                velocity_z -= dd14;

                dd15 = densityDistributions[id + 15*numOfCells];
                rho += dd15;
                velocity_y -= dd15;
                velocity_z += dd15;

                /* +++++++++++
                 * +++ DD4 +++
                 * +++++++++++
                 * dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
                 */
                dd16 = densityDistributions[id + 16*numOfCells];
                rho += dd16;
                velocity_z += dd16;

                dd17 = densityDistributions[id + 17*numOfCells];
                rho += dd17;
                velocity_z -= dd17;

                dd18 = densityDistributions[id + 18*numOfCells];
                rho += dd18;


                /**
                 * instead of storing the density distributions after modification,
                 * we store it during the modifications to hide memory waiting stalls
                 */

        #define tmp rho
                T dd_param; // modified rho as temporary variable
                switch(flags[id])
                {
                    case FLUID:    // this is the whole collision operator
                        vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
                        dd_param = rho - ((T)3/(T)2)*(vel2);

                        tmp = gravitation[0]*((T)1/(T)18)*rho;
                        vela2 = velocity_x*velocity_x;
                        dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
                        dd1 -= tmp;
                        densityDistributions[id + 0*numOfCells] = dd1;

                        dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
                        dd0 += tmp;
                        densityDistributions[id + 1*numOfCells] = dd0;

                        tmp = gravitation[1]*((T)-1/(T)18)*rho;
                        vela2 = velocity_y*velocity_y;
                        dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
                        dd3 -= tmp;
                        densityDistributions[id + 2*numOfCells] = dd3;

                        dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
                        dd2 += tmp;
                        densityDistributions[id + 3*numOfCells] = dd2;

        #define vela_velb_2 vela2
                        /***********************
                         * DD1
                         ***********************/
                        vela_velb = velocity_x+velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;

                        tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;

                        dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
                        dd5 -= tmp;
                        densityDistributions[id + 4*numOfCells] = dd5;

                        dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
                        dd4 += tmp;
                        densityDistributions[id + 5*numOfCells] = dd4;

                        vela_velb = velocity_x-velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;
                        dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
                        dd7 -= tmp;
                        densityDistributions[id + 6*numOfCells] = dd7;

                        dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
                        dd6 += tmp;
                        densityDistributions[id + 7*numOfCells] = dd6;

                        /***********************
                         * DD2
                         ***********************/
                        vela_velb = velocity_x+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;

                        dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
                        dd9 -= tmp;
                        densityDistributions[id + 8*numOfCells] = dd9;

                        dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
                        dd8 += tmp;
                        densityDistributions[id + 9*numOfCells] = dd8;

                        tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;
                        vela_velb = velocity_x-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
                        dd11 -= tmp;
                        densityDistributions[id + 10*numOfCells] = dd11;

                        dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
                        dd10 += tmp;
                        densityDistributions[id + 11*numOfCells] = dd10;

                        /***********************
                         * DD3
                         ***********************/
                        vela_velb = velocity_y+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;

                        tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;
                        dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
                        dd13 -= tmp;
                        densityDistributions[id + 12*numOfCells] = dd13;

                        dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
                        dd12 += tmp;
                        densityDistributions[id + 13*numOfCells] = dd12;

                        vela_velb = velocity_y-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;

                        dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
                        dd15 -= tmp;
                        densityDistributions[id + 14*numOfCells] = dd15;

                        dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
                        dd14 += tmp;
                        densityDistributions[id + 15*numOfCells] = dd14;

        #undef vela_velb_2
                        /***********************
                         * DD4
                         ***********************/
                        vela2 = velocity_z*velocity_z;

                        tmp = gravitation[2]*((T)1/(T)18)*rho;
                        dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
                        dd17 -= tmp;
                        densityDistributions[id + 16*numOfCells] = dd17;

                        dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
                        dd16 += tmp;
                        densityDistributions[id + 17*numOfCells] = dd16;

                        dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
                        densityDistributions[id + 18*numOfCells] = dd18;

                        break;

                    case OBSTACLE: // in case of an obstacle, we bounce back the values

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

                    case VELOCITY_INJECTION:   // this flag specifies the injection area of the fluid
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
                        densityDistributions[id + 0*numOfCells] = dd1;

                        dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
                        dd0 += tmp;
                        densityDistributions[id + 1*numOfCells] = dd0;

                        vela2 = velocity_y*velocity_y;
                        tmp = gravitation[1]*((T)-1/(T)18)*rho;

                        dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
                        dd3 -= tmp;
                        densityDistributions[id + 2*numOfCells] = dd3;

                        dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
                        dd2 += tmp;
                        densityDistributions[id + 3*numOfCells] = dd2;

        #define vela_velb_2 vela2
                        /***********************
                         * DD1
                         ***********************/
                        vela_velb = velocity_x+velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] - gravitation[1])*((T)1/(T)36)*rho;

                        dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd5 -= tmp;
                        densityDistributions[id + 4*numOfCells] = dd5;

                        dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd4 += tmp;
                        densityDistributions[id + 5*numOfCells] = dd4;

                        vela_velb = velocity_x-velocity_y;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[1])*((T)1/(T)36)*rho;

                        dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd7 -= tmp;
                        densityDistributions[id + 6*numOfCells] = dd7;

                        dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd6 += tmp;
                        densityDistributions[id + 7*numOfCells] = dd6;

                        /***********************
                         * DD2
                         ***********************/
                        vela_velb = velocity_x+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] + gravitation[2])*((T)1/(T)36)*rho;

                        dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd9 -= tmp;
                        densityDistributions[id + 8*numOfCells] = dd9;

                        dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd8 += tmp;
                        densityDistributions[id + 9*numOfCells] = dd8;

                        vela_velb = velocity_x-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[0] - gravitation[2])*((T)1/(T)36)*rho;

                        dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd11 -= tmp;
                        densityDistributions[id + 10*numOfCells] = dd11;

                        dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd10 += tmp;
                        densityDistributions[id + 11*numOfCells] = dd10;

                        /***********************
                         * DD3
                         ***********************/
                        vela_velb = velocity_y+velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[2] - gravitation[1])*((T)1/(T)36)*rho;

                        dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd13 -= tmp;
                        densityDistributions[id + 12*numOfCells] = dd13;

                        dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd12 += tmp;
                        densityDistributions[id + 13*numOfCells] = dd12;

                        vela_velb = velocity_y-velocity_z;
                        vela_velb_2 = vela_velb*vela_velb;
                        tmp = (gravitation[2] + gravitation[1])*((T)-1/(T)36)*rho;

                        dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
                        dd15 -= tmp;
                        densityDistributions[id + 14*numOfCells] = dd15;

                        dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
                        dd14 += tmp;
                        densityDistributions[id + 15*numOfCells] = dd14;
        #undef vela_velb_2
                        /***********************
                         * DD4
                         ***********************/
                        vela2 = velocity_z*velocity_z;

                        tmp = gravitation[2]*((T)1/(T)18)*rho;
                        dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
                        dd17 -= tmp;
                        densityDistributions[id + 16*numOfCells] = dd17;

                        dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
                        dd16 += tmp;
                        densityDistributions[id + 17*numOfCells] = dd16;

                        dd18 = eq_dd18(dd_param);
                        densityDistributions[id + 18*numOfCells] = dd18;
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

template class CLbmAlphaCPU<double>;
template class CLbmAlphaCPU<float>;
