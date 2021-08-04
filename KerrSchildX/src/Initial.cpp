/*
 *  KerrSchildX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KerrSchildX.
 *
 *  KerrSchildX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KerrSchildX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with KerrSchildX.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Initial.cpp
 *  Fill the ADM grid functions with Kerr initial data.
 */

#include "KerrSchildX.hpp"

/**************************
 * Std. lib. includes     *
 * and external libraries *
 **************************/
#include <cmath>

using Arith::vect;
using Loop::dim;
using Loop::GF3D2;
using Loop::GF3D2layout;
using Loop::loop_all;
using Loop::PointDesc;

using namespace std;

extern "C" void KerrSchildX::KerrSchildX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KerrSchildX_Initial;
  DECLARE_CCTK_PARAMETERS;

  // Grid layout -----------------------------------------------
  const vect<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  // Grid functions --------------------------------------------
  const GF3D2<CCTK_REAL> gxx_(layout, gxx);
  const GF3D2<CCTK_REAL> gxy_(layout, gxy);
  const GF3D2<CCTK_REAL> gxz_(layout, gxz);
  const GF3D2<CCTK_REAL> gyy_(layout, gyy);
  const GF3D2<CCTK_REAL> gyz_(layout, gyz);
  const GF3D2<CCTK_REAL> gzz_(layout, gzz);

  const GF3D2<CCTK_REAL> kxx_(layout, kxx);
  const GF3D2<CCTK_REAL> kxy_(layout, kxy);
  const GF3D2<CCTK_REAL> kxz_(layout, kxz);
  const GF3D2<CCTK_REAL> kyy_(layout, kyy);
  const GF3D2<CCTK_REAL> kyz_(layout, kyz);
  const GF3D2<CCTK_REAL> kzz_(layout, kzz);

  const GF3D2<CCTK_REAL> alp_(layout, alp);

  const GF3D2<CCTK_REAL> dtalp_(layout, dtalp);

  const GF3D2<CCTK_REAL> betax_(layout, betax);
  const GF3D2<CCTK_REAL> betay_(layout, betay);
  const GF3D2<CCTK_REAL> betaz_(layout, betaz);

  const GF3D2<CCTK_REAL> dtbetax_(layout, dtbetax);
  const GF3D2<CCTK_REAL> dtbetay_(layout, dtbetay);
  const GF3D2<CCTK_REAL> dtbetaz_(layout, dtbetaz);

  const CCTK_REAL t = cctk_time;

  auto id_lambda = [&](const PointDesc &p) {
    CCTK_REAL alpL = alp_(p.I);
    CCTK_REAL betaxL = betax_(p.I);
    CCTK_REAL betayL = betay_(p.I);
    CCTK_REAL betazL = betaz_(p.I);
    CCTK_REAL dtbetaxL = dtbetax_(p.I);
    CCTK_REAL dtbetayL = dtbetay_(p.I);
    CCTK_REAL dtbetazL = dtbetaz_(p.I);
    CCTK_REAL gxxL = gxx_(p.I);
    CCTK_REAL gxyL = gxy_(p.I);
    CCTK_REAL gxzL = gxz_(p.I);
    CCTK_REAL gyyL = gyy_(p.I);
    CCTK_REAL gyzL = gyz_(p.I);
    CCTK_REAL gzzL = gzz_(p.I);
    const CCTK_REAL xL = p.x;
    const CCTK_REAL yL = p.y;
    const CCTK_REAL zL = p.z;

    /* The following code was adapted from the Kranc-based Einstein Toolkit
     * thorn "EinsteinExact."
     */

    const CCTK_REAL xform1L00 = 1;

    const CCTK_REAL xform1L01 = 0;

    const CCTK_REAL xform1L02 = 0;

    const CCTK_REAL xform1L03 = 0;

    const CCTK_REAL xform1L10 = 0;

    const CCTK_REAL csetemp0 = cos(phi);

    const CCTK_REAL csetemp1 = cos(psi);

    const CCTK_REAL csetemp2 = cos(theta);

    const CCTK_REAL csetemp3 = sin(phi);

    const CCTK_REAL csetemp4 = sin(psi);

    const CCTK_REAL xform1L11 =
        csetemp0 * csetemp1 - csetemp2 * csetemp3 * csetemp4;

    const CCTK_REAL xform1L12 =
        csetemp1 * csetemp3 + csetemp0 * csetemp2 * csetemp4;

    const CCTK_REAL csetemp5 = sin(theta);

    const CCTK_REAL xform1L13 = csetemp4 * csetemp5;

    const CCTK_REAL xform1L20 = 0;

    const CCTK_REAL xform1L21 =
        -(csetemp1 * csetemp2 * csetemp3) - csetemp0 * csetemp4;

    const CCTK_REAL xform1L22 =
        csetemp0 * csetemp1 * csetemp2 - csetemp3 * csetemp4;

    const CCTK_REAL xform1L23 = csetemp1 * csetemp5;

    const CCTK_REAL xform1L30 = 0;

    const CCTK_REAL xform1L31 = csetemp3 * csetemp5;

    const CCTK_REAL xform1L32 = -(csetemp0 * csetemp5);

    const CCTK_REAL xform1L33 = csetemp2;

    const CCTK_REAL csetemp6 = boostx * boostx;

    const CCTK_REAL csetemp7 = boosty * boosty;

    const CCTK_REAL csetemp8 = boostz * boostz;

    const CCTK_REAL csetemp9 = 1.0 / lapsefactor;

    const CCTK_REAL xform2L00 =
        csetemp9 *
        (-1 + boostx * shiftaddx + boosty * shiftaddy + boostz * shiftaddz) *
        (1 - csetemp6 - csetemp7 - csetemp8 +
         pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L01 =
        boostx * csetemp9 *
        (1 - csetemp6 - csetemp7 - csetemp8 +
         pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L02 =
        boosty * csetemp9 *
        (1 - csetemp6 - csetemp7 - csetemp8 +
         pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L03 =
        boostz * csetemp9 *
        (1 - csetemp6 - csetemp7 - csetemp8 +
         pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L10 =
        (-(boostx * (-1 + csetemp6 + csetemp7 + csetemp8 +
                     (-1 + boosty * shiftaddy + boostz * shiftaddz) *
                         pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5))) +
         shiftaddx * (csetemp6 +
                      (-1 + csetemp7 + csetemp8) *
                          (1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)))) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L11 =
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        (-1 + csetemp7 + csetemp8 +
         csetemp6 * pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1));

    const CCTK_REAL xform2L12 =
        -(boostx * boosty *
          pow(-1 + csetemp6 + csetemp7 + csetemp8 -
                  pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5),
              -1));

    const CCTK_REAL xform2L13 =
        -(boostx * boostz *
          pow(-1 + csetemp6 + csetemp7 + csetemp8 -
                  pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5),
              -1));

    const CCTK_REAL xform2L20 =
        (-(boosty * (-1 + csetemp6 + csetemp7 + csetemp8 +
                     (-1 + boostx * shiftaddx + boostz * shiftaddz) *
                         pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5))) +
         shiftaddy * (csetemp7 +
                      (-1 + csetemp6 + csetemp8) *
                          (1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)))) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L21 =
        -(boostx * boosty *
          pow(-1 + csetemp6 + csetemp7 + csetemp8 -
                  pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5),
              -1));

    const CCTK_REAL xform2L22 =
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        (-1 + csetemp6 + csetemp8 +
         csetemp7 * pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1));

    const CCTK_REAL xform2L23 =
        -(boosty * boostz *
          pow(-1 + csetemp6 + csetemp7 + csetemp8 -
                  pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5),
              -1));

    const CCTK_REAL xform2L30 =
        (shiftaddz * (-1 + csetemp6 + csetemp7 + csetemp8 +
                      (-1 + csetemp6 + csetemp7) *
                          pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)) -
         boostz * (-1 + csetemp6 + csetemp7 + csetemp8 +
                   (-1 + boostx * shiftaddx + boosty * shiftaddy) *
                       pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5))) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xform2L31 =
        -(boostx * boostz *
          pow(-1 + csetemp6 + csetemp7 + csetemp8 -
                  pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5),
              -1));

    const CCTK_REAL xform2L32 =
        -(boosty * boostz *
          pow(-1 + csetemp6 + csetemp7 + csetemp8 -
                  pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5),
              -1));

    const CCTK_REAL xform2L33 =
        (-1 + csetemp6 + csetemp7 + csetemp8 +
         (-1 + csetemp6 + csetemp7) *
             pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5)) *
        pow(-1 + csetemp6 + csetemp7 + csetemp8, -1) *
        pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8, 0.5), -1);

    const CCTK_REAL xformL00 = xform1L00 * xform2L00 + xform1L01 * xform2L10 +
                               xform1L02 * xform2L20 + xform1L03 * xform2L30;

    const CCTK_REAL xformL01 = xform1L00 * xform2L01 + xform1L01 * xform2L11 +
                               xform1L02 * xform2L21 + xform1L03 * xform2L31;

    const CCTK_REAL xformL02 = xform1L00 * xform2L02 + xform1L01 * xform2L12 +
                               xform1L02 * xform2L22 + xform1L03 * xform2L32;

    const CCTK_REAL xformL03 = xform1L00 * xform2L03 + xform1L01 * xform2L13 +
                               xform1L02 * xform2L23 + xform1L03 * xform2L33;

    const CCTK_REAL xformL10 = xform1L10 * xform2L00 + xform1L11 * xform2L10 +
                               xform1L12 * xform2L20 + xform1L13 * xform2L30;

    const CCTK_REAL xformL11 = xform1L10 * xform2L01 + xform1L11 * xform2L11 +
                               xform1L12 * xform2L21 + xform1L13 * xform2L31;

    const CCTK_REAL xformL12 = xform1L10 * xform2L02 + xform1L11 * xform2L12 +
                               xform1L12 * xform2L22 + xform1L13 * xform2L32;

    const CCTK_REAL xformL13 = xform1L10 * xform2L03 + xform1L11 * xform2L13 +
                               xform1L12 * xform2L23 + xform1L13 * xform2L33;

    const CCTK_REAL xformL20 = xform1L20 * xform2L00 + xform1L21 * xform2L10 +
                               xform1L22 * xform2L20 + xform1L23 * xform2L30;

    const CCTK_REAL xformL21 = xform1L20 * xform2L01 + xform1L21 * xform2L11 +
                               xform1L22 * xform2L21 + xform1L23 * xform2L31;

    const CCTK_REAL xformL22 = xform1L20 * xform2L02 + xform1L21 * xform2L12 +
                               xform1L22 * xform2L22 + xform1L23 * xform2L32;

    const CCTK_REAL xformL23 = xform1L20 * xform2L03 + xform1L21 * xform2L13 +
                               xform1L22 * xform2L23 + xform1L23 * xform2L33;

    const CCTK_REAL xformL30 = xform1L30 * xform2L00 + xform1L31 * xform2L10 +
                               xform1L32 * xform2L20 + xform1L33 * xform2L30;

    const CCTK_REAL xformL31 = xform1L30 * xform2L01 + xform1L31 * xform2L11 +
                               xform1L32 * xform2L21 + xform1L33 * xform2L31;

    const CCTK_REAL xformL32 = xform1L30 * xform2L02 + xform1L31 * xform2L12 +
                               xform1L32 * xform2L22 + xform1L33 * xform2L32;

    const CCTK_REAL xformL33 = xform1L30 * xform2L03 + xform1L31 * xform2L13 +
                               xform1L32 * xform2L23 + xform1L33 * xform2L33;

    const CCTK_REAL xx0 = t;

    const CCTK_REAL xx1 = xL - positionx;

    const CCTK_REAL xx2 = yL - positiony;

    const CCTK_REAL xx3 = zL - positionz;

    const CCTK_REAL txx1 =
        xformL10 * xx0 + xformL11 * xx1 + xformL12 * xx2 + xformL13 * xx3;

    const CCTK_REAL txx2 =
        xformL20 * xx0 + xformL21 * xx1 + xformL22 * xx2 + xformL23 * xx3;

    const CCTK_REAL txx3 =
        xformL30 * xx0 + xformL31 * xx1 + xformL32 * xx2 + xformL33 * xx3;

    const CCTK_REAL X = txx1;

    const CCTK_REAL Y = txx2;

    const CCTK_REAL Z = txx3;

    const CCTK_REAL csetemp10 = pow(a, 2);

    const CCTK_REAL csetemp11 = pow(X, 2);

    const CCTK_REAL csetemp12 = pow(Y, 2);

    const CCTK_REAL csetemp13 = pow(Z, 2);

    const CCTK_REAL rXYZ =
        0.707106781186547524400844362105 *
        pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13 +
                pow(4 * csetemp10 * csetemp13 +
                        pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
                    0.5),
            0.5);

    const CCTK_REAL csetemp14 = pow(rXYZ, 3);

    const CCTK_REAL csetemp15 = pow(rXYZ, 4);

    const CCTK_REAL tg400 =
        -1 + 2 * csetemp14 * M * pow(csetemp10 * csetemp13 + csetemp15, -1);

    const CCTK_REAL csetemp16 = pow(rXYZ, 2);

    const CCTK_REAL csetemp17 = rXYZ * X;

    const CCTK_REAL csetemp18 = a * Y;

    const CCTK_REAL tg401 = 2 * csetemp14 * (csetemp17 + csetemp18) * M *
                            pow(csetemp10 * csetemp13 + csetemp15, -1) *
                            pow(csetemp10 + csetemp16, -1);

    const CCTK_REAL csetemp19 = a * X;

    const CCTK_REAL csetemp20 = -csetemp19;

    const CCTK_REAL csetemp21 = rXYZ * Y;

    const CCTK_REAL tg402 = 2 * csetemp14 * (csetemp20 + csetemp21) * M *
                            pow(csetemp10 * csetemp13 + csetemp15, -1) *
                            pow(csetemp10 + csetemp16, -1);

    const CCTK_REAL tg403 =
        2 * csetemp16 * M * Z * pow(csetemp10 * csetemp13 + csetemp15, -1);

    const CCTK_REAL tg411 =
        1 + 2 * csetemp14 * M * pow(csetemp10 * csetemp13 + csetemp15, -1) *
                pow(csetemp10 + csetemp16, -2) * pow(csetemp17 + csetemp18, 2);

    const CCTK_REAL tg412 = 2 * csetemp14 * (csetemp17 + csetemp18) *
                            (csetemp20 + csetemp21) * M *
                            pow(csetemp10 * csetemp13 + csetemp15, -1) *
                            pow(csetemp10 + csetemp16, -2);

    const CCTK_REAL tg413 = 2 * csetemp16 * (csetemp17 + csetemp18) * M * Z *
                            pow(csetemp10 * csetemp13 + csetemp15, -1) *
                            pow(csetemp10 + csetemp16, -1);

    const CCTK_REAL tg422 =
        1 + 2 * csetemp14 * M * pow(csetemp10 * csetemp13 + csetemp15, -1) *
                pow(csetemp10 + csetemp16, -2) * pow(csetemp20 + csetemp21, 2);

    const CCTK_REAL tg423 = 2 * csetemp16 * (csetemp20 + csetemp21) * M * Z *
                            pow(csetemp10 * csetemp13 + csetemp15, -1) *
                            pow(csetemp10 + csetemp16, -1);

    const CCTK_REAL tg433 = 1 + 2 * csetemp13 * M * rXYZ *
                                    pow(csetemp10 * csetemp13 + csetemp15, -1);

    const CCTK_REAL tdg4000 = 0;

    const CCTK_REAL csetemp22 = pow(rXYZ, 7);

    const CCTK_REAL tdg4001 =
        2 * (3 * csetemp10 * csetemp13 * csetemp14 - csetemp22) * M * X *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4002 =
        2 * (3 * csetemp10 * csetemp13 * csetemp14 - csetemp22) * M * Y *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL csetemp23 = pow(rXYZ, 6);

    const CCTK_REAL csetemp24 = pow(a, 4);

    const CCTK_REAL tdg4003 =
        2 *
        (csetemp10 *
             (-5 * csetemp15 +
              (-2 * csetemp10 + 2 * (csetemp11 + csetemp12) + 5 * csetemp13) *
                  csetemp16) -
         csetemp23 + 3 * csetemp13 * csetemp24) *
        M * rXYZ * Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4010 = 0;

    const CCTK_REAL csetemp25 = pow(rXYZ, 9);

    const CCTK_REAL csetemp26 = pow(a, 3);

    const CCTK_REAL csetemp27 = pow(a, 5);

    const CCTK_REAL csetemp28 = pow(rXYZ, 5);

    const CCTK_REAL tdg4011 =
        2 * csetemp14 * M *
        ((3 * csetemp10 - 3 * csetemp11 - csetemp12 - csetemp13) * csetemp22 +
         2 * csetemp25 +
         csetemp10 *
             ((3 * csetemp10 + csetemp11 - csetemp12 - csetemp13) * csetemp13 *
                  csetemp14 +
              (csetemp10 - csetemp11 - csetemp12 + csetemp13) * csetemp28) +
         (-3 * a * csetemp23 - csetemp15 * csetemp26 +
          3 * csetemp13 * csetemp27) *
             X * Y +
         csetemp13 * ((csetemp10 + 3 * csetemp11 - csetemp12 - csetemp13) *
                          csetemp24 * rXYZ +
                      csetemp16 * csetemp26 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL csetemp29 = pow(rXYZ, 8);

    const CCTK_REAL tdg4012 =
        2 * csetemp14 * M *
        (((csetemp10 - csetemp11 - 2 * csetemp12 + csetemp13) * csetemp15 +
          (3 * csetemp10 - csetemp11 - csetemp13) * csetemp13 * csetemp16) *
             csetemp26 +
         a * ((3 * csetemp10 - csetemp11 - 4 * csetemp12 - csetemp13) *
                  csetemp23 +
              2 * csetemp29) +
         (-2 * csetemp22 + 4 * csetemp13 * csetemp24 * rXYZ) * X * Y +
         csetemp13 *
             ((csetemp10 - csetemp11 + 2 * csetemp12 - csetemp13) * csetemp27 +
              2 * csetemp10 * csetemp14 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4013 =
        2 * M * rXYZ *
        ((-2 * csetemp22 +
          csetemp10 *
              (2 * (-csetemp10 + csetemp11 + csetemp12 + 2 * csetemp13) *
                   csetemp14 -
               4 * csetemp28) +
          4 * csetemp13 * csetemp24 * rXYZ) *
             X +
         (-3 * a * csetemp23 +
          (-5 * csetemp15 +
           (-2 * csetemp10 + 2 * (csetemp11 + csetemp12) + 3 * csetemp13) *
               csetemp16) *
              csetemp26 +
          3 * csetemp13 * csetemp27) *
             Y) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -1) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4020 = 0;

    const CCTK_REAL tdg4021 =
        -2 * csetemp14 * M *
        (((csetemp10 - 2 * csetemp11 - csetemp12 + csetemp13) * csetemp15 +
          (3 * csetemp10 - csetemp12 - csetemp13) * csetemp13 * csetemp16) *
             csetemp26 +
         a * ((3 * csetemp10 - 4 * csetemp11 - csetemp12 - csetemp13) *
                  csetemp23 +
              2 * csetemp29) +
         (2 * csetemp22 - 4 * csetemp13 * csetemp24 * rXYZ) * X * Y +
         csetemp13 *
             ((csetemp10 + 2 * csetemp11 - csetemp12 - csetemp13) * csetemp27 -
              2 * csetemp10 * csetemp14 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4022 =
        2 * csetemp14 * M *
        ((3 * csetemp10 - csetemp11 - 3 * csetemp12 - csetemp13) * csetemp22 +
         2 * csetemp25 +
         csetemp10 *
             ((3 * csetemp10 - csetemp11 + csetemp12 - csetemp13) * csetemp13 *
                  csetemp14 +
              (csetemp10 - csetemp11 - csetemp12 + csetemp13) * csetemp28) +
         (3 * a * csetemp23 + csetemp15 * csetemp26 -
          3 * csetemp13 * csetemp27) *
             X * Y +
         csetemp13 * ((csetemp10 - csetemp11 + 3 * csetemp12 - csetemp13) *
                          csetemp24 * rXYZ -
                      csetemp16 * csetemp26 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4023 =
        2 * M * rXYZ *
        ((3 * a * csetemp23 +
          (5 * csetemp15 +
           (2 * csetemp10 - 2 * (csetemp11 + csetemp12) - 3 * csetemp13) *
               csetemp16) *
              csetemp26) *
             X +
         (-2 * csetemp22 +
          csetemp10 *
              (2 * (-csetemp10 + csetemp11 + csetemp12 + 2 * csetemp13) *
                   csetemp14 -
               4 * csetemp28)) *
             Y +
         csetemp13 * (-3 * csetemp27 * X + 4 * csetemp24 * rXYZ * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -1) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4030 = 0;

    const CCTK_REAL tdg4031 =
        4 * (csetemp10 * csetemp13 * csetemp16 - csetemp23) * M * X * Z *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4032 =
        4 * (csetemp10 * csetemp13 * csetemp16 - csetemp23) * M * Y * Z *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4033 =
        2 * (csetemp10 * csetemp13 - csetemp15) *
        (2 * csetemp10 * csetemp13 - 2 * csetemp15 +
         (-csetemp10 + csetemp11 + csetemp12 + 3 * csetemp13) * csetemp16) *
        M * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4110 = 0;

    const CCTK_REAL tdg4111 =
        2 * csetemp14 * (csetemp17 + csetemp18) * M *
        ((6 * csetemp10 - 5 * csetemp11 - 2 * (csetemp12 + csetemp13)) *
             csetemp22 +
         4 * csetemp25 +
         csetemp10 *
             (csetemp13 *
                  (6 * csetemp10 - csetemp11 - 2 * (csetemp12 + csetemp13)) *
                  csetemp14 +
              (-csetemp11 - 2 * csetemp12 + 2 * (csetemp10 + csetemp13)) *
                  csetemp28) +
         (-5 * a * csetemp23 - csetemp15 * csetemp26 +
          3 * csetemp13 * csetemp27) *
             X * Y +
         csetemp13 *
             ((2 * csetemp10 + 3 * csetemp11 - 2 * (csetemp12 + csetemp13)) *
                  csetemp24 * rXYZ -
              csetemp16 * csetemp26 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -3) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4112 =
        2 * csetemp14 * (csetemp17 + csetemp18) * M *
        (((-2 * csetemp11 - 3 * csetemp12 + 2 * (csetemp10 + csetemp13)) *
              csetemp15 +
          csetemp13 *
              (6 * csetemp10 - 3 * csetemp12 - 2 * (csetemp11 + csetemp13)) *
              csetemp16) *
             csetemp26 +
         a * ((6 * csetemp10 - 7 * csetemp12 - 2 * (csetemp11 + csetemp13)) *
                  csetemp23 +
              4 * csetemp29) +
         (-3 * csetemp22 + csetemp10 * csetemp28 +
          5 * csetemp13 * csetemp24 * rXYZ) *
             X * Y +
         csetemp13 *
             ((2 * csetemp10 + csetemp12 - 2 * (csetemp11 + csetemp13)) *
                  csetemp27 +
              csetemp10 * csetemp14 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -3) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4113 =
        2 * (csetemp17 + csetemp18) * M * rXYZ *
        ((-3 * csetemp22 +
          csetemp10 *
              ((-2 * csetemp10 + 2 * (csetemp11 + csetemp12) + 3 * csetemp13) *
                   csetemp14 -
               3 * csetemp28) +
          5 * csetemp13 * csetemp24 * rXYZ) *
             X +
         (-5 * a * csetemp23 +
          (-5 * csetemp15 +
           (-2 * csetemp10 + 2 * (csetemp11 + csetemp12) + csetemp13) *
               csetemp16) *
              csetemp26 +
          3 * csetemp13 * csetemp27) *
             Y) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4120 = 0;

    const CCTK_REAL csetemp30 = pow(rXYZ, 10);

    const CCTK_REAL csetemp31 = pow(a, 6);

    const CCTK_REAL tdg4121 =
        2 * csetemp14 * M *
        ((a * (2 * (-3 * csetemp10 + 3 * csetemp11 - csetemp12 + csetemp13) *
                   csetemp22 -
               4 * csetemp25) +
          csetemp26 *
              (2 * csetemp13 *
                   (-3 * csetemp10 + csetemp11 + csetemp12 + csetemp13) *
                   csetemp14 -
               2 * (csetemp10 - csetemp11 - csetemp12 + csetemp13) *
                   csetemp28)) *
             X +
         (csetemp10 *
              (2 * (-csetemp10 + 3 * csetemp11 + csetemp13) * csetemp23 -
               csetemp15 *
                   (csetemp13 * (csetemp12 + csetemp13) -
                    csetemp10 * (csetemp12 + 2 * (csetemp11 + csetemp13)) +
                    csetemp24)) -
          (-csetemp10 + 4 * csetemp11 + csetemp12 + csetemp13) * csetemp29 +
          2 * csetemp30 +
          csetemp13 * (-csetemp10 - 2 * csetemp11 + csetemp12 + csetemp13) *
              csetemp31) *
             Y -
         2 * csetemp13 *
             ((csetemp10 + csetemp11 - 3 * csetemp12 - csetemp13) * csetemp27 *
                  rXYZ * X +
              (csetemp10 - 3 * csetemp11) * csetemp16 * csetemp24 * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -3) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4122 =
        -2 * csetemp14 * M *
        ((csetemp10 * csetemp15 *
              (csetemp13 * (csetemp11 + csetemp13) -
               csetemp10 * (csetemp11 + 2 * (csetemp12 + csetemp13)) +
               csetemp24) +
          (-csetemp10 + csetemp11 + 4 * csetemp12 + csetemp13) * csetemp29 -
          2 * csetemp30 +
          (csetemp10 - csetemp11 + 2 * csetemp12 - csetemp13) * csetemp13 *
              csetemp31) *
             X +
         (-4 * a * csetemp25 -
          2 * ((csetemp10 - csetemp11 - csetemp12 + csetemp13) * csetemp26 *
                   csetemp28 +
               (csetemp10 - 3 * csetemp11 + csetemp12 - csetemp13) * csetemp13 *
                   csetemp27 * rXYZ)) *
             Y +
         2 * ((csetemp10 * (csetemp10 - 3 * csetemp12 - csetemp13) * csetemp23 +
               (csetemp10 - 3 * csetemp12) * csetemp13 * csetemp16 *
                   csetemp24) *
                  X +
              (a * (-3 * csetemp10 - csetemp11 + 3 * csetemp12 + csetemp13) *
                   csetemp22 +
               csetemp13 *
                   (-3 * csetemp10 + csetemp11 + csetemp12 + csetemp13) *
                   csetemp14 * csetemp26) *
                  Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -3) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4123 =
        2 * M * rXYZ *
        ((csetemp11 - csetemp12) *
             (4 * a * csetemp22 +
              csetemp26 *
                  (-2 * (-csetemp10 + csetemp11 + csetemp12 + csetemp13) *
                       csetemp14 +
                   4 * csetemp28)) +
         (csetemp10 *
              ((2 * (csetemp11 + csetemp12) + 3 * (csetemp10 + csetemp13)) *
                   csetemp15 +
               2 * csetemp23) +
          2 * (csetemp10 - csetemp11 - csetemp12 + 2 * csetemp13) * csetemp16 *
              csetemp24 -
          3 * csetemp29) *
             X * Y +
         csetemp13 * (4 * (-csetemp11 + csetemp12) * csetemp27 * rXYZ -
                      3 * csetemp31 * X * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4130 = 0;

    const CCTK_REAL tdg4131 =
        -2 * csetemp16 * M *
        ((-3 * csetemp10 + 4 * csetemp11 + csetemp12 + csetemp13) * csetemp22 +
         csetemp10 *
             (csetemp13 * (-3 * csetemp10 + csetemp12 + csetemp13) * csetemp14 -
              (csetemp10 - 2 * csetemp11 - csetemp12 + csetemp13) * csetemp28) +
         csetemp13 * (-csetemp10 - 2 * csetemp11 + csetemp12 + csetemp13) *
             csetemp24 * rXYZ +
         (4 * a * csetemp23 + 2 * csetemp15 * csetemp26) * X * Y -
         2 * (csetemp25 + csetemp13 * csetemp27 * X * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4132 =
        -2 * csetemp16 * M *
        ((-((csetemp10 - csetemp11 - 3 * csetemp12 + csetemp13) * csetemp15) +
          csetemp13 * (-3 * csetemp10 + csetemp11 + csetemp12 + csetemp13) *
              csetemp16) *
             csetemp26 +
         a * ((-3 * csetemp10 + csetemp11 + 5 * csetemp12 + csetemp13) *
                  csetemp23 -
              2 * csetemp29) +
         (3 * csetemp22 + csetemp10 * csetemp28 -
          3 * csetemp13 * csetemp24 * rXYZ) *
             X * Y +
         csetemp13 *
             ((-csetemp10 + csetemp11 - csetemp12 + csetemp13) * csetemp27 -
              csetemp10 * csetemp14 * X * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL csetemp32 = pow(Z, 4);

    const CCTK_REAL tdg4133 =
        2 * M *
        ((-((-csetemp10 + csetemp11 + csetemp12 + 4 * csetemp13) * csetemp22) +
          2 * csetemp25 +
          csetemp10 * csetemp13 *
              ((-csetemp10 + csetemp11 + csetemp12 + 2 * csetemp13) *
                   csetemp14 -
               3 * csetemp28) +
          3 * csetemp24 * csetemp32 * rXYZ) *
             X +
         (a * (csetemp10 - csetemp11 - csetemp12 - 5 * csetemp13) * csetemp23 +
          csetemp13 *
              (-4 * csetemp15 +
               (-csetemp10 + csetemp11 + csetemp12 + csetemp13) * csetemp16) *
              csetemp26 +
          2 * (a * csetemp29 + csetemp27 * csetemp32)) *
             Y) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -1) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4220 = 0;

    const CCTK_REAL csetemp33 = -csetemp21;

    const CCTK_REAL tdg4221 =
        2 * csetemp14 * (csetemp19 + csetemp33) * M *
        (((-3 * csetemp11 - 2 * csetemp12 + 2 * (csetemp10 + csetemp13)) *
              csetemp15 +
          csetemp13 *
              (6 * csetemp10 - 3 * csetemp11 - 2 * (csetemp12 + csetemp13)) *
              csetemp16) *
             csetemp26 +
         a * ((6 * csetemp10 - 7 * csetemp11 - 2 * (csetemp12 + csetemp13)) *
                  csetemp23 +
              4 * csetemp29) +
         (3 * csetemp22 - csetemp10 * csetemp28 -
          5 * csetemp13 * csetemp24 * rXYZ) *
             X * Y +
         csetemp13 *
             ((2 * csetemp10 + csetemp11 - 2 * (csetemp12 + csetemp13)) *
                  csetemp27 -
              csetemp10 * csetemp14 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -3) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4222 =
        2 * csetemp14 * (csetemp20 + csetemp21) * M *
        ((6 * csetemp10 - 5 * csetemp12 - 2 * (csetemp11 + csetemp13)) *
             csetemp22 +
         4 * csetemp25 +
         csetemp10 *
             (csetemp13 *
                  (6 * csetemp10 - csetemp12 - 2 * (csetemp11 + csetemp13)) *
                  csetemp14 +
              (-2 * csetemp11 - csetemp12 + 2 * (csetemp10 + csetemp13)) *
                  csetemp28) +
         (5 * a * csetemp23 + csetemp15 * csetemp26 -
          3 * csetemp13 * csetemp27) *
             X * Y +
         csetemp13 *
             ((2 * csetemp10 + 3 * csetemp12 - 2 * (csetemp11 + csetemp13)) *
                  csetemp24 * rXYZ +
              csetemp16 * csetemp26 * X * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -3) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4223 =
        2 * (csetemp19 + csetemp33) * M * rXYZ *
        ((-5 * csetemp15 +
          (-2 * csetemp10 + 2 * (csetemp11 + csetemp12) + csetemp13) *
              csetemp16) *
             csetemp26 * X +
         csetemp10 *
             ((2 * csetemp10 - 2 * (csetemp11 + csetemp12) - 3 * csetemp13) *
                  csetemp14 +
              3 * csetemp28) *
             Y +
         3 * (csetemp13 * csetemp27 * X + csetemp22 * Y) -
         5 * (a * csetemp23 * X + csetemp13 * csetemp24 * rXYZ * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4230 = 0;

    const CCTK_REAL tdg4231 =
        2 * csetemp16 * M *
        ((-((csetemp10 - 3 * csetemp11 - csetemp12 + csetemp13) * csetemp15) +
          csetemp13 * (-3 * csetemp10 + csetemp11 + csetemp12 + csetemp13) *
              csetemp16) *
             csetemp26 +
         a * ((-3 * csetemp10 + 5 * csetemp11 + csetemp12 + csetemp13) *
                  csetemp23 -
              2 * csetemp29) +
         (-3 * csetemp22 - csetemp10 * csetemp28 +
          3 * csetemp13 * csetemp24 * rXYZ) *
             X * Y +
         csetemp13 *
             ((-csetemp10 - csetemp11 + csetemp12 + csetemp13) * csetemp27 +
              csetemp10 * csetemp14 * X * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4232 =
        -2 * csetemp16 * M *
        ((-3 * csetemp10 + csetemp11 + 4 * csetemp12 + csetemp13) * csetemp22 +
         csetemp10 *
             (csetemp13 * (-3 * csetemp10 + csetemp11 + csetemp13) * csetemp14 -
              (csetemp10 - csetemp11 - 2 * csetemp12 + csetemp13) * csetemp28) -
         4 * a * csetemp23 * X * Y -
         2 * (csetemp25 + csetemp15 * csetemp26 * X * Y) +
         csetemp13 * ((-csetemp10 + csetemp11 - 2 * csetemp12 + csetemp13) *
                          csetemp24 * rXYZ +
                      2 * csetemp27 * X * Y)) *
        Z * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4233 =
        -2 * M *
        (a *
             ((csetemp10 - csetemp11 - csetemp12 - 5 * csetemp13) * csetemp23 +
              2 * csetemp29) *
             X +
         ((-csetemp10 + csetemp11 + csetemp12 + 4 * csetemp13) * csetemp22 -
          2 * csetemp25 + 3 * csetemp10 * csetemp13 * csetemp28) *
             Y +
         csetemp13 *
             ((-4 * csetemp15 +
               (-csetemp10 + csetemp11 + csetemp12 + csetemp13) * csetemp16) *
                  csetemp26 * X +
              csetemp10 * (csetemp10 - csetemp11 - csetemp12 - 2 * csetemp13) *
                  csetemp14 * Y) +
         csetemp32 * (2 * csetemp27 * X - 3 * csetemp24 * rXYZ * Y)) *
        pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(csetemp10 + csetemp16, -1) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4330 = 0;

    const CCTK_REAL tdg4331 =
        2 * csetemp13 * M * (-3 * csetemp28 + csetemp10 * csetemp13 * rXYZ) *
        X * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4332 =
        2 * csetemp13 * M * (-3 * csetemp28 + csetemp10 * csetemp13 * rXYZ) *
        Y * pow(csetemp10 * csetemp13 + csetemp15, -2) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL tdg4333 =
        2 *
        ((2 * csetemp10 - 2 * (csetemp11 + csetemp12) - 5 * csetemp13) *
             csetemp23 +
         4 * csetemp29 + csetemp24 * csetemp32 +
         csetemp10 * (-3 * csetemp13 * csetemp15 + csetemp16 * csetemp32)) *
        M * Z * pow(csetemp10 * csetemp13 + csetemp15, -2) * pow(rXYZ, -1) *
        pow(4 * csetemp10 * csetemp13 +
                pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13, 2),
            -0.5);

    const CCTK_REAL g400 =
        2 * (tg423 * xformL20 * xformL30 +
             xformL00 *
                 (tg401 * xformL10 + tg402 * xformL20 + tg403 * xformL30) +
             xformL10 * (tg412 * xformL20 + tg413 * xformL30)) +
        tg400 * pow(xformL00, 2) + tg411 * pow(xformL10, 2) +
        tg422 * pow(xformL20, 2) + tg433 * pow(xformL30, 2);

    const CCTK_REAL csetemp34 = tg400 * xformL01;

    const CCTK_REAL csetemp35 = tg401 * xformL11;

    const CCTK_REAL csetemp36 = tg402 * xformL21;

    const CCTK_REAL csetemp37 = tg403 * xformL31;

    const CCTK_REAL csetemp38 = tg401 * xformL01;

    const CCTK_REAL csetemp39 = tg411 * xformL11;

    const CCTK_REAL csetemp40 = tg412 * xformL21;

    const CCTK_REAL csetemp41 = tg413 * xformL31;

    const CCTK_REAL csetemp42 = tg402 * xformL01;

    const CCTK_REAL csetemp43 = tg412 * xformL11;

    const CCTK_REAL csetemp44 = tg422 * xformL21;

    const CCTK_REAL csetemp45 = tg423 * xformL31;

    const CCTK_REAL csetemp46 = tg403 * xformL01;

    const CCTK_REAL csetemp47 = tg413 * xformL11;

    const CCTK_REAL csetemp48 = tg423 * xformL21;

    const CCTK_REAL csetemp49 = tg433 * xformL31;

    const CCTK_REAL g401 =
        (csetemp34 + csetemp35 + csetemp36 + csetemp37) * xformL00 +
        (csetemp38 + csetemp39 + csetemp40 + csetemp41) * xformL10 +
        (csetemp42 + csetemp43 + csetemp44 + csetemp45) * xformL20 +
        (csetemp46 + csetemp47 + csetemp48 + csetemp49) * xformL30;

    const CCTK_REAL csetemp50 = tg400 * xformL02;

    const CCTK_REAL csetemp51 = tg401 * xformL12;

    const CCTK_REAL csetemp52 = tg402 * xformL22;

    const CCTK_REAL csetemp53 = tg403 * xformL32;

    const CCTK_REAL csetemp54 = tg401 * xformL02;

    const CCTK_REAL csetemp55 = tg411 * xformL12;

    const CCTK_REAL csetemp56 = tg412 * xformL22;

    const CCTK_REAL csetemp57 = tg413 * xformL32;

    const CCTK_REAL csetemp58 = tg402 * xformL02;

    const CCTK_REAL csetemp59 = tg412 * xformL12;

    const CCTK_REAL csetemp60 = tg422 * xformL22;

    const CCTK_REAL csetemp61 = tg423 * xformL32;

    const CCTK_REAL csetemp62 = tg403 * xformL02;

    const CCTK_REAL csetemp63 = tg413 * xformL12;

    const CCTK_REAL csetemp64 = tg423 * xformL22;

    const CCTK_REAL csetemp65 = tg433 * xformL32;

    const CCTK_REAL g402 =
        (csetemp50 + csetemp51 + csetemp52 + csetemp53) * xformL00 +
        (csetemp54 + csetemp55 + csetemp56 + csetemp57) * xformL10 +
        (csetemp58 + csetemp59 + csetemp60 + csetemp61) * xformL20 +
        (csetemp62 + csetemp63 + csetemp64 + csetemp65) * xformL30;

    const CCTK_REAL csetemp66 = tg400 * xformL03;

    const CCTK_REAL csetemp67 = tg401 * xformL13;

    const CCTK_REAL csetemp68 = tg402 * xformL23;

    const CCTK_REAL csetemp69 = tg403 * xformL33;

    const CCTK_REAL csetemp70 = tg401 * xformL03;

    const CCTK_REAL csetemp71 = tg411 * xformL13;

    const CCTK_REAL csetemp72 = tg412 * xformL23;

    const CCTK_REAL csetemp73 = tg413 * xformL33;

    const CCTK_REAL csetemp74 = tg402 * xformL03;

    const CCTK_REAL csetemp75 = tg412 * xformL13;

    const CCTK_REAL csetemp76 = tg422 * xformL23;

    const CCTK_REAL csetemp77 = tg423 * xformL33;

    const CCTK_REAL csetemp78 = tg403 * xformL03;

    const CCTK_REAL csetemp79 = tg413 * xformL13;

    const CCTK_REAL csetemp80 = tg423 * xformL23;

    const CCTK_REAL csetemp81 = tg433 * xformL33;

    const CCTK_REAL g403 =
        (csetemp66 + csetemp67 + csetemp68 + csetemp69) * xformL00 +
        (csetemp70 + csetemp71 + csetemp72 + csetemp73) * xformL10 +
        (csetemp74 + csetemp75 + csetemp76 + csetemp77) * xformL20 +
        (csetemp78 + csetemp79 + csetemp80 + csetemp81) * xformL30;

    const CCTK_REAL g411 =
        (csetemp34 + csetemp35 + csetemp36 + csetemp37) * xformL01 +
        (csetemp38 + csetemp39 + csetemp40 + csetemp41) * xformL11 +
        (csetemp42 + csetemp43 + csetemp44 + csetemp45) * xformL21 +
        (csetemp46 + csetemp47 + csetemp48 + csetemp49) * xformL31;

    const CCTK_REAL g412 =
        (csetemp50 + csetemp51 + csetemp52 + csetemp53) * xformL01 +
        (csetemp54 + csetemp55 + csetemp56 + csetemp57) * xformL11 +
        (csetemp58 + csetemp59 + csetemp60 + csetemp61) * xformL21 +
        (csetemp62 + csetemp63 + csetemp64 + csetemp65) * xformL31;

    const CCTK_REAL g413 =
        (csetemp66 + csetemp67 + csetemp68 + csetemp69) * xformL01 +
        (csetemp70 + csetemp71 + csetemp72 + csetemp73) * xformL11 +
        (csetemp74 + csetemp75 + csetemp76 + csetemp77) * xformL21 +
        (csetemp78 + csetemp79 + csetemp80 + csetemp81) * xformL31;

    const CCTK_REAL g422 =
        (csetemp50 + csetemp51 + csetemp52 + csetemp53) * xformL02 +
        (csetemp54 + csetemp55 + csetemp56 + csetemp57) * xformL12 +
        (csetemp58 + csetemp59 + csetemp60 + csetemp61) * xformL22 +
        (csetemp62 + csetemp63 + csetemp64 + csetemp65) * xformL32;

    const CCTK_REAL g423 =
        (csetemp66 + csetemp67 + csetemp68 + csetemp69) * xformL02 +
        (csetemp70 + csetemp71 + csetemp72 + csetemp73) * xformL12 +
        (csetemp74 + csetemp75 + csetemp76 + csetemp77) * xformL22 +
        (csetemp78 + csetemp79 + csetemp80 + csetemp81) * xformL32;

    const CCTK_REAL g433 =
        (csetemp66 + csetemp67 + csetemp68 + csetemp69) * xformL03 +
        (csetemp70 + csetemp71 + csetemp72 + csetemp73) * xformL13 +
        (csetemp74 + csetemp75 + csetemp76 + csetemp77) * xformL23 +
        (csetemp78 + csetemp79 + csetemp80 + csetemp81) * xformL33;

    const CCTK_REAL csetemp82 = tdg4000 * xformL00;

    const CCTK_REAL csetemp83 = tdg4001 * xformL10;

    const CCTK_REAL csetemp84 = tdg4002 * xformL20;

    const CCTK_REAL csetemp85 = tdg4003 * xformL30;

    const CCTK_REAL csetemp86 = tdg4010 * xformL00;

    const CCTK_REAL csetemp87 = tdg4011 * xformL10;

    const CCTK_REAL csetemp88 = tdg4012 * xformL20;

    const CCTK_REAL csetemp89 = tdg4013 * xformL30;

    const CCTK_REAL csetemp90 = tdg4020 * xformL00;

    const CCTK_REAL csetemp91 = tdg4021 * xformL10;

    const CCTK_REAL csetemp92 = tdg4022 * xformL20;

    const CCTK_REAL csetemp93 = tdg4023 * xformL30;

    const CCTK_REAL csetemp94 = tdg4030 * xformL00;

    const CCTK_REAL csetemp95 = tdg4031 * xformL10;

    const CCTK_REAL csetemp96 = tdg4032 * xformL20;

    const CCTK_REAL csetemp97 = tdg4033 * xformL30;

    const CCTK_REAL csetemp98 = tdg4110 * xformL00;

    const CCTK_REAL csetemp99 = tdg4111 * xformL10;

    const CCTK_REAL csetemp100 = tdg4112 * xformL20;

    const CCTK_REAL csetemp101 = tdg4113 * xformL30;

    const CCTK_REAL csetemp102 = tdg4120 * xformL00;

    const CCTK_REAL csetemp103 = tdg4121 * xformL10;

    const CCTK_REAL csetemp104 = tdg4122 * xformL20;

    const CCTK_REAL csetemp105 = tdg4123 * xformL30;

    const CCTK_REAL csetemp106 = tdg4130 * xformL00;

    const CCTK_REAL csetemp107 = tdg4131 * xformL10;

    const CCTK_REAL csetemp108 = tdg4132 * xformL20;

    const CCTK_REAL csetemp109 = tdg4133 * xformL30;

    const CCTK_REAL csetemp110 = tdg4220 * xformL00;

    const CCTK_REAL csetemp111 = tdg4221 * xformL10;

    const CCTK_REAL csetemp112 = tdg4222 * xformL20;

    const CCTK_REAL csetemp113 = tdg4223 * xformL30;

    const CCTK_REAL csetemp114 = tdg4230 * xformL00;

    const CCTK_REAL csetemp115 = tdg4231 * xformL10;

    const CCTK_REAL csetemp116 = tdg4232 * xformL20;

    const CCTK_REAL csetemp117 = tdg4233 * xformL30;

    const CCTK_REAL csetemp118 = tdg4330 * xformL00;

    const CCTK_REAL csetemp119 = tdg4331 * xformL10;

    const CCTK_REAL csetemp120 = tdg4332 * xformL20;

    const CCTK_REAL csetemp121 = tdg4333 * xformL30;

    const CCTK_REAL dg4000 =
        2 * ((csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL20 *
                 xformL30 +
             xformL10 * ((csetemp102 + csetemp103 + csetemp104 + csetemp105) *
                             xformL20 +
                         (csetemp106 + csetemp107 + csetemp108 + csetemp109) *
                             xformL30) +
             xformL00 *
                 ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL10 +
                  (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL20 +
                  (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL30)) +
        (csetemp82 + csetemp83 + csetemp84 + csetemp85) * pow(xformL00, 2) +
        (csetemp100 + csetemp101 + csetemp98 + csetemp99) * pow(xformL10, 2) +
        (csetemp110 + csetemp111 + csetemp112 + csetemp113) * pow(xformL20, 2) +
        (csetemp118 + csetemp119 + csetemp120 + csetemp121) * pow(xformL30, 2);

    const CCTK_REAL dg4010 =
        xformL10 *
            ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL01 +
             (csetemp100 + csetemp101 + csetemp98 + csetemp99) * xformL11 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL21 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL31) +
        xformL20 *
            ((csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL01 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL11 +
             (csetemp110 + csetemp111 + csetemp112 + csetemp113) * xformL21 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL31) +
        xformL30 *
            ((csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL01 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL11 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL21 +
             (csetemp118 + csetemp119 + csetemp120 + csetemp121) * xformL31) +
        xformL00 * ((csetemp82 + csetemp83 + csetemp84 + csetemp85) * xformL01 +
                    (csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL11 +
                    (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL21 +
                    (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL31);

    const CCTK_REAL csetemp122 = tdg4000 * xformL01;

    const CCTK_REAL csetemp123 = tdg4001 * xformL11;

    const CCTK_REAL csetemp124 = tdg4002 * xformL21;

    const CCTK_REAL csetemp125 = tdg4003 * xformL31;

    const CCTK_REAL csetemp126 = tdg4010 * xformL01;

    const CCTK_REAL csetemp127 = tdg4011 * xformL11;

    const CCTK_REAL csetemp128 = tdg4012 * xformL21;

    const CCTK_REAL csetemp129 = tdg4013 * xformL31;

    const CCTK_REAL csetemp130 = tdg4020 * xformL01;

    const CCTK_REAL csetemp131 = tdg4021 * xformL11;

    const CCTK_REAL csetemp132 = tdg4022 * xformL21;

    const CCTK_REAL csetemp133 = tdg4023 * xformL31;

    const CCTK_REAL csetemp134 = tdg4030 * xformL01;

    const CCTK_REAL csetemp135 = tdg4031 * xformL11;

    const CCTK_REAL csetemp136 = tdg4032 * xformL21;

    const CCTK_REAL csetemp137 = tdg4033 * xformL31;

    const CCTK_REAL csetemp138 = tdg4110 * xformL01;

    const CCTK_REAL csetemp139 = tdg4111 * xformL11;

    const CCTK_REAL csetemp140 = tdg4112 * xformL21;

    const CCTK_REAL csetemp141 = tdg4113 * xformL31;

    const CCTK_REAL csetemp142 = tdg4120 * xformL01;

    const CCTK_REAL csetemp143 = tdg4121 * xformL11;

    const CCTK_REAL csetemp144 = tdg4122 * xformL21;

    const CCTK_REAL csetemp145 = tdg4123 * xformL31;

    const CCTK_REAL csetemp146 = tdg4130 * xformL01;

    const CCTK_REAL csetemp147 = tdg4131 * xformL11;

    const CCTK_REAL csetemp148 = tdg4132 * xformL21;

    const CCTK_REAL csetemp149 = tdg4133 * xformL31;

    const CCTK_REAL csetemp150 = tdg4220 * xformL01;

    const CCTK_REAL csetemp151 = tdg4221 * xformL11;

    const CCTK_REAL csetemp152 = tdg4222 * xformL21;

    const CCTK_REAL csetemp153 = tdg4223 * xformL31;

    const CCTK_REAL csetemp154 = tdg4230 * xformL01;

    const CCTK_REAL csetemp155 = tdg4231 * xformL11;

    const CCTK_REAL csetemp156 = tdg4232 * xformL21;

    const CCTK_REAL csetemp157 = tdg4233 * xformL31;

    const CCTK_REAL csetemp158 = tdg4330 * xformL01;

    const CCTK_REAL csetemp159 = tdg4331 * xformL11;

    const CCTK_REAL csetemp160 = tdg4332 * xformL21;

    const CCTK_REAL csetemp161 = tdg4333 * xformL31;

    const CCTK_REAL dg4011 =
        xformL00 *
            ((csetemp122 + csetemp123 + csetemp124 + csetemp125) * xformL01 +
             (csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL11 +
             (csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL21 +
             (csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL31) +
        xformL10 *
            ((csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL01 +
             (csetemp138 + csetemp139 + csetemp140 + csetemp141) * xformL11 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL21 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL31) +
        xformL20 *
            ((csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL01 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL11 +
             (csetemp150 + csetemp151 + csetemp152 + csetemp153) * xformL21 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL31) +
        xformL30 *
            ((csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL01 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL11 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL21 +
             (csetemp158 + csetemp159 + csetemp160 + csetemp161) * xformL31);

    const CCTK_REAL csetemp162 = tdg4000 * xformL02;

    const CCTK_REAL csetemp163 = tdg4001 * xformL12;

    const CCTK_REAL csetemp164 = tdg4002 * xformL22;

    const CCTK_REAL csetemp165 = tdg4003 * xformL32;

    const CCTK_REAL csetemp166 = tdg4010 * xformL02;

    const CCTK_REAL csetemp167 = tdg4011 * xformL12;

    const CCTK_REAL csetemp168 = tdg4012 * xformL22;

    const CCTK_REAL csetemp169 = tdg4013 * xformL32;

    const CCTK_REAL csetemp170 = tdg4020 * xformL02;

    const CCTK_REAL csetemp171 = tdg4021 * xformL12;

    const CCTK_REAL csetemp172 = tdg4022 * xformL22;

    const CCTK_REAL csetemp173 = tdg4023 * xformL32;

    const CCTK_REAL csetemp174 = tdg4030 * xformL02;

    const CCTK_REAL csetemp175 = tdg4031 * xformL12;

    const CCTK_REAL csetemp176 = tdg4032 * xformL22;

    const CCTK_REAL csetemp177 = tdg4033 * xformL32;

    const CCTK_REAL csetemp178 = tdg4110 * xformL02;

    const CCTK_REAL csetemp179 = tdg4111 * xformL12;

    const CCTK_REAL csetemp180 = tdg4112 * xformL22;

    const CCTK_REAL csetemp181 = tdg4113 * xformL32;

    const CCTK_REAL csetemp182 = tdg4120 * xformL02;

    const CCTK_REAL csetemp183 = tdg4121 * xformL12;

    const CCTK_REAL csetemp184 = tdg4122 * xformL22;

    const CCTK_REAL csetemp185 = tdg4123 * xformL32;

    const CCTK_REAL csetemp186 = tdg4130 * xformL02;

    const CCTK_REAL csetemp187 = tdg4131 * xformL12;

    const CCTK_REAL csetemp188 = tdg4132 * xformL22;

    const CCTK_REAL csetemp189 = tdg4133 * xformL32;

    const CCTK_REAL csetemp190 = tdg4220 * xformL02;

    const CCTK_REAL csetemp191 = tdg4221 * xformL12;

    const CCTK_REAL csetemp192 = tdg4222 * xformL22;

    const CCTK_REAL csetemp193 = tdg4223 * xformL32;

    const CCTK_REAL csetemp194 = tdg4230 * xformL02;

    const CCTK_REAL csetemp195 = tdg4231 * xformL12;

    const CCTK_REAL csetemp196 = tdg4232 * xformL22;

    const CCTK_REAL csetemp197 = tdg4233 * xformL32;

    const CCTK_REAL csetemp198 = tdg4330 * xformL02;

    const CCTK_REAL csetemp199 = tdg4331 * xformL12;

    const CCTK_REAL csetemp200 = tdg4332 * xformL22;

    const CCTK_REAL csetemp201 = tdg4333 * xformL32;

    const CCTK_REAL dg4012 =
        xformL00 *
            ((csetemp162 + csetemp163 + csetemp164 + csetemp165) * xformL01 +
             (csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL11 +
             (csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL21 +
             (csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL31) +
        xformL10 *
            ((csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL01 +
             (csetemp178 + csetemp179 + csetemp180 + csetemp181) * xformL11 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL21 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL31) +
        xformL20 *
            ((csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL01 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL11 +
             (csetemp190 + csetemp191 + csetemp192 + csetemp193) * xformL21 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL31) +
        xformL30 *
            ((csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL01 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL11 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL21 +
             (csetemp198 + csetemp199 + csetemp200 + csetemp201) * xformL31);

    const CCTK_REAL csetemp202 = tdg4000 * xformL03;

    const CCTK_REAL csetemp203 = tdg4001 * xformL13;

    const CCTK_REAL csetemp204 = tdg4002 * xformL23;

    const CCTK_REAL csetemp205 = tdg4003 * xformL33;

    const CCTK_REAL csetemp206 = tdg4010 * xformL03;

    const CCTK_REAL csetemp207 = tdg4011 * xformL13;

    const CCTK_REAL csetemp208 = tdg4012 * xformL23;

    const CCTK_REAL csetemp209 = tdg4013 * xformL33;

    const CCTK_REAL csetemp210 = tdg4020 * xformL03;

    const CCTK_REAL csetemp211 = tdg4021 * xformL13;

    const CCTK_REAL csetemp212 = tdg4022 * xformL23;

    const CCTK_REAL csetemp213 = tdg4023 * xformL33;

    const CCTK_REAL csetemp214 = tdg4030 * xformL03;

    const CCTK_REAL csetemp215 = tdg4031 * xformL13;

    const CCTK_REAL csetemp216 = tdg4032 * xformL23;

    const CCTK_REAL csetemp217 = tdg4033 * xformL33;

    const CCTK_REAL csetemp218 = tdg4110 * xformL03;

    const CCTK_REAL csetemp219 = tdg4111 * xformL13;

    const CCTK_REAL csetemp220 = tdg4112 * xformL23;

    const CCTK_REAL csetemp221 = tdg4113 * xformL33;

    const CCTK_REAL csetemp222 = tdg4120 * xformL03;

    const CCTK_REAL csetemp223 = tdg4121 * xformL13;

    const CCTK_REAL csetemp224 = tdg4122 * xformL23;

    const CCTK_REAL csetemp225 = tdg4123 * xformL33;

    const CCTK_REAL csetemp226 = tdg4130 * xformL03;

    const CCTK_REAL csetemp227 = tdg4131 * xformL13;

    const CCTK_REAL csetemp228 = tdg4132 * xformL23;

    const CCTK_REAL csetemp229 = tdg4133 * xformL33;

    const CCTK_REAL csetemp230 = tdg4220 * xformL03;

    const CCTK_REAL csetemp231 = tdg4221 * xformL13;

    const CCTK_REAL csetemp232 = tdg4222 * xformL23;

    const CCTK_REAL csetemp233 = tdg4223 * xformL33;

    const CCTK_REAL csetemp234 = tdg4230 * xformL03;

    const CCTK_REAL csetemp235 = tdg4231 * xformL13;

    const CCTK_REAL csetemp236 = tdg4232 * xformL23;

    const CCTK_REAL csetemp237 = tdg4233 * xformL33;

    const CCTK_REAL csetemp238 = tdg4330 * xformL03;

    const CCTK_REAL csetemp239 = tdg4331 * xformL13;

    const CCTK_REAL csetemp240 = tdg4332 * xformL23;

    const CCTK_REAL csetemp241 = tdg4333 * xformL33;

    const CCTK_REAL dg4013 =
        xformL00 *
            ((csetemp202 + csetemp203 + csetemp204 + csetemp205) * xformL01 +
             (csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL11 +
             (csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL21 +
             (csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL31) +
        xformL10 *
            ((csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL01 +
             (csetemp218 + csetemp219 + csetemp220 + csetemp221) * xformL11 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL21 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL31) +
        xformL20 *
            ((csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL01 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL11 +
             (csetemp230 + csetemp231 + csetemp232 + csetemp233) * xformL21 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL31) +
        xformL30 *
            ((csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL01 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL11 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL21 +
             (csetemp238 + csetemp239 + csetemp240 + csetemp241) * xformL31);

    const CCTK_REAL dg4020 =
        xformL10 *
            ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL02 +
             (csetemp100 + csetemp101 + csetemp98 + csetemp99) * xformL12 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL22 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL32) +
        xformL20 *
            ((csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL02 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL12 +
             (csetemp110 + csetemp111 + csetemp112 + csetemp113) * xformL22 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL32) +
        xformL30 *
            ((csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL02 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL12 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL22 +
             (csetemp118 + csetemp119 + csetemp120 + csetemp121) * xformL32) +
        xformL00 * ((csetemp82 + csetemp83 + csetemp84 + csetemp85) * xformL02 +
                    (csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL12 +
                    (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL22 +
                    (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL32);

    const CCTK_REAL dg4021 =
        xformL00 *
            ((csetemp122 + csetemp123 + csetemp124 + csetemp125) * xformL02 +
             (csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL12 +
             (csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL22 +
             (csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL32) +
        xformL10 *
            ((csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL02 +
             (csetemp138 + csetemp139 + csetemp140 + csetemp141) * xformL12 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL22 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL32) +
        xformL20 *
            ((csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL02 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL12 +
             (csetemp150 + csetemp151 + csetemp152 + csetemp153) * xformL22 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL32) +
        xformL30 *
            ((csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL02 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL12 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL22 +
             (csetemp158 + csetemp159 + csetemp160 + csetemp161) * xformL32);

    const CCTK_REAL dg4022 =
        xformL00 *
            ((csetemp162 + csetemp163 + csetemp164 + csetemp165) * xformL02 +
             (csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL12 +
             (csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL22 +
             (csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL32) +
        xformL10 *
            ((csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL02 +
             (csetemp178 + csetemp179 + csetemp180 + csetemp181) * xformL12 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL22 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL32) +
        xformL20 *
            ((csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL02 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL12 +
             (csetemp190 + csetemp191 + csetemp192 + csetemp193) * xformL22 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL32) +
        xformL30 *
            ((csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL02 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL12 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL22 +
             (csetemp198 + csetemp199 + csetemp200 + csetemp201) * xformL32);

    const CCTK_REAL dg4023 =
        xformL00 *
            ((csetemp202 + csetemp203 + csetemp204 + csetemp205) * xformL02 +
             (csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL12 +
             (csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL22 +
             (csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL32) +
        xformL10 *
            ((csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL02 +
             (csetemp218 + csetemp219 + csetemp220 + csetemp221) * xformL12 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL22 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL32) +
        xformL20 *
            ((csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL02 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL12 +
             (csetemp230 + csetemp231 + csetemp232 + csetemp233) * xformL22 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL32) +
        xformL30 *
            ((csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL02 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL12 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL22 +
             (csetemp238 + csetemp239 + csetemp240 + csetemp241) * xformL32);

    const CCTK_REAL dg4030 =
        xformL10 *
            ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL03 +
             (csetemp100 + csetemp101 + csetemp98 + csetemp99) * xformL13 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL23 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL33) +
        xformL20 *
            ((csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL03 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL13 +
             (csetemp110 + csetemp111 + csetemp112 + csetemp113) * xformL23 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL33) +
        xformL30 *
            ((csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL03 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL13 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL23 +
             (csetemp118 + csetemp119 + csetemp120 + csetemp121) * xformL33) +
        xformL00 * ((csetemp82 + csetemp83 + csetemp84 + csetemp85) * xformL03 +
                    (csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL13 +
                    (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL23 +
                    (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL33);

    const CCTK_REAL dg4031 =
        xformL00 *
            ((csetemp122 + csetemp123 + csetemp124 + csetemp125) * xformL03 +
             (csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL13 +
             (csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL23 +
             (csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL33) +
        xformL10 *
            ((csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL03 +
             (csetemp138 + csetemp139 + csetemp140 + csetemp141) * xformL13 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL23 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL33) +
        xformL20 *
            ((csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL03 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL13 +
             (csetemp150 + csetemp151 + csetemp152 + csetemp153) * xformL23 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL33) +
        xformL30 *
            ((csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL03 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL13 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL23 +
             (csetemp158 + csetemp159 + csetemp160 + csetemp161) * xformL33);

    const CCTK_REAL dg4032 =
        xformL00 *
            ((csetemp162 + csetemp163 + csetemp164 + csetemp165) * xformL03 +
             (csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL13 +
             (csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL23 +
             (csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL33) +
        xformL10 *
            ((csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL03 +
             (csetemp178 + csetemp179 + csetemp180 + csetemp181) * xformL13 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL23 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL33) +
        xformL20 *
            ((csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL03 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL13 +
             (csetemp190 + csetemp191 + csetemp192 + csetemp193) * xformL23 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL33) +
        xformL30 *
            ((csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL03 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL13 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL23 +
             (csetemp198 + csetemp199 + csetemp200 + csetemp201) * xformL33);

    const CCTK_REAL dg4033 =
        xformL00 *
            ((csetemp202 + csetemp203 + csetemp204 + csetemp205) * xformL03 +
             (csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL13 +
             (csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL23 +
             (csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL33) +
        xformL10 *
            ((csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL03 +
             (csetemp218 + csetemp219 + csetemp220 + csetemp221) * xformL13 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL23 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL33) +
        xformL20 *
            ((csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL03 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL13 +
             (csetemp230 + csetemp231 + csetemp232 + csetemp233) * xformL23 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL33) +
        xformL30 *
            ((csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL03 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL13 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL23 +
             (csetemp238 + csetemp239 + csetemp240 + csetemp241) * xformL33);

    const CCTK_REAL dg4110 =
        2 * ((csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL21 *
                 xformL31 +
             xformL11 * ((csetemp102 + csetemp103 + csetemp104 + csetemp105) *
                             xformL21 +
                         (csetemp106 + csetemp107 + csetemp108 + csetemp109) *
                             xformL31) +
             xformL01 *
                 ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL11 +
                  (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL21 +
                  (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL31)) +
        (csetemp82 + csetemp83 + csetemp84 + csetemp85) * pow(xformL01, 2) +
        (csetemp100 + csetemp101 + csetemp98 + csetemp99) * pow(xformL11, 2) +
        (csetemp110 + csetemp111 + csetemp112 + csetemp113) * pow(xformL21, 2) +
        (csetemp118 + csetemp119 + csetemp120 + csetemp121) * pow(xformL31, 2);

    const CCTK_REAL dg4111 =
        2 * ((csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL21 *
                 xformL31 +
             xformL01 * ((csetemp126 + csetemp127 + csetemp128 + csetemp129) *
                             xformL11 +
                         (csetemp130 + csetemp131 + csetemp132 + csetemp133) *
                             xformL21 +
                         (csetemp134 + csetemp135 + csetemp136 + csetemp137) *
                             xformL31) +
             xformL11 * ((csetemp142 + csetemp143 + csetemp144 + csetemp145) *
                             xformL21 +
                         (csetemp146 + csetemp147 + csetemp148 + csetemp149) *
                             xformL31)) +
        (csetemp122 + csetemp123 + csetemp124 + csetemp125) * pow(xformL01, 2) +
        (csetemp138 + csetemp139 + csetemp140 + csetemp141) * pow(xformL11, 2) +
        (csetemp150 + csetemp151 + csetemp152 + csetemp153) * pow(xformL21, 2) +
        (csetemp158 + csetemp159 + csetemp160 + csetemp161) * pow(xformL31, 2);

    const CCTK_REAL dg4112 =
        2 * ((csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL21 *
                 xformL31 +
             xformL01 * ((csetemp166 + csetemp167 + csetemp168 + csetemp169) *
                             xformL11 +
                         (csetemp170 + csetemp171 + csetemp172 + csetemp173) *
                             xformL21 +
                         (csetemp174 + csetemp175 + csetemp176 + csetemp177) *
                             xformL31) +
             xformL11 * ((csetemp182 + csetemp183 + csetemp184 + csetemp185) *
                             xformL21 +
                         (csetemp186 + csetemp187 + csetemp188 + csetemp189) *
                             xformL31)) +
        (csetemp162 + csetemp163 + csetemp164 + csetemp165) * pow(xformL01, 2) +
        (csetemp178 + csetemp179 + csetemp180 + csetemp181) * pow(xformL11, 2) +
        (csetemp190 + csetemp191 + csetemp192 + csetemp193) * pow(xformL21, 2) +
        (csetemp198 + csetemp199 + csetemp200 + csetemp201) * pow(xformL31, 2);

    const CCTK_REAL dg4113 =
        2 * ((csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL21 *
                 xformL31 +
             xformL01 * ((csetemp206 + csetemp207 + csetemp208 + csetemp209) *
                             xformL11 +
                         (csetemp210 + csetemp211 + csetemp212 + csetemp213) *
                             xformL21 +
                         (csetemp214 + csetemp215 + csetemp216 + csetemp217) *
                             xformL31) +
             xformL11 * ((csetemp222 + csetemp223 + csetemp224 + csetemp225) *
                             xformL21 +
                         (csetemp226 + csetemp227 + csetemp228 + csetemp229) *
                             xformL31)) +
        (csetemp202 + csetemp203 + csetemp204 + csetemp205) * pow(xformL01, 2) +
        (csetemp218 + csetemp219 + csetemp220 + csetemp221) * pow(xformL11, 2) +
        (csetemp230 + csetemp231 + csetemp232 + csetemp233) * pow(xformL21, 2) +
        (csetemp238 + csetemp239 + csetemp240 + csetemp241) * pow(xformL31, 2);

    const CCTK_REAL dg4120 =
        xformL11 *
            ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL02 +
             (csetemp100 + csetemp101 + csetemp98 + csetemp99) * xformL12 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL22 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL32) +
        xformL21 *
            ((csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL02 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL12 +
             (csetemp110 + csetemp111 + csetemp112 + csetemp113) * xformL22 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL32) +
        xformL31 *
            ((csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL02 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL12 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL22 +
             (csetemp118 + csetemp119 + csetemp120 + csetemp121) * xformL32) +
        xformL01 * ((csetemp82 + csetemp83 + csetemp84 + csetemp85) * xformL02 +
                    (csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL12 +
                    (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL22 +
                    (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL32);

    const CCTK_REAL dg4121 =
        xformL01 *
            ((csetemp122 + csetemp123 + csetemp124 + csetemp125) * xformL02 +
             (csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL12 +
             (csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL22 +
             (csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL32) +
        xformL11 *
            ((csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL02 +
             (csetemp138 + csetemp139 + csetemp140 + csetemp141) * xformL12 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL22 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL32) +
        xformL21 *
            ((csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL02 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL12 +
             (csetemp150 + csetemp151 + csetemp152 + csetemp153) * xformL22 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL32) +
        xformL31 *
            ((csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL02 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL12 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL22 +
             (csetemp158 + csetemp159 + csetemp160 + csetemp161) * xformL32);

    const CCTK_REAL dg4122 =
        xformL01 *
            ((csetemp162 + csetemp163 + csetemp164 + csetemp165) * xformL02 +
             (csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL12 +
             (csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL22 +
             (csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL32) +
        xformL11 *
            ((csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL02 +
             (csetemp178 + csetemp179 + csetemp180 + csetemp181) * xformL12 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL22 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL32) +
        xformL21 *
            ((csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL02 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL12 +
             (csetemp190 + csetemp191 + csetemp192 + csetemp193) * xformL22 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL32) +
        xformL31 *
            ((csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL02 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL12 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL22 +
             (csetemp198 + csetemp199 + csetemp200 + csetemp201) * xformL32);

    const CCTK_REAL dg4123 =
        xformL01 *
            ((csetemp202 + csetemp203 + csetemp204 + csetemp205) * xformL02 +
             (csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL12 +
             (csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL22 +
             (csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL32) +
        xformL11 *
            ((csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL02 +
             (csetemp218 + csetemp219 + csetemp220 + csetemp221) * xformL12 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL22 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL32) +
        xformL21 *
            ((csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL02 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL12 +
             (csetemp230 + csetemp231 + csetemp232 + csetemp233) * xformL22 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL32) +
        xformL31 *
            ((csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL02 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL12 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL22 +
             (csetemp238 + csetemp239 + csetemp240 + csetemp241) * xformL32);

    const CCTK_REAL dg4130 =
        xformL11 *
            ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL03 +
             (csetemp100 + csetemp101 + csetemp98 + csetemp99) * xformL13 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL23 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL33) +
        xformL21 *
            ((csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL03 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL13 +
             (csetemp110 + csetemp111 + csetemp112 + csetemp113) * xformL23 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL33) +
        xformL31 *
            ((csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL03 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL13 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL23 +
             (csetemp118 + csetemp119 + csetemp120 + csetemp121) * xformL33) +
        xformL01 * ((csetemp82 + csetemp83 + csetemp84 + csetemp85) * xformL03 +
                    (csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL13 +
                    (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL23 +
                    (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL33);

    const CCTK_REAL dg4131 =
        xformL01 *
            ((csetemp122 + csetemp123 + csetemp124 + csetemp125) * xformL03 +
             (csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL13 +
             (csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL23 +
             (csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL33) +
        xformL11 *
            ((csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL03 +
             (csetemp138 + csetemp139 + csetemp140 + csetemp141) * xformL13 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL23 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL33) +
        xformL21 *
            ((csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL03 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL13 +
             (csetemp150 + csetemp151 + csetemp152 + csetemp153) * xformL23 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL33) +
        xformL31 *
            ((csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL03 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL13 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL23 +
             (csetemp158 + csetemp159 + csetemp160 + csetemp161) * xformL33);

    const CCTK_REAL dg4132 =
        xformL01 *
            ((csetemp162 + csetemp163 + csetemp164 + csetemp165) * xformL03 +
             (csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL13 +
             (csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL23 +
             (csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL33) +
        xformL11 *
            ((csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL03 +
             (csetemp178 + csetemp179 + csetemp180 + csetemp181) * xformL13 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL23 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL33) +
        xformL21 *
            ((csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL03 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL13 +
             (csetemp190 + csetemp191 + csetemp192 + csetemp193) * xformL23 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL33) +
        xformL31 *
            ((csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL03 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL13 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL23 +
             (csetemp198 + csetemp199 + csetemp200 + csetemp201) * xformL33);

    const CCTK_REAL dg4133 =
        xformL01 *
            ((csetemp202 + csetemp203 + csetemp204 + csetemp205) * xformL03 +
             (csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL13 +
             (csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL23 +
             (csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL33) +
        xformL11 *
            ((csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL03 +
             (csetemp218 + csetemp219 + csetemp220 + csetemp221) * xformL13 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL23 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL33) +
        xformL21 *
            ((csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL03 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL13 +
             (csetemp230 + csetemp231 + csetemp232 + csetemp233) * xformL23 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL33) +
        xformL31 *
            ((csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL03 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL13 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL23 +
             (csetemp238 + csetemp239 + csetemp240 + csetemp241) * xformL33);

    const CCTK_REAL dg4220 =
        2 * ((csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL22 *
                 xformL32 +
             xformL12 * ((csetemp102 + csetemp103 + csetemp104 + csetemp105) *
                             xformL22 +
                         (csetemp106 + csetemp107 + csetemp108 + csetemp109) *
                             xformL32) +
             xformL02 *
                 ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL12 +
                  (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL22 +
                  (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL32)) +
        (csetemp82 + csetemp83 + csetemp84 + csetemp85) * pow(xformL02, 2) +
        (csetemp100 + csetemp101 + csetemp98 + csetemp99) * pow(xformL12, 2) +
        (csetemp110 + csetemp111 + csetemp112 + csetemp113) * pow(xformL22, 2) +
        (csetemp118 + csetemp119 + csetemp120 + csetemp121) * pow(xformL32, 2);

    const CCTK_REAL dg4221 =
        2 * ((csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL22 *
                 xformL32 +
             xformL02 * ((csetemp126 + csetemp127 + csetemp128 + csetemp129) *
                             xformL12 +
                         (csetemp130 + csetemp131 + csetemp132 + csetemp133) *
                             xformL22 +
                         (csetemp134 + csetemp135 + csetemp136 + csetemp137) *
                             xformL32) +
             xformL12 * ((csetemp142 + csetemp143 + csetemp144 + csetemp145) *
                             xformL22 +
                         (csetemp146 + csetemp147 + csetemp148 + csetemp149) *
                             xformL32)) +
        (csetemp122 + csetemp123 + csetemp124 + csetemp125) * pow(xformL02, 2) +
        (csetemp138 + csetemp139 + csetemp140 + csetemp141) * pow(xformL12, 2) +
        (csetemp150 + csetemp151 + csetemp152 + csetemp153) * pow(xformL22, 2) +
        (csetemp158 + csetemp159 + csetemp160 + csetemp161) * pow(xformL32, 2);

    const CCTK_REAL dg4222 =
        2 * ((csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL22 *
                 xformL32 +
             xformL02 * ((csetemp166 + csetemp167 + csetemp168 + csetemp169) *
                             xformL12 +
                         (csetemp170 + csetemp171 + csetemp172 + csetemp173) *
                             xformL22 +
                         (csetemp174 + csetemp175 + csetemp176 + csetemp177) *
                             xformL32) +
             xformL12 * ((csetemp182 + csetemp183 + csetemp184 + csetemp185) *
                             xformL22 +
                         (csetemp186 + csetemp187 + csetemp188 + csetemp189) *
                             xformL32)) +
        (csetemp162 + csetemp163 + csetemp164 + csetemp165) * pow(xformL02, 2) +
        (csetemp178 + csetemp179 + csetemp180 + csetemp181) * pow(xformL12, 2) +
        (csetemp190 + csetemp191 + csetemp192 + csetemp193) * pow(xformL22, 2) +
        (csetemp198 + csetemp199 + csetemp200 + csetemp201) * pow(xformL32, 2);

    const CCTK_REAL dg4223 =
        2 * ((csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL22 *
                 xformL32 +
             xformL02 * ((csetemp206 + csetemp207 + csetemp208 + csetemp209) *
                             xformL12 +
                         (csetemp210 + csetemp211 + csetemp212 + csetemp213) *
                             xformL22 +
                         (csetemp214 + csetemp215 + csetemp216 + csetemp217) *
                             xformL32) +
             xformL12 * ((csetemp222 + csetemp223 + csetemp224 + csetemp225) *
                             xformL22 +
                         (csetemp226 + csetemp227 + csetemp228 + csetemp229) *
                             xformL32)) +
        (csetemp202 + csetemp203 + csetemp204 + csetemp205) * pow(xformL02, 2) +
        (csetemp218 + csetemp219 + csetemp220 + csetemp221) * pow(xformL12, 2) +
        (csetemp230 + csetemp231 + csetemp232 + csetemp233) * pow(xformL22, 2) +
        (csetemp238 + csetemp239 + csetemp240 + csetemp241) * pow(xformL32, 2);

    const CCTK_REAL dg4230 =
        xformL12 *
            ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL03 +
             (csetemp100 + csetemp101 + csetemp98 + csetemp99) * xformL13 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL23 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL33) +
        xformL22 *
            ((csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL03 +
             (csetemp102 + csetemp103 + csetemp104 + csetemp105) * xformL13 +
             (csetemp110 + csetemp111 + csetemp112 + csetemp113) * xformL23 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL33) +
        xformL32 *
            ((csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL03 +
             (csetemp106 + csetemp107 + csetemp108 + csetemp109) * xformL13 +
             (csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL23 +
             (csetemp118 + csetemp119 + csetemp120 + csetemp121) * xformL33) +
        xformL02 * ((csetemp82 + csetemp83 + csetemp84 + csetemp85) * xformL03 +
                    (csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL13 +
                    (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL23 +
                    (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL33);

    const CCTK_REAL dg4231 =
        xformL02 *
            ((csetemp122 + csetemp123 + csetemp124 + csetemp125) * xformL03 +
             (csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL13 +
             (csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL23 +
             (csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL33) +
        xformL12 *
            ((csetemp126 + csetemp127 + csetemp128 + csetemp129) * xformL03 +
             (csetemp138 + csetemp139 + csetemp140 + csetemp141) * xformL13 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL23 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL33) +
        xformL22 *
            ((csetemp130 + csetemp131 + csetemp132 + csetemp133) * xformL03 +
             (csetemp142 + csetemp143 + csetemp144 + csetemp145) * xformL13 +
             (csetemp150 + csetemp151 + csetemp152 + csetemp153) * xformL23 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL33) +
        xformL32 *
            ((csetemp134 + csetemp135 + csetemp136 + csetemp137) * xformL03 +
             (csetemp146 + csetemp147 + csetemp148 + csetemp149) * xformL13 +
             (csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL23 +
             (csetemp158 + csetemp159 + csetemp160 + csetemp161) * xformL33);

    const CCTK_REAL dg4232 =
        xformL02 *
            ((csetemp162 + csetemp163 + csetemp164 + csetemp165) * xformL03 +
             (csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL13 +
             (csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL23 +
             (csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL33) +
        xformL12 *
            ((csetemp166 + csetemp167 + csetemp168 + csetemp169) * xformL03 +
             (csetemp178 + csetemp179 + csetemp180 + csetemp181) * xformL13 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL23 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL33) +
        xformL22 *
            ((csetemp170 + csetemp171 + csetemp172 + csetemp173) * xformL03 +
             (csetemp182 + csetemp183 + csetemp184 + csetemp185) * xformL13 +
             (csetemp190 + csetemp191 + csetemp192 + csetemp193) * xformL23 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL33) +
        xformL32 *
            ((csetemp174 + csetemp175 + csetemp176 + csetemp177) * xformL03 +
             (csetemp186 + csetemp187 + csetemp188 + csetemp189) * xformL13 +
             (csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL23 +
             (csetemp198 + csetemp199 + csetemp200 + csetemp201) * xformL33);

    const CCTK_REAL dg4233 =
        xformL02 *
            ((csetemp202 + csetemp203 + csetemp204 + csetemp205) * xformL03 +
             (csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL13 +
             (csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL23 +
             (csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL33) +
        xformL12 *
            ((csetemp206 + csetemp207 + csetemp208 + csetemp209) * xformL03 +
             (csetemp218 + csetemp219 + csetemp220 + csetemp221) * xformL13 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL23 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL33) +
        xformL22 *
            ((csetemp210 + csetemp211 + csetemp212 + csetemp213) * xformL03 +
             (csetemp222 + csetemp223 + csetemp224 + csetemp225) * xformL13 +
             (csetemp230 + csetemp231 + csetemp232 + csetemp233) * xformL23 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL33) +
        xformL32 *
            ((csetemp214 + csetemp215 + csetemp216 + csetemp217) * xformL03 +
             (csetemp226 + csetemp227 + csetemp228 + csetemp229) * xformL13 +
             (csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL23 +
             (csetemp238 + csetemp239 + csetemp240 + csetemp241) * xformL33);

    const CCTK_REAL dg4330 =
        2 * ((csetemp114 + csetemp115 + csetemp116 + csetemp117) * xformL23 *
                 xformL33 +
             xformL13 * ((csetemp102 + csetemp103 + csetemp104 + csetemp105) *
                             xformL23 +
                         (csetemp106 + csetemp107 + csetemp108 + csetemp109) *
                             xformL33) +
             xformL03 *
                 ((csetemp86 + csetemp87 + csetemp88 + csetemp89) * xformL13 +
                  (csetemp90 + csetemp91 + csetemp92 + csetemp93) * xformL23 +
                  (csetemp94 + csetemp95 + csetemp96 + csetemp97) * xformL33)) +
        (csetemp82 + csetemp83 + csetemp84 + csetemp85) * pow(xformL03, 2) +
        (csetemp100 + csetemp101 + csetemp98 + csetemp99) * pow(xformL13, 2) +
        (csetemp110 + csetemp111 + csetemp112 + csetemp113) * pow(xformL23, 2) +
        (csetemp118 + csetemp119 + csetemp120 + csetemp121) * pow(xformL33, 2);

    const CCTK_REAL dg4331 =
        2 * ((csetemp154 + csetemp155 + csetemp156 + csetemp157) * xformL23 *
                 xformL33 +
             xformL03 * ((csetemp126 + csetemp127 + csetemp128 + csetemp129) *
                             xformL13 +
                         (csetemp130 + csetemp131 + csetemp132 + csetemp133) *
                             xformL23 +
                         (csetemp134 + csetemp135 + csetemp136 + csetemp137) *
                             xformL33) +
             xformL13 * ((csetemp142 + csetemp143 + csetemp144 + csetemp145) *
                             xformL23 +
                         (csetemp146 + csetemp147 + csetemp148 + csetemp149) *
                             xformL33)) +
        (csetemp122 + csetemp123 + csetemp124 + csetemp125) * pow(xformL03, 2) +
        (csetemp138 + csetemp139 + csetemp140 + csetemp141) * pow(xformL13, 2) +
        (csetemp150 + csetemp151 + csetemp152 + csetemp153) * pow(xformL23, 2) +
        (csetemp158 + csetemp159 + csetemp160 + csetemp161) * pow(xformL33, 2);

    const CCTK_REAL dg4332 =
        2 * ((csetemp194 + csetemp195 + csetemp196 + csetemp197) * xformL23 *
                 xformL33 +
             xformL03 * ((csetemp166 + csetemp167 + csetemp168 + csetemp169) *
                             xformL13 +
                         (csetemp170 + csetemp171 + csetemp172 + csetemp173) *
                             xformL23 +
                         (csetemp174 + csetemp175 + csetemp176 + csetemp177) *
                             xformL33) +
             xformL13 * ((csetemp182 + csetemp183 + csetemp184 + csetemp185) *
                             xformL23 +
                         (csetemp186 + csetemp187 + csetemp188 + csetemp189) *
                             xformL33)) +
        (csetemp162 + csetemp163 + csetemp164 + csetemp165) * pow(xformL03, 2) +
        (csetemp178 + csetemp179 + csetemp180 + csetemp181) * pow(xformL13, 2) +
        (csetemp190 + csetemp191 + csetemp192 + csetemp193) * pow(xformL23, 2) +
        (csetemp198 + csetemp199 + csetemp200 + csetemp201) * pow(xformL33, 2);

    const CCTK_REAL dg4333 =
        2 * ((csetemp234 + csetemp235 + csetemp236 + csetemp237) * xformL23 *
                 xformL33 +
             xformL03 * ((csetemp206 + csetemp207 + csetemp208 + csetemp209) *
                             xformL13 +
                         (csetemp210 + csetemp211 + csetemp212 + csetemp213) *
                             xformL23 +
                         (csetemp214 + csetemp215 + csetemp216 + csetemp217) *
                             xformL33) +
             xformL13 * ((csetemp222 + csetemp223 + csetemp224 + csetemp225) *
                             xformL23 +
                         (csetemp226 + csetemp227 + csetemp228 + csetemp229) *
                             xformL33)) +
        (csetemp202 + csetemp203 + csetemp204 + csetemp205) * pow(xformL03, 2) +
        (csetemp218 + csetemp219 + csetemp220 + csetemp221) * pow(xformL13, 2) +
        (csetemp230 + csetemp231 + csetemp232 + csetemp233) * pow(xformL23, 2) +
        (csetemp238 + csetemp239 + csetemp240 + csetemp241) * pow(xformL33, 2);

    const CCTK_REAL betal1 = g401;

    const CCTK_REAL betal2 = g402;

    const CCTK_REAL betal3 = g403;

    gxxL = g411;

    gxyL = g412;

    gxzL = g413;

    gyyL = g422;

    gyzL = g423;

    gzzL = g433;

    const CCTK_REAL csetemp242 = pow(gxzL, 2);

    const CCTK_REAL csetemp243 = pow(gyzL, 2);

    const CCTK_REAL csetemp244 = pow(gxyL, 2);

    const CCTK_REAL detg = 2 * gxyL * gxzL * gyzL +
                           gyyL * (gxxL * gzzL - csetemp242) -
                           gxxL * csetemp243 - gzzL * csetemp244;

    const CCTK_REAL csetemp245 = pow(detg, -1);

    const CCTK_REAL gu11 = (gyyL * gzzL - csetemp243) * csetemp245;

    const CCTK_REAL gu12 = (gxzL * gyzL - gxyL * gzzL) * csetemp245;

    const CCTK_REAL gu13 = (-(gxzL * gyyL) + gxyL * gyzL) * csetemp245;

    const CCTK_REAL gu22 = (gxxL * gzzL - csetemp242) * csetemp245;

    const CCTK_REAL gu23 = (gxyL * gxzL - gxxL * gyzL) * csetemp245;

    const CCTK_REAL gu33 = (gxxL * gyyL - csetemp244) * csetemp245;

    betaxL = betal1 * gu11 + betal2 * gu12 + betal3 * gu13;

    betayL = betal1 * gu12 + betal2 * gu22 + betal3 * gu23;

    betazL = betal1 * gu13 + betal2 * gu23 + betal3 * gu33;

    const CCTK_REAL betasq =
        betaxL * betal1 + betayL * betal2 + betazL * betal3;

    alpL = pow(betasq - g400, 0.5);

    const CCTK_REAL dtg11 = dg4110;

    const CCTK_REAL dtg12 = dg4120;

    const CCTK_REAL dtg13 = dg4130;

    const CCTK_REAL dtg22 = dg4220;

    const CCTK_REAL dtg23 = dg4230;

    const CCTK_REAL dtg33 = dg4330;

    const CCTK_REAL dg111 = dg4111;

    const CCTK_REAL dg112 = dg4112;

    const CCTK_REAL dg113 = dg4113;

    const CCTK_REAL dg121 = dg4121;

    const CCTK_REAL dg122 = dg4122;

    const CCTK_REAL dg123 = dg4123;

    const CCTK_REAL dg131 = dg4131;

    const CCTK_REAL dg132 = dg4132;

    const CCTK_REAL dg133 = dg4133;

    const CCTK_REAL dg221 = dg4221;

    const CCTK_REAL dg222 = dg4222;

    const CCTK_REAL dg223 = dg4223;

    const CCTK_REAL dg231 = dg4231;

    const CCTK_REAL dg232 = dg4232;

    const CCTK_REAL dg233 = dg4233;

    const CCTK_REAL dg331 = dg4331;

    const CCTK_REAL dg332 = dg4332;

    const CCTK_REAL dg333 = dg4333;

    const CCTK_REAL csetemp246 = dtg11 * gu11;

    const CCTK_REAL csetemp247 = dtg12 * gu12;

    const CCTK_REAL csetemp248 = dtg13 * gu13;

    const CCTK_REAL csetemp249 = dtg12 * gu11;

    const CCTK_REAL csetemp250 = dtg22 * gu12;

    const CCTK_REAL csetemp251 = dtg23 * gu13;

    const CCTK_REAL csetemp252 = dtg13 * gu11;

    const CCTK_REAL csetemp253 = dtg23 * gu12;

    const CCTK_REAL csetemp254 = dtg33 * gu13;

    const CCTK_REAL dtgu11 = -((csetemp246 + csetemp247 + csetemp248) * gu11) -
                             (csetemp249 + csetemp250 + csetemp251) * gu12 -
                             (csetemp252 + csetemp253 + csetemp254) * gu13;

    const CCTK_REAL dtgu12 = -((csetemp246 + csetemp247 + csetemp248) * gu12) -
                             (csetemp249 + csetemp250 + csetemp251) * gu22 -
                             (csetemp252 + csetemp253 + csetemp254) * gu23;

    const CCTK_REAL dtgu13 = -((csetemp246 + csetemp247 + csetemp248) * gu13) -
                             (csetemp249 + csetemp250 + csetemp251) * gu23 -
                             (csetemp252 + csetemp253 + csetemp254) * gu33;

    const CCTK_REAL csetemp255 = dtg11 * gu12;

    const CCTK_REAL csetemp256 = dtg12 * gu22;

    const CCTK_REAL csetemp257 = dtg13 * gu23;

    const CCTK_REAL csetemp258 = dtg22 * gu22;

    const CCTK_REAL csetemp259 = dtg23 * gu23;

    const CCTK_REAL csetemp260 = dtg13 * gu12;

    const CCTK_REAL csetemp261 = dtg23 * gu22;

    const CCTK_REAL csetemp262 = dtg33 * gu23;

    const CCTK_REAL dtgu22 = -((csetemp255 + csetemp256 + csetemp257) * gu12) -
                             (csetemp247 + csetemp258 + csetemp259) * gu22 -
                             (csetemp260 + csetemp261 + csetemp262) * gu23;

    const CCTK_REAL dtgu23 = -((csetemp255 + csetemp256 + csetemp257) * gu13) -
                             (csetemp247 + csetemp258 + csetemp259) * gu23 -
                             (csetemp260 + csetemp261 + csetemp262) * gu33;

    const CCTK_REAL dtgu33 =
        -(gu13 * (dtg11 * gu13 + dtg12 * gu23 + dtg13 * gu33)) -
        gu23 * (dtg12 * gu13 + dtg22 * gu23 + dtg23 * gu33) -
        gu33 * (csetemp248 + csetemp259 + dtg33 * gu33);

    const CCTK_REAL csetemp263 = dg111 * gu11;

    const CCTK_REAL csetemp264 = dg121 * gu12;

    const CCTK_REAL csetemp265 = dg131 * gu13;

    const CCTK_REAL csetemp266 = dg121 * gu11;

    const CCTK_REAL csetemp267 = dg221 * gu12;

    const CCTK_REAL csetemp268 = dg231 * gu13;

    const CCTK_REAL csetemp269 = dg131 * gu11;

    const CCTK_REAL csetemp270 = dg231 * gu12;

    const CCTK_REAL csetemp271 = dg331 * gu13;

    const CCTK_REAL dgu111 = -((csetemp263 + csetemp264 + csetemp265) * gu11) -
                             (csetemp266 + csetemp267 + csetemp268) * gu12 -
                             (csetemp269 + csetemp270 + csetemp271) * gu13;

    const CCTK_REAL dgu121 = -((csetemp263 + csetemp264 + csetemp265) * gu12) -
                             (csetemp266 + csetemp267 + csetemp268) * gu22 -
                             (csetemp269 + csetemp270 + csetemp271) * gu23;

    const CCTK_REAL dgu131 = -((csetemp263 + csetemp264 + csetemp265) * gu13) -
                             (csetemp266 + csetemp267 + csetemp268) * gu23 -
                             (csetemp269 + csetemp270 + csetemp271) * gu33;

    const CCTK_REAL csetemp272 = dg111 * gu12;

    const CCTK_REAL csetemp273 = dg121 * gu22;

    const CCTK_REAL csetemp274 = dg131 * gu23;

    const CCTK_REAL csetemp275 = dg221 * gu22;

    const CCTK_REAL csetemp276 = dg231 * gu23;

    const CCTK_REAL csetemp277 = dg131 * gu12;

    const CCTK_REAL csetemp278 = dg231 * gu22;

    const CCTK_REAL csetemp279 = dg331 * gu23;

    const CCTK_REAL dgu221 = -((csetemp272 + csetemp273 + csetemp274) * gu12) -
                             (csetemp264 + csetemp275 + csetemp276) * gu22 -
                             (csetemp277 + csetemp278 + csetemp279) * gu23;

    const CCTK_REAL dgu231 = -((csetemp272 + csetemp273 + csetemp274) * gu13) -
                             (csetemp264 + csetemp275 + csetemp276) * gu23 -
                             (csetemp277 + csetemp278 + csetemp279) * gu33;

    const CCTK_REAL dgu331 =
        -(gu13 * (dg111 * gu13 + dg121 * gu23 + dg131 * gu33)) -
        gu23 * (dg121 * gu13 + dg221 * gu23 + dg231 * gu33) -
        gu33 * (csetemp265 + csetemp276 + dg331 * gu33);

    const CCTK_REAL csetemp280 = dg112 * gu11;

    const CCTK_REAL csetemp281 = dg122 * gu12;

    const CCTK_REAL csetemp282 = dg132 * gu13;

    const CCTK_REAL csetemp283 = dg122 * gu11;

    const CCTK_REAL csetemp284 = dg222 * gu12;

    const CCTK_REAL csetemp285 = dg232 * gu13;

    const CCTK_REAL csetemp286 = dg132 * gu11;

    const CCTK_REAL csetemp287 = dg232 * gu12;

    const CCTK_REAL csetemp288 = dg332 * gu13;

    const CCTK_REAL dgu112 = -((csetemp280 + csetemp281 + csetemp282) * gu11) -
                             (csetemp283 + csetemp284 + csetemp285) * gu12 -
                             (csetemp286 + csetemp287 + csetemp288) * gu13;

    const CCTK_REAL dgu122 = -((csetemp280 + csetemp281 + csetemp282) * gu12) -
                             (csetemp283 + csetemp284 + csetemp285) * gu22 -
                             (csetemp286 + csetemp287 + csetemp288) * gu23;

    const CCTK_REAL dgu132 = -((csetemp280 + csetemp281 + csetemp282) * gu13) -
                             (csetemp283 + csetemp284 + csetemp285) * gu23 -
                             (csetemp286 + csetemp287 + csetemp288) * gu33;

    const CCTK_REAL csetemp289 = dg112 * gu12;

    const CCTK_REAL csetemp290 = dg122 * gu22;

    const CCTK_REAL csetemp291 = dg132 * gu23;

    const CCTK_REAL csetemp292 = dg222 * gu22;

    const CCTK_REAL csetemp293 = dg232 * gu23;

    const CCTK_REAL csetemp294 = dg132 * gu12;

    const CCTK_REAL csetemp295 = dg232 * gu22;

    const CCTK_REAL csetemp296 = dg332 * gu23;

    const CCTK_REAL dgu222 = -((csetemp289 + csetemp290 + csetemp291) * gu12) -
                             (csetemp281 + csetemp292 + csetemp293) * gu22 -
                             (csetemp294 + csetemp295 + csetemp296) * gu23;

    const CCTK_REAL dgu232 = -((csetemp289 + csetemp290 + csetemp291) * gu13) -
                             (csetemp281 + csetemp292 + csetemp293) * gu23 -
                             (csetemp294 + csetemp295 + csetemp296) * gu33;

    const CCTK_REAL dgu332 =
        -(gu13 * (dg112 * gu13 + dg122 * gu23 + dg132 * gu33)) -
        gu23 * (dg122 * gu13 + dg222 * gu23 + dg232 * gu33) -
        gu33 * (csetemp282 + csetemp293 + dg332 * gu33);

    const CCTK_REAL csetemp297 = dg113 * gu11;

    const CCTK_REAL csetemp298 = dg123 * gu12;

    const CCTK_REAL csetemp299 = dg133 * gu13;

    const CCTK_REAL csetemp300 = dg123 * gu11;

    const CCTK_REAL csetemp301 = dg223 * gu12;

    const CCTK_REAL csetemp302 = dg233 * gu13;

    const CCTK_REAL csetemp303 = dg133 * gu11;

    const CCTK_REAL csetemp304 = dg233 * gu12;

    const CCTK_REAL csetemp305 = dg333 * gu13;

    const CCTK_REAL dgu113 = -((csetemp297 + csetemp298 + csetemp299) * gu11) -
                             (csetemp300 + csetemp301 + csetemp302) * gu12 -
                             (csetemp303 + csetemp304 + csetemp305) * gu13;

    const CCTK_REAL dgu123 = -((csetemp297 + csetemp298 + csetemp299) * gu12) -
                             (csetemp300 + csetemp301 + csetemp302) * gu22 -
                             (csetemp303 + csetemp304 + csetemp305) * gu23;

    const CCTK_REAL dgu133 = -((csetemp297 + csetemp298 + csetemp299) * gu13) -
                             (csetemp300 + csetemp301 + csetemp302) * gu23 -
                             (csetemp303 + csetemp304 + csetemp305) * gu33;

    const CCTK_REAL csetemp306 = dg113 * gu12;

    const CCTK_REAL csetemp307 = dg123 * gu22;

    const CCTK_REAL csetemp308 = dg133 * gu23;

    const CCTK_REAL csetemp309 = dg223 * gu22;

    const CCTK_REAL csetemp310 = dg233 * gu23;

    const CCTK_REAL csetemp311 = dg133 * gu12;

    const CCTK_REAL csetemp312 = dg233 * gu22;

    const CCTK_REAL csetemp313 = dg333 * gu23;

    const CCTK_REAL dgu223 = -((csetemp306 + csetemp307 + csetemp308) * gu12) -
                             (csetemp298 + csetemp309 + csetemp310) * gu22 -
                             (csetemp311 + csetemp312 + csetemp313) * gu23;

    const CCTK_REAL dgu233 = -((csetemp306 + csetemp307 + csetemp308) * gu13) -
                             (csetemp298 + csetemp309 + csetemp310) * gu23 -
                             (csetemp311 + csetemp312 + csetemp313) * gu33;

    const CCTK_REAL dgu333 =
        -(gu13 * (dg113 * gu13 + dg123 * gu23 + dg133 * gu33)) -
        gu23 * (dg123 * gu13 + dg223 * gu23 + dg233 * gu33) -
        gu33 * (csetemp299 + csetemp310 + dg333 * gu33);

    const CCTK_REAL dtbetal1 = dg4010;

    const CCTK_REAL dtbetal2 = dg4020;

    const CCTK_REAL dtbetal3 = dg4030;

    const CCTK_REAL dbetal11 = dg4011;

    const CCTK_REAL dbetal12 = dg4012;

    const CCTK_REAL dbetal13 = dg4013;

    const CCTK_REAL dbetal21 = dg4021;

    const CCTK_REAL dbetal22 = dg4022;

    const CCTK_REAL dbetal23 = dg4023;

    const CCTK_REAL dbetal31 = dg4031;

    const CCTK_REAL dbetal32 = dg4032;

    const CCTK_REAL dbetal33 = dg4033;

    dtbetaxL = betal1 * dtgu11 + betal2 * dtgu12 + betal3 * dtgu13 +
               dtbetal1 * gu11 + dtbetal2 * gu12 + dtbetal3 * gu13;

    dtbetayL = betal1 * dtgu12 + betal2 * dtgu22 + betal3 * dtgu23 +
               dtbetal1 * gu12 + dtbetal2 * gu22 + dtbetal3 * gu23;

    dtbetazL = betal1 * dtgu13 + betal2 * dtgu23 + betal3 * dtgu33 +
               dtbetal1 * gu13 + dtbetal2 * gu23 + dtbetal3 * gu33;

    const CCTK_REAL dbeta11 = betal1 * dgu111 + betal2 * dgu121 +
                              betal3 * dgu131 + dbetal11 * gu11 +
                              dbetal21 * gu12 + dbetal31 * gu13;

    const CCTK_REAL dbeta21 = betal1 * dgu121 + betal2 * dgu221 +
                              betal3 * dgu231 + dbetal11 * gu12 +
                              dbetal21 * gu22 + dbetal31 * gu23;

    const CCTK_REAL dbeta31 = betal1 * dgu131 + betal2 * dgu231 +
                              betal3 * dgu331 + dbetal11 * gu13 +
                              dbetal21 * gu23 + dbetal31 * gu33;

    const CCTK_REAL dbeta12 = betal1 * dgu112 + betal2 * dgu122 +
                              betal3 * dgu132 + dbetal12 * gu11 +
                              dbetal22 * gu12 + dbetal32 * gu13;

    const CCTK_REAL dbeta22 = betal1 * dgu122 + betal2 * dgu222 +
                              betal3 * dgu232 + dbetal12 * gu12 +
                              dbetal22 * gu22 + dbetal32 * gu23;

    const CCTK_REAL dbeta32 = betal1 * dgu132 + betal2 * dgu232 +
                              betal3 * dgu332 + dbetal12 * gu13 +
                              dbetal22 * gu23 + dbetal32 * gu33;

    const CCTK_REAL dbeta13 = betal1 * dgu113 + betal2 * dgu123 +
                              betal3 * dgu133 + dbetal13 * gu11 +
                              dbetal23 * gu12 + dbetal33 * gu13;

    const CCTK_REAL dbeta23 = betal1 * dgu123 + betal2 * dgu223 +
                              betal3 * dgu233 + dbetal13 * gu12 +
                              dbetal23 * gu22 + dbetal33 * gu23;

    const CCTK_REAL dbeta33 = betal1 * dgu133 + betal2 * dgu233 +
                              betal3 * dgu333 + dbetal13 * gu13 +
                              dbetal23 * gu23 + dbetal33 * gu33;

    const CCTK_REAL dtbetasq = dtbetaxL * betal1 + dtbetayL * betal2 +
                               dtbetazL * betal3 + betaxL * dtbetal1 +
                               betayL * dtbetal2 + betazL * dtbetal3;

    const CCTK_REAL csetemp314 = pow(alpL, -1);

    const CCTK_REAL dtalpL = 0.5 * csetemp314 * (-dg4000 + dtbetasq);

    const CCTK_REAL kxxL =
        0.5 * csetemp314 *
        (2 * (gxxL * dbeta11 + gxyL * dbeta21 + gxzL * dbeta31) +
         betaxL * dg111 + betayL * dg112 + betazL * dg113 - dtg11);

    const CCTK_REAL kxyL =
        0.5 * csetemp314 *
        (gxxL * dbeta12 + gyyL * dbeta21 + gxyL * (dbeta11 + dbeta22) +
         gyzL * dbeta31 + gxzL * dbeta32 + betaxL * dg121 + betayL * dg122 +
         betazL * dg123 - dtg12);

    const CCTK_REAL kxzL =
        0.5 * csetemp314 *
        (gxxL * dbeta13 + gyzL * dbeta21 + gxyL * dbeta23 + gzzL * dbeta31 +
         gxzL * (dbeta11 + dbeta33) + betaxL * dg131 + betayL * dg132 +
         betazL * dg133 - dtg13);

    const CCTK_REAL kyyL =
        0.5 * csetemp314 *
        (2 * (gxyL * dbeta12 + gyyL * dbeta22 + gyzL * dbeta32) +
         betaxL * dg221 + betayL * dg222 + betazL * dg223 - dtg22);

    const CCTK_REAL kyzL =
        0.5 * csetemp314 *
        (gxzL * dbeta12 + gxyL * dbeta13 + gyyL * dbeta23 + gzzL * dbeta32 +
         gyzL * (dbeta22 + dbeta33) + betaxL * dg231 + betayL * dg232 +
         betazL * dg233 - dtg23);

    const CCTK_REAL kzzL =
        0.5 * csetemp314 *
        (2 * (gxzL * dbeta13 + gyzL * dbeta23 + gzzL * dbeta33) +
         betaxL * dg331 + betayL * dg332 + betazL * dg333 - dtg33);

    alp_(p.I) = alpL;
    betax_(p.I) = betaxL;
    betay_(p.I) = betayL;
    betaz_(p.I) = betazL;
    dtalp_(p.I) = dtalpL;
    dtbetax_(p.I) = dtbetaxL;
    dtbetay_(p.I) = dtbetayL;
    dtbetaz_(p.I) = dtbetazL;
    gxx_(p.I) = gxxL;
    gxy_(p.I) = gxyL;
    gxz_(p.I) = gxzL;
    gyy_(p.I) = gyyL;
    gyz_(p.I) = gyzL;
    gzz_(p.I) = gzzL;
    kxx_(p.I) = kxxL;
    kxy_(p.I) = kxyL;
    kxz_(p.I) = kxzL;
    kyy_(p.I) = kyyL;
    kyz_(p.I) = kyzL;
    kzz_(p.I) = kzzL;
  };

  loop_all<0, 0, 0>(cctkGH, id_lambda);
}
