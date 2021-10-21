/*
 *  MinkowskiX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of MinkowskiX.
 *
 *  MinkowskiX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MinkowskiX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MinkowskiX.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Initial.cpp
 *  Fill the ADM grid functions with Minkowski initial data.
 */

#include "MinkowskiX.hpp"

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

extern "C" void MinkowskiX::MinkowskiX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MinkowskiX_Initial;
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

  const GF3D2<CCTK_REAL> betax_(layout, betax);
  const GF3D2<CCTK_REAL> betay_(layout, betay);
  const GF3D2<CCTK_REAL> betaz_(layout, betaz);

  const GF3D2<CCTK_REAL> dtalp_(layout, dtalp);
  const GF3D2<CCTK_REAL> dt2alp_(layout, dt2alp);

  const GF3D2<CCTK_REAL> dtbetax_(layout, dtbetax);
  const GF3D2<CCTK_REAL> dtbetay_(layout, dtbetay);
  const GF3D2<CCTK_REAL> dtbetaz_(layout, dtbetaz);

  const GF3D2<CCTK_REAL> dt2betax_(layout, dt2betax);
  const GF3D2<CCTK_REAL> dt2betay_(layout, dt2betay);
  const GF3D2<CCTK_REAL> dt2betaz_(layout, dt2betaz);

  const GF3D2<CCTK_REAL> dtkxx_(layout, dtkxx);
  const GF3D2<CCTK_REAL> dtkxy_(layout, dtkxy);
  const GF3D2<CCTK_REAL> dtkxz_(layout, dtkxz);
  const GF3D2<CCTK_REAL> dtkyy_(layout, dtkyy);
  const GF3D2<CCTK_REAL> dtkyz_(layout, dtkyz);
  const GF3D2<CCTK_REAL> dtkzz_(layout, dtkzz);

  auto id_lambda = [&](const PointDesc &p) {
    alp_(p.I) = 1;

    betax_(p.I) = 0.0;
    betay_(p.I) = 0.0;
    betaz_(p.I) = 0.0;

    gxx_(p.I) = 1.0;
    gxy_(p.I) = 0.0;
    gxz_(p.I) = 0.0;
    gyy_(p.I) = 1.0;
    gyz_(p.I) = 0.0;
    gzz_(p.I) = 1.0;

    kxx_(p.I) = 0.0;
    kxy_(p.I) = 0.0;
    kxz_(p.I) = 0.0;
    kyy_(p.I) = 0.0;
    kyz_(p.I) = 0.0;
    kzz_(p.I) = 0.0;

    dtalp_(p.I) = 0.0;

    dtbetax_(p.I) = 0.0;
    dtbetay_(p.I) = 0.0;
    dtbetaz_(p.I) = 0.0;

    dtkxx_(p.I) = 0.0;
    dtkxy_(p.I) = 0.0;
    dtkxz_(p.I) = 0.0;
    dtkyy_(p.I) = 0.0;
    dtkyz_(p.I) = 0.0;
    dtkzz_(p.I) = 0.0;

    dt2alp_(p.I) = 0.0;

    dt2betax_(p.I) = 0.0;
    dt2betay_(p.I) = 0.0;
    dt2betaz_(p.I) = 0.0;
  };

  loop_all<0, 0, 0>(cctkGH, id_lambda);
}
