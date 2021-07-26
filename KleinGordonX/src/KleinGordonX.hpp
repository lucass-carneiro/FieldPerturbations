/*
 *  KleinGordonX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordonX.
 *
 *  KleinGordonX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordonX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with KleinGordonX.  If not, see <https://www.gnu.org/licenses/>.
 *
 * KleinGordonX.hpp
 * Common includes and prototypes for the KleinGordonX thorn
 */

#ifndef KLEINGORDONX_HPP
#define KLEINGORDONX_HPP

/***********************************
 * CarpetX includes part 1         *
 * This include order is important *
 ***********************************/
#include <fixmath.hxx>

/*******************
 * Cactus includes *
 *******************/
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

/***********************************
 * CarpetX includes part 2         *
 * This include order is important *
 ***********************************/
#include <loop.hxx>
#include <vect.hxx>

namespace KleinGordonX {

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Spline with compact support of radius 1 and volume 1
template <typename T> T spline(T r) {
  if (r >= 1.0)
    return 0.0;
  constexpr CCTK_REAL f =
      Loop::dim == 1 ? 1.0
                     : Loop::dim == 2 ? 24.0 / 7.0 / M_PI
                                      : Loop::dim == 3 ? 4.0 / M_PI : -1;
  const T r2 = pow(r, 2);
  return f * (r <= 0.5 ? 1 - 2 * r2 : 2 + r * (-4 + 2 * r));
}

// The potential for the spline
template <typename T> T spline_potential(T r) {
  // \Laplace u = 4 \pi \rho
  if (r >= 1.0)
    return -1 / r;
  static_assert(Loop::dim == 3, "");
  const T r2 = pow(r, 2);
  return r <= 0.5
             ? -7 / T(3) + r2 * (8 / T(3) - 8 / T(5) * r2)
             : (1 / T(15) +
                r * (-40 / T(15) +
                     r2 * (80 / T(15) + r * (-80 / T(15) + r * 24 / T(15))))) /
                   r;
}

// Time derivative
// TODO: Use dual numbers for derivative
template <typename F, typename T> auto timederiv(const F &f, T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z) - f(t - dt, x, y, z)) / dt;
  };
}

// Gradient
template <typename T> auto xderiv(T f(T t, T x, T y, T z), T dx) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x + dx, y, z) - f(t, x - dx, y, z)) / (2 * dx);
  };
}
template <typename T> auto yderiv(T f(T t, T x, T y, T z), T dy) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y + dy, z) - f(t, x, y - dy, z)) / (2 * dy);
  };
}
template <typename T> auto zderiv(T f(T t, T x, T y, T z), T dz) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z + dz) - f(t, x, y, z - dz)) / (2 * dz);
  };
}

// Central potential
CCTK_REAL central_potential(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);

/****************************************************************
 * KleinGordonX_Startup()                                       *
 *                                                              *
 * This functions registers a banner in Cactus's startup        *
 * startup.                                                     *
 *                                                              *
 * Input: Nothing                                               *
 *                                                              *
 * Output: 0 on success                                         *
 ****************************************************************/
extern "C" int KleinGordonX_Startup(void);

/****************************************************************
 * KleinGordonX_CheckParameters(CCTK_ARGUMENTS)                 *
 *                                                              *
 * This function checks the parsed parameters to                *
 * assert that they will produce valid results.                 *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from               *
 * interface.ccl                                                *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_CheckParameters(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_Initialize(CCTK_ARGUMENTS)                      *
 *                                                              *
 * This function provies the scalar fields with                 *
 * initial data to begin it's time evolution.                   *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from               *
 * interface.ccl                                                *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_Initialize(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_Initialize(CCTK_ARGUMENTS)                      *
 *                                                              *
 * Trigger a sync of the grid functions. This is achieved by    *
 * specifying a sync in the schedule and thus this function     *
 * does not need to do anything                                 *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 ************* **************************************************/
extern "C" void KleinGordonX_Sync(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_Boundaries(CCTK_ARGUMENTS)                      *
 *                                                              *
 * This function updates the grid functions at the outer        *
 * boundary points with boundary conditions.                    *
 *                                                              *
 * Users can choose which BC they wish to apply via parameters  *
 * Input: CCTK_ARGUMENTS (the grid functions                    *
 * from interface.ccl                                           *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_Boundaries(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_EstimateError(CCTK_ARGUMENTS)                   *
 *                                                              *
 * This functions computes the error measure and store it the   *
 * CarpetX regrid GF. If the error estimate is above a certain  *
 * treshold, regrid occurs.                                     *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_EstimateError(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_RHS(CCTK_ARGUMENTS)                             *
 *                                                              *
 * This function computes te RHS of the wave equation.          *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_RHS(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_RHSSync(CCTK_ARGUMENTS)                         *
 *                                                              *
 * Trigger a sync of the grid functions. This is achieved by    *
 * specifying a sync in the schedule and thus this function     *
 * does not need to anything.                                   *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_RHSSync(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_RHSBoundaries(CCTK_ARGUMENTS)                   *
 *                                                              *
 * This function updates the grid functions at the outer        *
 * boundary points with boundary conditions.                    *
 *                                                              *
 * Users can choose radiating or reflecting boundary conditions *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 **************** ***********************************************/
extern "C" void KleinGordonX_RHSBoundaries(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_RHS(CCTK_ARGUMENTS)                             *
 *                                                              *
 * This function computes the energy density of the wave eq.    *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl)*
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_Energy(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordonX_Error(CCTK_ARGUMENTS)                           *
 *                                                              *
 * This function computes the error of the numerical solution   *
 * with respect to the analytical solution                      *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl)*
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
extern "C" void KleinGordonX_Error(CCTK_ARGUMENTS);

/****************************************************************
 * gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) *
 *                                                              *
 * This function computes the gaussian initial condition for a  *
 * 4D point.                                                    *
 *                                                              *
 * Input: The 4D coordinates of the point                       *
 *                                                              *
 * Output: The gaussian at the point.                           *
 ************ ***************************************************/
CCTK_REAL exact_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL);

/****************************************************************
 * dt_exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,     *
                     CCTK_REAL z)                               *
 *                                                              *
 * Computes the time derivative exect result of a gaussian      *
 * field in a Minkowski background for a given time and         *
 * position                                                     *
 *                                                              *
 * Input: The 4-D point where the time derivative if the        *
 * gaussian should be computed                                  *
 *                                                              *
 * Output: The result of the time derivative of the gaussian at *
 * the given point                                              *
 ****************************************************************/
CCTK_REAL dt_exact_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL);

} // namespace KleinGordonX

#endif // KLEINGORDONX_HPP
