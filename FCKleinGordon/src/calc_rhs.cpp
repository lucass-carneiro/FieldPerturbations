/*
 *  FCKleinGordon - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of FCKleinGordon.
 *
 *  FCKleinGordon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FCKleinGordon is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with FCKleinGordon.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  calc_rhs.cpp
 *  Computes the RHS of the KG evolution equations.
 */

#include "derivatives.hpp"

#include <array>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_calc_rhs(CCTK_ARGUMENTS) {
  using namespace fckg;

  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_calc_rhs);
  DECLARE_CCTK_PARAMETERS;

  // Inverse metric
  auto det_gamma = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    return -(gxz[ijk] * gxz[ijk] * gyy[ijk]) + 2 * gxy[ijk] * gxz[ijk] * gyz[ijk]
           - gxx[ijk] * gyz[ijk] * gyz[ijk] - gxy[ijk] * gxy[ijk] * gzz[ijk]
           + gxx[ijk] * gyy[ijk] * gzz[ijk];
  };

  auto igxx = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return (-gyz[ijk] * gyz[ijk] + gyy[ijk] * gzz[ijk]) / det_gamma(i, j, k);
  };

  auto igxy = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return (gxz[ijk] * gyz[ijk] - gxy[ijk] * gzz[ijk]) / det_gamma(i, j, k);
  };

  auto igxz = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return (gxy[ijk] * gyz[ijk] - (gxz[ijk] * gyy[ijk])) / det_gamma(i, j, k);
  };

  auto igyy = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return (-gxz[ijk] * gxz[ijk] + gxx[ijk] * gzz[ijk]) / det_gamma(i, j, k);
  };

  auto igyz = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return (gxy[ijk] * gxz[ijk] - gxx[ijk] * gyz[ijk]) / det_gamma(i, j, k);
  };

  auto igzz = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return (-gxy[ijk] * gxy[ijk] + gxx[ijk] * gyy[ijk]) / det_gamma(i, j, k);
  };

  // Flux vectors
  auto F_x = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return alp[ijk] * det_gamma(i, j, k)
               * (igxx(i, j, k) * Psi_x[ijk] + igxy(i, j, k) * Psi_y[ijk]
                  + igxz(i, j, k) * Psi_z[ijk])
           - betax[ijk] * Pi[ijk];
  };

  auto F_y = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return alp[ijk] * det_gamma(i, j, k)
               * (igxy(i, j, k) * Psi_x[ijk] + igyy(i, j, k) * Psi_y[ijk]
                  + igyz(i, j, k) * Psi_z[ijk])
           - betay[ijk] * Pi[ijk];
  };

  auto F_z = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return alp[ijk] * det_gamma(i, j, k)
               * (igxz(i, j, k) * Psi_x[ijk] + igyz(i, j, k) * Psi_y[ijk]
                  + igzz(i, j, k) * Psi_z[ijk])
           - betaz[ijk] * Pi[ijk];
  };

  // RHS of Phi's evolution equation
  auto dt_Phi = [=](CCTK_INT i, CCTK_INT j, CCTK_INT k) {
    using std::sqrt;
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
    return betax[ijk] * Psi_x[ijk] + betay[ijk] * Psi_y[ijk] + betaz[ijk] * Psi_z[ijk]
           - alp[ijk] * Pi[ijk] / sqrt(det_gamma(i, j, k));
  };

  const CCTK_INT nx = cctk_lsh[0];
  const CCTK_INT ny = cctk_lsh[1];
  const CCTK_INT nz = cctk_lsh[2];

  const CCTK_REAL hx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL hy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL hz = CCTK_DELTA_SPACE(2);

  // Get SBP coefficients
  sbp_coefficients<derivative_direction::x> coeffs_x(cctkGH, nx);
  sbp_coefficients<derivative_direction::y> coeffs_y(cctkGH, ny);
  sbp_coefficients<derivative_direction::z> coeffs_z(cctkGH, nz);

#pragma omp parallel
  CCTK_LOOP3_INT(loop_rhs, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    std::array<std::array<CCTK_REAL, 3>, 3> jac{};
    jac[0][0] = J11[ijk];
    jac[1][0] = J21[ijk];
    jac[2][0] = J31[ijk];

    jac[0][1] = J12[ijk];
    jac[1][1] = J22[ijk];
    jac[2][1] = J32[ijk];

    jac[0][2] = J13[ijk];
    jac[1][2] = J23[ijk];
    jac[2][2] = J33[ijk];

    Pi_rhs[ijk] = alp[ijk] * det_gamma(i, j, k) * field_mass * field_mass * Phi[ijk]
                  - global_d<derivative_direction::x>(coeffs_x, coeffs_y, coeffs_z, F_x, jac, i, j,
                                                      k, hx, hy, hz)
                  - global_d<derivative_direction::y>(coeffs_x, coeffs_y, coeffs_z, F_y, jac, i, j,
                                                      k, hx, hy, hz)
                  - global_d<derivative_direction::z>(coeffs_x, coeffs_y, coeffs_z, F_z, jac, i, j,
                                                      k, hx, hy, hz);

    Psi_x_rhs[ijk] = global_d<derivative_direction::x>(coeffs_x, coeffs_y, coeffs_z, dt_Phi, jac, i,
                                                       j, k, hx, hy, hz);
    Psi_y_rhs[ijk] = global_d<derivative_direction::y>(coeffs_x, coeffs_y, coeffs_z, dt_Phi, jac, i,
                                                       j, k, hx, hy, hz);
    Psi_z_rhs[ijk] = global_d<derivative_direction::z>(coeffs_x, coeffs_y, coeffs_z, dt_Phi, jac, i,
                                                       j, k, hx, hy, hz);

    Phi_rhs[ijk] = dt_Phi(i, j, k);
  }
  CCTK_ENDLOOP3_INT(loop_rhs);
}
