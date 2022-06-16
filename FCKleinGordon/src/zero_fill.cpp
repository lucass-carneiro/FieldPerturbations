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
 *  zero_fill.cpp
 *  Fill grid functions with zeros.
 */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_zero_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_zero_rhs);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_ALL(loop_zero_rhs, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    Pi_rhs[ijk] = 0.0;
    Psi_x_rhs[ijk] = 0.0;
    Psi_y_rhs[ijk] = 0.0;
    Psi_z_rhs[ijk] = 0.0;
    Phi_rhs[ijk] = 0.0;
  }
  CCTK_ENDLOOP3_ALL(loop_zero_rhs);
}

extern "C" void FCKleinGordon_zero_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_zero_error);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_ALL(loop_zero_rhs, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    Pi_error[ijk] = 0.0;
    Psi_x_error[ijk] = 0.0;
    Psi_y_error[ijk] = 0.0;
    Psi_z_error[ijk] = 0.0;
    Phi_error[ijk] = 0.0;
  }
  CCTK_ENDLOOP3_ALL(loop_zero_rhs);
}

extern "C" void FCKleinGordon_zero_multipatch_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_zero_multipatch_error);
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel
  CCTK_LOOP3_ALL(loop_zero_rhs, cctkGH, i, j, k) {
    const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    Pi_multipatch_error[ijk] = 0.0;
    Psi_x_multipatch_error[ijk] = 0.0;
    Psi_y_multipatch_error[ijk] = 0.0;
    Psi_z_multipatch_error[ijk] = 0.0;
    Phi_multipatch_error[ijk] = 0.0;
  }
  CCTK_ENDLOOP3_ALL(loop_zero_rhs);
}