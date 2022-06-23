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
 *  boundary.cpp
 *  Implements boundary conditions
 */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_outer_boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_outer_boundaries);
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(bc_type, "zero")) {
#pragma omp parallel
    CCTK_LOOP3_INTBND(loop_zero, cctkGH, i, j, k, ni, nj, nk) {
      const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      Pi[ijk] = 0.0;
      Psi_x[ijk] = 0.0;
      Psi_y[ijk] = 0.0;
      Psi_z[ijk] = 0.0;
      Phi[ijk] = 0.0;
    }
    CCTK_ENDLOOP3_INTBND(loop_zero);
  }
}

extern "C" void FCKleinGordon_rhs_outer_boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_rhs_outer_boundaries);
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(bc_type, "NewRad")) {
    CCTK_INT ierr = 0;

    ierr += NewRad_Apply(cctkGH, Pi, Pi_rhs, 0.0, 1.0, 2.0);
    ierr += NewRad_Apply(cctkGH, Psi_x, Psi_x_rhs, 0.0, 1.0, 2.0);
    ierr += NewRad_Apply(cctkGH, Psi_y, Psi_y_rhs, 0.0, 1.0, 2.0);
    ierr += NewRad_Apply(cctkGH, Psi_z, Psi_z_rhs, 0.0, 1.0, 2.0);
    ierr += NewRad_Apply(cctkGH, Phi, Phi_rhs, 0.0, 1.0, 2.0);

    if (ierr < 0)
      CCTK_ERROR("Failed to register NewRad boundary conditions");
  }
}

extern "C" void FCKleinGordon_boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_boundaries);
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_IsFunctionAliased("Boundary_SelectGroupForBC")) {
    int ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "FCKleinGordon::evolved_group", "none");
    if (ierr)
      CCTK_ERROR("Error applaying BCs in KleinGordon::evolved_group");
  }
}
