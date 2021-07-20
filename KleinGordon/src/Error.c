/*
 *  KleinGordon - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordon.
 *
 *  KleinGordon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordon is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Error.c
 * Calculate the wave equation's solution error.
 * This error measure only makes sense when evolving an "exact_gaussian"
 * pulse in a Minkowski background.
 */

/*************************
 * This thorn's includes *
 *************************/
#include "KleinGordon.h"
#include "Derivatives.h"

void KleinGordon_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Loop indexes */
  CCTK_INT i = 0, j = 0, k = 0, ijk = 0;

#pragma omp parallel for
  for (k = 0; k < cctk_lsh[2]; k++) {
    for (j = 0; j < cctk_lsh[1]; j++) {
      for (i = 0; i < cctk_lsh[0]; i++) {
        ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        Phi_err[ijk] = 0;
        K_Phi_err[ijk] = 0.0;
      }
    }
  }
}
