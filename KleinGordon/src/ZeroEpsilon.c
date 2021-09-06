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
 * ZeroEpsilon.c
 * Zeros the energy density variable to prevent spurious NaNs.
 */

/*************************
 * This thorn's includes *
 *************************/
#include "KleinGordon.h"

void KleinGordon_ZeroEpsilon(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_LOOP3_ALL(loop_epsilon, cctkGH, i, j, k) {
        const CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        epsilon[ijk] = 0.0;
    }
    CCTK_ENDLOOP3_ALL(loop_epsilon);
}
