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
 * CheckParameters.c
 * Check if the parameters parsed by Cactus are valid.
 */

/*******************
 * Cactus includes *
 *******************/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/**************
 * Prototypes *
 **************/
void KleinGordon_CheckParameters(CCTK_ARGUMENTS);

/**************************************************
 * KleinGordon_CheckParameters(CCTK_ARGUMENTS)  *
 *                                                *
 * This function checks the parsed parameters to  *
 * assert that they will produce valid results.   *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void KleinGordon_CheckParameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (gaussian_sigma * gaussian_sigma < 1.0e-3)
    CCTK_PARAMWARN("The gaussian parameter sigma is too small. Increase it in "
                   "order to avoid singularities.");

  switch (fd_order) {
  case 4: {
    CCTK_INFO("Using 4th order finite differencing. Make sure that you have at "
              "least 2 ghost zones");
    break;
  }

  case 6: {
    CCTK_INFO("Using 6th order finite differencing. Make sure that you have at "
              "least 3 ghost zones");
    break;
  }

  case 8: {
    CCTK_INFO("Using 8th order finite differencing. Make sure that you have at "
              "least 4 ghost zones");
    break;
  }
  }
}
