/*
 *  ADMScalarWave - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of ADMScalarWave.
 *
 *  ADMScalarWave is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADMScalarWave is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Startup.c
 * Code executed at Cactus startup.
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
int ADMScalarWave_Startup(void);

/**************************************************
 * ADMScalarWave_Startup(void)                    *
 *                                                *
 * This functions registers a banner in Cactus's  *
 * startup.                                       *
 *                                                *
 * Input: Nothing                                 *
 *                                                *
 * Output: 0 on success                           *
 **************************************************/
int ADMScalarWave_Startup(void) {
  const char *banner = "ADMScalarWave: Evolutions of a Scalar Field over an "
                       "arbitrary ADM Background";
  CCTK_RegisterBanner(banner);

  return 0;
}
