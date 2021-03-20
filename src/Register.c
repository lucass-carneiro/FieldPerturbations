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
 * Register.c
 * Register vriables to be evolved with the MoL thorn.
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
void ADMScalarWave_MoLRegister(CCTK_ARGUMENTS);

/**************************************************
 * ADMScalarWave_MoLRegister(CCTK_ARGUMENTS)      *
 *                                                *
 * This function registers variables within the   *
 * MoL thorn. Specifically it tells it which      *
 * variables are to be evolved in time and        *
 * which variables contain the RHS to be computed *
 * during the evolution                           *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void ADMScalarWave_MoLRegister(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;

  CCTK_INT ierr = 0;

  const int evolved_idx = CCTK_GroupIndex("ADMScalarWave::evolved_group");
  const int rhs_idx = CCTK_GroupIndex("ADMScalarWave::rhs_group");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup")) {
    ierr += MoLRegisterEvolvedGroup(evolved_idx, rhs_idx);
  } else {
    CCTK_WARN(CCTK_WARN_ABORT, "MoLRegisterEvolvedGroup not aliased !");
    ++ierr;
  }

  if (ierr != 0) {
    CCTK_WARN(CCTK_WARN_ABORT,
              "Error registering the ADM wave equation with MoL. Aborting.");
  }
}
