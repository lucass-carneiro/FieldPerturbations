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
 *  Symmetries.c
 *  Register field symmetries.
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
void ADMScalarWave_Symmetries(CCTK_ARGUMENTS);
void ADMScalarWave_RHSSymmetries(CCTK_ARGUMENTS);

/**************************************************
 * ADMScalarWave_Symmetries(CCTK_ARGUMENTS)       *
 *                                                *
 * This function registers symmetries for the     *
 * evolved functions                              *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void ADMScalarWave_Symmetries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0;

  const CCTK_INT sym[3] = {1, 1, 1};

  ierr += SetCartSymVN(cctkGH, sym, "ADMScalarWave::Phi");
  ierr += SetCartSymVN(cctkGH, sym, "ADMScalarWave::K_Phi");

  if (ierr != 0)
    CCTK_WARN(CCTK_WARN_ABORT,
              "Error registering symmetries for evolved variables. Aborting.");
}

/**************************************************
 * ADMScalarWave_RHSSymmetries(CCTK_ARGUMENTS)    *
 *                                                *
 * This function registers symmetries for the     *
 * RHS functions                                  *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void ADMScalarWave_RHSSymmetries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0;

  const CCTK_INT sym[3] = {1, 1, 1};

  ierr += SetCartSymVN(cctkGH, sym, "ADMScalarWave::Phi_rhs");
  ierr += SetCartSymVN(cctkGH, sym, "ADMScalarWave::K_Phi_rhs");

  if (ierr != 0)
    CCTK_WARN(CCTK_WARN_ABORT,
              "Error registering symmetries for RHS variables. Aborting.");
}
