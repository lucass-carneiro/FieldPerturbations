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

extern "C" void FCKleinGordon_MoL_register(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_MoL_register);
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0;

  /*
   * Save and restore variables are those that a thorn depends on but
   * does not set or evolve.
   */
  const CCTK_INT lapse_group_idx = CCTK_GroupIndex("ADMBase::lapse");
  const CCTK_INT shift_group_idx = CCTK_GroupIndex("ADMBase::shift");
  const CCTK_INT metric_group_idx = CCTK_GroupIndex("ADMBase::metric");
  const CCTK_INT curv_group_idx = CCTK_GroupIndex("ADMBase::curv");

  ierr += MoLRegisterSaveAndRestoreGroup(lapse_group_idx);
  ierr += MoLRegisterSaveAndRestoreGroup(shift_group_idx);
  ierr += MoLRegisterSaveAndRestoreGroup(metric_group_idx);
  ierr += MoLRegisterSaveAndRestoreGroup(curv_group_idx);

  /**
   * Constrained variables are those which your thorn sets but does not evolve.
   */
  if (compute_error) {
    const CCTK_INT error_group_idx = CCTK_GroupIndex("FCKleinGordon::error_group");
    ierr += MoLRegisterConstrainedGroup(error_group_idx);
  }

  /*
   * Here we register the evolved variables.
   */
  const CCTK_INT evolved_group_idx = CCTK_GroupIndex("FCKleinGordon::evolved_group");
  const CCTK_INT rhs_group_idx = CCTK_GroupIndex("FCKleinGordon::rhs_group");

  ierr += MoLRegisterEvolvedGroup(evolved_group_idx, rhs_group_idx);

  if (ierr != 0)
    CCTK_WARN(CCTK_WARN_ABORT, "Error registering variables within MoL. Aborting.");
}
