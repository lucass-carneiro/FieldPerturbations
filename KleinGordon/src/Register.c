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
 *  Register.c
 *  Register vriables to be evolved with the MoL thorn.
 */

/*************************
 * This thorn's includes *
 *************************/
#include "KleinGordon.h"

void KleinGordon_MoLRegister(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0;

  /*
   * The ADM variables must be set as "Save and restore" within MoL.
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
   * The energy momentum tensor and error variables shoud be registered as "Constrained"
   * within mol. Constrained variables are those which your thorn sets but does not evolve.
   */
  if (compute_Tmunu) {
    const CCTK_INT energy_scalar_group_idx = CCTK_GroupIndex("TmunuBase::stress_energy_scalar");
    const CCTK_INT energy_vector_group_idx = CCTK_GroupIndex("TmunuBase::stress_energy_vector");
    const CCTK_INT energy_tensor_group_idx = CCTK_GroupIndex("TmunuBase::stress_energy_tensor");
    ierr += MoLRegisterConstrainedGroup(energy_scalar_group_idx);
    ierr += MoLRegisterConstrainedGroup(energy_vector_group_idx);
    ierr += MoLRegisterConstrainedGroup(energy_tensor_group_idx);
  }

  if (compute_error) {
    const CCTK_INT error_group_idx = CCTK_GroupIndex("KleinGordon::error_group");
    ierr += MoLRegisterConstrainedGroup(error_group_idx);
  }

  if (compute_energy_density) {
    const CCTK_INT energy_density_group_idx = CCTK_GroupIndex("KleinGordon::energy_density_group");
    ierr += MoLRegisterConstrainedGroup(energy_density_group_idx);
  }

  /*
   * Here we register the evolved variables, the field and it's
   * conjugate momentum
   */
  const CCTK_INT evolved_group_idx = CCTK_GroupIndex("KleinGordon::evolved_group");
  const CCTK_INT rhs_group_idx = CCTK_GroupIndex("KleinGordon::rhs_group");

  ierr += MoLRegisterEvolvedGroup(evolved_group_idx, rhs_group_idx);

  if (ierr != 0)
    CCTK_WARN(CCTK_WARN_ABORT, "Error registering variables within MoL. Aborting.");
}
