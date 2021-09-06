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
     * These are the variables that the Thorn sets but doesn't evolve.
     */
    const CCTK_INT lapse_idx = CCTK_GroupIndex("ADMBase::lapse");
    const CCTK_INT shift_idx = CCTK_GroupIndex("ADMBase::shift");
    const CCTK_INT metric_idx = CCTK_GroupIndex("ADMBase::metric");
    const CCTK_INT curv_idx = CCTK_GroupIndex("ADMBase::curv");

    ierr += MoLRegisterSaveAndRestoreGroup(lapse_idx);
    ierr += MoLRegisterSaveAndRestoreGroup(shift_idx);
    ierr += MoLRegisterSaveAndRestoreGroup(metric_idx);
    ierr += MoLRegisterSaveAndRestoreGroup(curv_idx);

    /*
     * Here we register the evolved variables, the field and it's
     * conjugate momentum
     */
    const CCTK_INT evolved_idx = CCTK_GroupIndex("KleinGordon::evolved_group");
    const CCTK_INT rhs_idx = CCTK_GroupIndex("KleinGordon::rhs_group");

    ierr += MoLRegisterEvolvedGroup(evolved_idx, rhs_idx);

    if (ierr != 0)
        CCTK_WARN(CCTK_WARN_ABORT, "Error registering variables within MoL. Aborting.");
}
