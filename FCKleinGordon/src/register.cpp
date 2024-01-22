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

  /*
   * Here we register the evolved variables.
   */
  const CCTK_INT state_idx = CCTK_GroupIndex("FCKleinGordon::state");
  const CCTK_INT rhs_group_idx = CCTK_GroupIndex("FCKleinGordon::rhs");

  ierr += MoLRegisterEvolvedGroup(state_idx, rhs_group_idx);

  if (ierr != 0)
    CCTK_WARN(CCTK_WARN_ABORT, "Error registering variables within MoL. Aborting.");
}
