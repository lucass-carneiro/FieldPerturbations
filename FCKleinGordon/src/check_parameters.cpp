#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_check_parameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_check_parameters);
  DECLARE_CCTK_PARAMETERS;

  if (W * W < 1.0e-3) {
    CCTK_PARAMWARN("The gaussian width parameter (W) is too small. Increase it in "
                   "order to avoid singularities.");
  }
}
