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
 *  check_parameters.cpp
 *  Make sures that passed parameters are valid
 */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

#include <cassert>

extern "C" void FCKleinGordon_check_parameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_check_parameters);
  DECLARE_CCTK_PARAMETERS;

  if (sigma * sigma < 1.0e-3)
    CCTK_PARAMWARN("The gaussian parameter sigma is too small. Increase it in "
                   "order to avoid singularities.");

  if (compute_error && !CCTK_Equals(initial_data, "exact_gaussian")) {
    CCTK_PARAMWARN("Error computing was requested with an initial condition other than "
                   "\"exact_gaussian\". The error estimate is only significant when "
                   "evolving \"exact_gaussian\" data on top of a Minkowski background.");
  }

  int type = 0;
  void const *ptr = CCTK_ParameterGet("cctk_itlast", "Cactus", &type);
  assert(ptr != nullptr);
  assert(type == PARAMETER_INTEGER);
  CCTK_INT const cctk_itlast = *static_cast<CCTK_INT const *>(ptr);

  if (test_multipatch && (!CCTK_Equals(initial_data, "plane_wave") || cctk_itlast != 0)) {
    CCTK_PARAMWARN("To perform a multipatch test, use plane wave data and make sure the last "
                   "cactus iteration is 0 by adding\n\n"
                   "Cactus::terminate   = \"iteration\"\n"
                   "Cactus::cctk_itlast = 0\n\n"
                   "to the parameter file.");
  }
}
