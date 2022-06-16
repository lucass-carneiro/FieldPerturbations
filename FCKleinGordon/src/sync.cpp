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
 *  sync.cpp
 *  Define no-op sync functions.
 */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#  define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

extern "C" void FCKleinGordon_sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_CHECKED(FCKleinGordon_sync);
  DECLARE_CCTK_PARAMETERS;
  return;
}