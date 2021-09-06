/*
 *  MinkowskiX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of MinkowskiX.
 *
 *  MinkowskiX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MinkowskiX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MinkowskiX. If not, see <https://www.gnu.org/licenses/>.
 *
 *  Startup.cpp
 *  Code executed at Cactus startup.
 */

#include "MinkowskiX.hpp"

extern "C" int MinkowskiX::MinkowskiX_Startup(void) {
  const char *banner = "MinkowskiX: Minkowski's flat spacetime";
  CCTK_RegisterBanner(banner);
  return 0;
}
