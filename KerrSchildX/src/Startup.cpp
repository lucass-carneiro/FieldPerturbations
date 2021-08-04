/*
 *  KerrSchildX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KerrSchildX.
 *
 *  KerrSchildX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KerrSchildX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with KerrSchildX. If not, see <https://www.gnu.org/licenses/>.
 *
 *  Startup.cpp
 *  Code executed at Cactus startup.
 */

#include "KerrSchildX.hpp"

extern "C" int KerrSchildX::KerrSchildX_Startup(void) {
  const char *banner =
      "KerrSchildX: Kerr black hole initial data in Kerr-Schild coordinates";
  CCTK_RegisterBanner(banner);

  return 0;
}
