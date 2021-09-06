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
 *  along with MinkowskiX.  If not, see <https://www.gnu.org/licenses/>.
 *
 * MinkowskiX.hpp
 * Common includes and prototypes for the MinkowskiX thorn
 */

#ifndef MINKOWSKIX_HPP
#define MINKOWSKIX_HPP

/***********************************
 * CarpetX includes part 1         *
 * This include order is important *
 ***********************************/
#include <fixmath.hxx>

/*******************
 * Cactus includes *
 *******************/
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

/***********************************
 * CarpetX includes part 2         *
 ***********************************/
#include <loop.hxx>
#include <vect.hxx>

namespace MinkowskiX {

extern "C" void MinkowskiX_Initial(CCTK_ARGUMENTS);
extern "C" int MinkowskiX_Startup(void);

} // namespace MinkowskiX

#endif // MINKOWSKIX_HPP
