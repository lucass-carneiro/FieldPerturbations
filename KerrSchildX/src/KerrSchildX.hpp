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
 *  along with KerrSchildX.  If not, see <https://www.gnu.org/licenses/>.
 *
 * KerrSchildX.hpp
 * Common includes and prototypes for the KerrSchildX thorn
 */

#ifndef KERRSCHILDX_HPP
#define KERRSCHILDX_HPP

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
 * This include order is important *
 ***********************************/
#include <loop.hxx>

namespace KerrSchildX {

extern "C" void KerrSchildX_Initial(CCTK_ARGUMENTS);
extern "C" int KerrSchildX_Startup(void);

} // namespace KerrSchildX

#endif // KERRSCHILDX_HPP
