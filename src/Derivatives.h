/*
 *  ADMScalarWave - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of ADMScalarWave.
 *
 *  ADMScalarWave is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADMScalarWave is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Derivatives.h
 * Define macros for finite differece derivatives
 */

#ifndef DERIVATIVES_H
#define DERIVATIVES_H

/* Shorthand for Cactus's index flattening */
#define I(i_, j_, k_) CCTK_GFINDEX3D(cctkGH, i_, j_, k_)

/**************************
 * FD order (accuracy): 4 *
 * Derivative Order: 2    *
 **************************/

#define D4xx(f)                                                                \
  (-(f[I(-2 + i, j, k)] - 16 * f[I(-1 + i, j, k)] + 30 * f[I(i, j, k)] -       \
     16 * f[I(1 + i, j, k)] + f[I(2 + i, j, k)]) /                             \
   (12 * dx * dx))
#define D4yy(f)                                                                \
  (-(f[I(i, -2 + j, k)] - 16 * f[I(i, -1 + j, k)] + 30 * f[I(i, j, k)] -       \
     16 * f[I(i, 1 + j, k)] + f[I(i, 2 + j, k)]) /                             \
   (12 * dy * dy))
#define D4zz(f)                                                                \
  (-(f[I(i, j, -2 + k)] - 16 * f[I(i, j, -1 + k)] + 30 * f[I(i, j, k)] -       \
     16 * f[I(i, j, 1 + k)] + f[I(i, j, 2 + k)]) /                             \
   (12 * dz * dz))

#define D4xy(f)                                                                \
  ((f[I(-2 + i, -2 + j, k)] - 8 * f[I(-2 + i, -1 + j, k)] +                    \
    8 * f[I(-2 + i, 1 + j, k)] - f[I(-2 + i, 2 + j, k)] -                      \
    8 * f[I(-1 + i, -2 + j, k)] + 64 * f[I(-1 + i, -1 + j, k)] -               \
    64 * f[I(-1 + i, 1 + j, k)] + 8 * f[I(-1 + i, 2 + j, k)] +                 \
    8 * f[I(1 + i, -2 + j, k)] - 64 * f[I(1 + i, -1 + j, k)] +                 \
    64 * f[I(1 + i, 1 + j, k)] - 8 * f[I(1 + i, 2 + j, k)] -                   \
    f[I(2 + i, -2 + j, k)] + 8 * f[I(2 + i, -1 + j, k)] -                      \
    8 * f[I(2 + i, 1 + j, k)] + f[I(2 + i, 2 + j, k)]) /                       \
   (144 * dx * dy))
#define D4xz(f)                                                                \
  ((f[I(-2 + i, j, -2 + k)] - 8 * f[I(-2 + i, j, -1 + k)] +                    \
    8 * f[I(-2 + i, j, 1 + k)] - f[I(-2 + i, j, 2 + k)] -                      \
    8 * f[I(-1 + i, j, -2 + k)] + 64 * f[I(-1 + i, j, -1 + k)] -               \
    64 * f[I(-1 + i, j, 1 + k)] + 8 * f[I(-1 + i, j, 2 + k)] +                 \
    8 * f[I(1 + i, j, -2 + k)] - 64 * f[I(1 + i, j, -1 + k)] +                 \
    64 * f[I(1 + i, j, 1 + k)] - 8 * f[I(1 + i, j, 2 + k)] -                   \
    f[I(2 + i, j, -2 + k)] + 8 * f[I(2 + i, j, -1 + k)] -                      \
    8 * f[I(2 + i, j, 1 + k)] + f[I(2 + i, j, 2 + k)]) /                       \
   (144 * dx * dz))
#define D4yz(f)                                                                \
  ((f[I(i, -2 + j, -2 + k)] - 8 * f[I(i, -2 + j, -1 + k)] +                    \
    8 * f[I(i, -2 + j, 1 + k)] - f[I(i, -2 + j, 2 + k)] -                      \
    8 * f[I(i, -1 + j, -2 + k)] + 64 * f[I(i, -1 + j, -1 + k)] -               \
    64 * f[I(i, -1 + j, 1 + k)] + 8 * f[I(i, -1 + j, 2 + k)] +                 \
    8 * f[I(i, 1 + j, -2 + k)] - 64 * f[I(i, 1 + j, -1 + k)] +                 \
    64 * f[I(i, 1 + j, 1 + k)] - 8 * f[I(i, 1 + j, 2 + k)] -                   \
    f[I(i, 2 + j, -2 + k)] + 8 * f[I(i, 2 + j, -1 + k)] -                      \
    8 * f[I(i, 2 + j, 1 + k)] + f[I(i, 2 + j, 2 + k)]) /                       \
   (144 * dy * dz))

/**************************
 * FD order (accuracy): 4 *
 * Derivative Order: 1    *
 **************************/

#define D4x(f)                                                                 \
  ((f[I(-2 + i, j, k)] - 8 * f[I(-1 + i, j, k)] + 8 * f[I(1 + i, j, k)] -      \
    f[I(2 + i, j, k)]) /                                                       \
   (12 * dx))
#define D4y(f)                                                                 \
  ((f[I(i, -2 + j, k)] - 8 * f[I(i, -1 + j, k)] + 8 * f[I(i, 1 + j, k)] -      \
    f[I(i, 2 + j, k)]) /                                                       \
   (12 * dy))
#define D4z(f)                                                                 \
  ((f[I(i, j, -2 + k)] - 8 * f[I(i, j, -1 + k)] + 8 * f[I(i, j, 1 + k)] -      \
    f[I(i, j, 2 + k)]) /                                                       \
   (12 * dz))

#endif /* DERIVATIVES_H */
