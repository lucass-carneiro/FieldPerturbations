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
 *  Derivatives.h
 *  Define macros for finite differece derivatives
 */

#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#define DECLARE_DERIVATIVE_FACTORS                                             \
  const CCTK_REAL dx12 = 12.0 * CCTK_DELTA_SPACE(0),                           \
                  dy12 = 12.0 * CCTK_DELTA_SPACE(1),                           \
                  dz12 = 12.0 * CCTK_DELTA_SPACE(2),                           \
                  dxsq12 = 12.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(0),   \
                  dysq12 = 12.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(1),   \
                  dzsq12 = 12.0 * CCTK_DELTA_SPACE(2) * CCTK_DELTA_SPACE(2),   \
                  dxdy144 = 144.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(1), \
                  dxdz144 = 144.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(2), \
                  dydz144 = 144.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(2)

/* Shorthand for Cactus's index flattening */
#define I(i_, j_, k_) CCTK_GFINDEX3D(cctkGH, i_, j_, k_)

/**************************
 * FD order (accuracy): 4 *
 * Derivative Order: 2    *
 **************************/

#define D4xx(f)                                                                \
  ((-f[I(i + 2, j, k)] + 16 * f[I(i + 1, j, k)] - 30 * f[I(i, j, k)] +         \
    16 * f[I(i - 1, j, k)] - f[I(i - 2, j, k)]) /                              \
   dxsq12)

#define D4yy(f)                                                                \
  ((-f[I(i, j + 2, k)] + 16 * f[I(i, j + 1, k)] - 30 * f[I(i, j, k)] +         \
    16 * f[I(i, j - 1, k)] - f[I(i, j - 2, k)]) /                              \
   dysq12)

#define D4zz(f)                                                                \
  ((-f[I(i, j, k + 2)] + 16 * f[I(i, j, k + 1)] - 30 * f[I(i, j, k)] +         \
    16 * f[I(i, j, k - 1)] - f[I(i, j, k - 2)]) /                              \
   dzsq12)

#define D4xy(f)                                                                \
  ((-f[I(i - 2, j + 2, k)] + 8 * f[I(i - 1, j + 2, k)] -                       \
    8 * f[I(i + 1, j + 2, k)] + f[I(i + 2, j + 2, k)] +                        \
    8 * f[I(i - 2, j + 1, k)] - 64 * f[I(i - 1, j + 1, k)] +                   \
    64 * f[I(i + 1, j + 1, k)] - 8 * f[I(i + 2, j + 1, k)] -                   \
    8 * f[I(i - 2, j - 1, k)] + 64 * f[I(i - 1, j - 1, k)] -                   \
    64 * f[I(i + 1, j - 1, k)] + 8 * f[I(i + 2, j - 1, k)] +                   \
    f[I(i - 2, j - 2, k)] - 8 * f[I(i - 1, j - 2, k)] +                        \
    8 * f[I(i + 1, j - 2, k)] - f[I(i + 2, j - 2, k)]) /                       \
   dxdy144)

#define D4xz(f)                                                                \
  ((-f[I(i - 2, j, k + 2)] + 8 * f[I(i - 1, j, k + 2)] -                       \
    8 * f[I(i + 1, j, k + 2)] + f[I(i + 2, j, k + 2)] +                        \
    8 * f[I(i - 2, j, k + 1)] - 64 * f[I(i - 1, j, k + 1)] +                   \
    64 * f[I(i + 1, j, k + 1)] - 8 * f[I(i + 2, j, k + 1)] -                   \
    8 * f[I(i - 2, j, k - 1)] + 64 * f[I(i - 1, j, k - 1)] -                   \
    64 * f[I(i + 1, j, k - 1)] + 8 * f[I(i + 2, j, k - 1)] +                   \
    f[I(i - 2, j, k - 2)] - 8 * f[I(i - 1, j, k - 2)] +                        \
    8 * f[I(i + 1, j, k - 2)] - f[I(i + 2, j, k - 2)]) /                       \
   dxdz144)

#define D4yz(f)                                                                \
  ((-f[I(i, j - 2, k + 2)] + 8 * f[I(i, j - 1, k + 2)] -                       \
    8 * f[I(i, j + 1, k + 2)] + f[I(i, j + 2, k + 2)] +                        \
    8 * f[I(i, j - 2, k + 1)] - 64 * f[I(i, j - 1, k + 1)] +                   \
    64 * f[I(i, j + 1, k + 1)] - 8 * f[I(i, j + 2, k + 1)] -                   \
    8 * f[I(i, j - 2, k - 1)] + 64 * f[I(i, j - 1, k - 1)] -                   \
    64 * f[I(i, j + 1, k - 1)] + 8 * f[I(i, j + 2, k - 1)] +                   \
    f[I(i, j - 2, k - 2)] - 8 * f[I(i, j - 1, k - 2)] +                        \
    8 * f[I(i, j + 1, k - 2)] - f[I(i, j + 2, k - 2)]) /                       \
   dydz144)

/**************************
 * FD order (accuracy): 4 *
 * Derivative Order: 1    *
 **************************/

#define D4x(f)                                                                 \
  ((-f[I(i + 2, j, k)] + 8 * f[I(i + 1, j, k)] - 8 * f[I(i - 1, j, k)] +       \
    f[I(i - 2, j, k)]) /                                                       \
   dx12)

#define D4y(f)                                                                 \
  ((-f[I(i, j + 2, k)] + 8 * f[I(i, j + 1, k)] - 8 * f[I(i, j - 1, k)] +       \
    f[I(i, j - 2, k)]) /                                                       \
   dy12)

#define D4z(f)                                                                 \
  ((-f[I(i, j, k + 2)] + 8 * f[I(i, j, k + 1)] - 8 * f[I(i, j, k - 1)] +       \
    f[I(i, j, k - 2)]) /                                                       \
   dz12)

#endif /* DERIVATIVES_H */
