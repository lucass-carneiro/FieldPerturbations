/*
 *  KleinGordon - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordon.
 *
 *  KleinGordon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordon is distributed in the hope that it will be useful,
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

#define DECLARE_DERIVATIVE_FACTORS_4                                                                 \
    const CCTK_REAL dx12 = 12.0 * CCTK_DELTA_SPACE(0), dy12 = 12.0 * CCTK_DELTA_SPACE(1),            \
                    dz12 = 12.0 * CCTK_DELTA_SPACE(2),                                               \
                    dxsq12 = 12.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(0),                       \
                    dysq12 = 12.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(1),                       \
                    dzsq12 = 12.0 * CCTK_DELTA_SPACE(2) * CCTK_DELTA_SPACE(2),                       \
                    dxdy144 = 144.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(1),                     \
                    dxdz144 = 144.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(2),                     \
                    dydz144 = 144.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(2)

#define DECLARE_FIRST_DERIVATIVE_FACTORS_4                                                           \
    const CCTK_REAL dx12 = 12.0 * CCTK_DELTA_SPACE(0), dy12 = 12.0 * CCTK_DELTA_SPACE(1),            \
                    dz12 = 12.0 * CCTK_DELTA_SPACE(2)

#define DECLARE_DERIVATIVE_FACTORS_6                                                                 \
    const CCTK_REAL dx60 = 1.0 / (60.0 * CCTK_DELTA_SPACE(0)),                                       \
                    dy60 = 1.0 / (60.0 * CCTK_DELTA_SPACE(1)),                                       \
                    dz60 = 1.0 / (60.0 * CCTK_DELTA_SPACE(2)),                                       \
                    dxsq180 = 1.0 / (180.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(0)),             \
                    dysq180 = 1.0 / (180.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(1)),             \
                    dzsq180 = 1.0 / (180.0 * CCTK_DELTA_SPACE(2) * CCTK_DELTA_SPACE(2)),             \
                    dxdy3600 = 1.0 / (3600.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(1)),           \
                    dxdz3600 = 1.0 / (3600.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(2)),           \
                    dydz3600 = 1.0 / (3600.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(2))

#define DECLARE_FIRST_DERIVATIVE_FACTORS_6                                                           \
    const CCTK_REAL dx60 = 1.0 / (60.0 * CCTK_DELTA_SPACE(0)),                                       \
                    dy60 = 1.0 / (60.0 * CCTK_DELTA_SPACE(1)),                                       \
                    dz60 = 1.0 / (60.0 * CCTK_DELTA_SPACE(2))

#define DECLARE_DERIVATIVE_FACTORS_8                                                                 \
    const CCTK_REAL dx840 = 1.0 / (840.0 * CCTK_DELTA_SPACE(0)),                                     \
                    dy840 = 1.0 / (840.0 * CCTK_DELTA_SPACE(1)),                                     \
                    dz840 = 1.0 / (840.0 * CCTK_DELTA_SPACE(2)),                                     \
                    dxsq5040 = 1.0 / (5040.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(0)),           \
                    dysq5040 = 1.0 / (5040.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(1)),           \
                    dzsq5040 = 1.0 / (5040.0 * CCTK_DELTA_SPACE(2) * CCTK_DELTA_SPACE(2)),           \
                    dxdy705600 = 1.0 / (705600.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(1)),       \
                    dxdz705600 = 1.0 / (705600.0 * CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(2)),       \
                    dydz705600 = 1.0 / (705600.0 * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(2))

#define DECLARE_FIRST_DERIVATIVE_FACTORS_8                                                           \
    const CCTK_REAL dx840 = 1.0 / (840.0 * CCTK_DELTA_SPACE(0)),                                     \
                    dy840 = 1.0 / (840.0 * CCTK_DELTA_SPACE(1)),                                     \
                    dz840 = 1.0 / (840.0 * CCTK_DELTA_SPACE(2))

/* Shorthand for Cactus's index flattening */
#define I(i_, j_, k_) CCTK_GFINDEX3D(cctkGH, i_, j_, k_)

/**************************
 * FD order (accuracy): 4 *
 * Derivative Order: 2    *
 **************************/

#define D4xx(f)                                                                                      \
    ((-f[I(i + 2, j, k)] + 16 * f[I(i + 1, j, k)] - 30 * f[I(i, j, k)] + 16 * f[I(i - 1, j, k)] -    \
      f[I(i - 2, j, k)]) /                                                                           \
     dxsq12)

#define D4yy(f)                                                                                      \
    ((-f[I(i, j + 2, k)] + 16 * f[I(i, j + 1, k)] - 30 * f[I(i, j, k)] + 16 * f[I(i, j - 1, k)] -    \
      f[I(i, j - 2, k)]) /                                                                           \
     dysq12)

#define D4zz(f)                                                                                      \
    ((-f[I(i, j, k + 2)] + 16 * f[I(i, j, k + 1)] - 30 * f[I(i, j, k)] + 16 * f[I(i, j, k - 1)] -    \
      f[I(i, j, k - 2)]) /                                                                           \
     dzsq12)

#define D4xy(f)                                                                                      \
    ((-f[I(i - 2, j + 2, k)] + 8 * f[I(i - 1, j + 2, k)] - 8 * f[I(i + 1, j + 2, k)] +               \
      f[I(i + 2, j + 2, k)] + 8 * f[I(i - 2, j + 1, k)] - 64 * f[I(i - 1, j + 1, k)] +               \
      64 * f[I(i + 1, j + 1, k)] - 8 * f[I(i + 2, j + 1, k)] - 8 * f[I(i - 2, j - 1, k)] +           \
      64 * f[I(i - 1, j - 1, k)] - 64 * f[I(i + 1, j - 1, k)] + 8 * f[I(i + 2, j - 1, k)] +          \
      f[I(i - 2, j - 2, k)] - 8 * f[I(i - 1, j - 2, k)] + 8 * f[I(i + 1, j - 2, k)] -                \
      f[I(i + 2, j - 2, k)]) /                                                                       \
     dxdy144)

#define D4xz(f)                                                                                      \
    ((-f[I(i - 2, j, k + 2)] + 8 * f[I(i - 1, j, k + 2)] - 8 * f[I(i + 1, j, k + 2)] +               \
      f[I(i + 2, j, k + 2)] + 8 * f[I(i - 2, j, k + 1)] - 64 * f[I(i - 1, j, k + 1)] +               \
      64 * f[I(i + 1, j, k + 1)] - 8 * f[I(i + 2, j, k + 1)] - 8 * f[I(i - 2, j, k - 1)] +           \
      64 * f[I(i - 1, j, k - 1)] - 64 * f[I(i + 1, j, k - 1)] + 8 * f[I(i + 2, j, k - 1)] +          \
      f[I(i - 2, j, k - 2)] - 8 * f[I(i - 1, j, k - 2)] + 8 * f[I(i + 1, j, k - 2)] -                \
      f[I(i + 2, j, k - 2)]) /                                                                       \
     dxdz144)

#define D4yz(f)                                                                                      \
    ((-f[I(i, j - 2, k + 2)] + 8 * f[I(i, j - 1, k + 2)] - 8 * f[I(i, j + 1, k + 2)] +               \
      f[I(i, j + 2, k + 2)] + 8 * f[I(i, j - 2, k + 1)] - 64 * f[I(i, j - 1, k + 1)] +               \
      64 * f[I(i, j + 1, k + 1)] - 8 * f[I(i, j + 2, k + 1)] - 8 * f[I(i, j - 2, k - 1)] +           \
      64 * f[I(i, j - 1, k - 1)] - 64 * f[I(i, j + 1, k - 1)] + 8 * f[I(i, j + 2, k - 1)] +          \
      f[I(i, j - 2, k - 2)] - 8 * f[I(i, j - 1, k - 2)] + 8 * f[I(i, j + 1, k - 2)] -                \
      f[I(i, j + 2, k - 2)]) /                                                                       \
     dydz144)

/**************************
 * FD order (accuracy): 4 *
 * Derivative Order: 1    *
 **************************/

#define D4x(f)                                                                                       \
    ((-f[I(i + 2, j, k)] + 8 * f[I(i + 1, j, k)] - 8 * f[I(i - 1, j, k)] + f[I(i - 2, j, k)]) / dx12)

#define D4y(f)                                                                                       \
    ((-f[I(i, j + 2, k)] + 8 * f[I(i, j + 1, k)] - 8 * f[I(i, j - 1, k)] + f[I(i, j - 2, k)]) / dy12)

#define D4z(f)                                                                                       \
    ((-f[I(i, j, k + 2)] + 8 * f[I(i, j, k + 1)] - 8 * f[I(i, j, k - 1)] + f[I(i, j, k - 2)]) / dz12)

/**************************
 * FD order (accuracy): 6 *
 * Derivative Order: 2    *
 **************************/

#define D6xx(f)                                                                                      \
    ((2 * f[I(-3 + i, j, k)] - 27 * f[I(-2 + i, j, k)] + 270 * f[I(-1 + i, j, k)] -                  \
      490 * f[I(i, j, k)] + 270 * f[I(1 + i, j, k)] - 27 * f[I(2 + i, j, k)] +                       \
      2 * f[I(3 + i, j, k)]) *                                                                       \
     dxsq180)

#define D6yy(f)                                                                                      \
    ((2 * f[I(i, -3 + j, k)] - 27 * f[I(i, -2 + j, k)] + 270 * f[I(i, -1 + j, k)] -                  \
      490 * f[I(i, j, k)] + 270 * f[I(i, 1 + j, k)] - 27 * f[I(i, 2 + j, k)] +                       \
      2 * f[I(i, 3 + j, k)]) *                                                                       \
     dysq180)

#define D6zz(f)                                                                                      \
    ((2 * f[I(i, j, -3 + k)] - 27 * f[I(i, j, -2 + k)] + 270 * f[I(i, j, -1 + k)] -                  \
      490 * f[I(i, j, k)] + 270 * f[I(i, j, 1 + k)] - 27 * f[I(i, j, 2 + k)] +                       \
      2 * f[I(i, j, 3 + k)]) *                                                                       \
     dzsq180)

#define D6xy(f)                                                                                      \
    ((f[I(-3 + i, -3 + j, k)] - 9 * f[I(-3 + i, -2 + j, k)] + 45 * f[I(-3 + i, -1 + j, k)] -         \
      45 * f[I(-3 + i, 1 + j, k)] + 9 * f[I(-3 + i, 2 + j, k)] - f[I(-3 + i, 3 + j, k)] -            \
      9 * f[I(-2 + i, -3 + j, k)] + 81 * f[I(-2 + i, -2 + j, k)] - 405 * f[I(-2 + i, -1 + j, k)] +   \
      405 * f[I(-2 + i, 1 + j, k)] - 81 * f[I(-2 + i, 2 + j, k)] + 9 * f[I(-2 + i, 3 + j, k)] +      \
      45 * f[I(-1 + i, -3 + j, k)] - 405 * f[I(-1 + i, -2 + j, k)] +                                 \
      2025 * f[I(-1 + i, -1 + j, k)] - 2025 * f[I(-1 + i, 1 + j, k)] +                               \
      405 * f[I(-1 + i, 2 + j, k)] - 45 * f[I(-1 + i, 3 + j, k)] - 45 * f[I(1 + i, -3 + j, k)] +     \
      405 * f[I(1 + i, -2 + j, k)] - 2025 * f[I(1 + i, -1 + j, k)] + 2025 * f[I(1 + i, 1 + j, k)] -  \
      405 * f[I(1 + i, 2 + j, k)] + 45 * f[I(1 + i, 3 + j, k)] + 9 * f[I(2 + i, -3 + j, k)] -        \
      81 * f[I(2 + i, -2 + j, k)] + 405 * f[I(2 + i, -1 + j, k)] - 405 * f[I(2 + i, 1 + j, k)] +     \
      81 * f[I(2 + i, 2 + j, k)] - 9 * f[I(2 + i, 3 + j, k)] - f[I(3 + i, -3 + j, k)] +              \
      9 * f[I(3 + i, -2 + j, k)] - 45 * f[I(3 + i, -1 + j, k)] + 45 * f[I(3 + i, 1 + j, k)] -        \
      9 * f[I(3 + i, 2 + j, k)] + f[I(3 + i, 3 + j, k)]) *                                           \
     dxdy3600)

#define D6xz(f)                                                                                      \
    ((f[I(-3 + i, j, -3 + k)] - 9 * f[I(-3 + i, j, -2 + k)] + 45 * f[I(-3 + i, j, -1 + k)] -         \
      45 * f[I(-3 + i, j, 1 + k)] + 9 * f[I(-3 + i, j, 2 + k)] - f[I(-3 + i, j, 3 + k)] -            \
      9 * f[I(-2 + i, j, -3 + k)] + 81 * f[I(-2 + i, j, -2 + k)] - 405 * f[I(-2 + i, j, -1 + k)] +   \
      405 * f[I(-2 + i, j, 1 + k)] - 81 * f[I(-2 + i, j, 2 + k)] + 9 * f[I(-2 + i, j, 3 + k)] +      \
      45 * f[I(-1 + i, j, -3 + k)] - 405 * f[I(-1 + i, j, -2 + k)] +                                 \
      2025 * f[I(-1 + i, j, -1 + k)] - 2025 * f[I(-1 + i, j, 1 + k)] +                               \
      405 * f[I(-1 + i, j, 2 + k)] - 45 * f[I(-1 + i, j, 3 + k)] - 45 * f[I(1 + i, j, -3 + k)] +     \
      405 * f[I(1 + i, j, -2 + k)] - 2025 * f[I(1 + i, j, -1 + k)] + 2025 * f[I(1 + i, j, 1 + k)] -  \
      405 * f[I(1 + i, j, 2 + k)] + 45 * f[I(1 + i, j, 3 + k)] + 9 * f[I(2 + i, j, -3 + k)] -        \
      81 * f[I(2 + i, j, -2 + k)] + 405 * f[I(2 + i, j, -1 + k)] - 405 * f[I(2 + i, j, 1 + k)] +     \
      81 * f[I(2 + i, j, 2 + k)] - 9 * f[I(2 + i, j, 3 + k)] - f[I(3 + i, j, -3 + k)] +              \
      9 * f[I(3 + i, j, -2 + k)] - 45 * f[I(3 + i, j, -1 + k)] + 45 * f[I(3 + i, j, 1 + k)] -        \
      9 * f[I(3 + i, j, 2 + k)] + f[I(3 + i, j, 3 + k)]) *                                           \
     dxdz3600)

#define D6yz(f)                                                                                      \
    ((f[I(i, -3 + j, -3 + k)] - 9 * f[I(i, -3 + j, -2 + k)] + 45 * f[I(i, -3 + j, -1 + k)] -         \
      45 * f[I(i, -3 + j, 1 + k)] + 9 * f[I(i, -3 + j, 2 + k)] - f[I(i, -3 + j, 3 + k)] -            \
      9 * f[I(i, -2 + j, -3 + k)] + 81 * f[I(i, -2 + j, -2 + k)] - 405 * f[I(i, -2 + j, -1 + k)] +   \
      405 * f[I(i, -2 + j, 1 + k)] - 81 * f[I(i, -2 + j, 2 + k)] + 9 * f[I(i, -2 + j, 3 + k)] +      \
      45 * f[I(i, -1 + j, -3 + k)] - 405 * f[I(i, -1 + j, -2 + k)] +                                 \
      2025 * f[I(i, -1 + j, -1 + k)] - 2025 * f[I(i, -1 + j, 1 + k)] +                               \
      405 * f[I(i, -1 + j, 2 + k)] - 45 * f[I(i, -1 + j, 3 + k)] - 45 * f[I(i, 1 + j, -3 + k)] +     \
      405 * f[I(i, 1 + j, -2 + k)] - 2025 * f[I(i, 1 + j, -1 + k)] + 2025 * f[I(i, 1 + j, 1 + k)] -  \
      405 * f[I(i, 1 + j, 2 + k)] + 45 * f[I(i, 1 + j, 3 + k)] + 9 * f[I(i, 2 + j, -3 + k)] -        \
      81 * f[I(i, 2 + j, -2 + k)] + 405 * f[I(i, 2 + j, -1 + k)] - 405 * f[I(i, 2 + j, 1 + k)] +     \
      81 * f[I(i, 2 + j, 2 + k)] - 9 * f[I(i, 2 + j, 3 + k)] - f[I(i, 3 + j, -3 + k)] +              \
      9 * f[I(i, 3 + j, -2 + k)] - 45 * f[I(i, 3 + j, -1 + k)] + 45 * f[I(i, 3 + j, 1 + k)] -        \
      9 * f[I(i, 3 + j, 2 + k)] + f[I(i, 3 + j, 3 + k)]) *                                           \
     dydz3600)

/**************************
 * FD order (accuracy): 6 *
 * Derivative Order: 1    *
 **************************/

#define D6x(f)                                                                                       \
    ((-f[I(-3 + i, j, k)] + 9 * f[I(-2 + i, j, k)] - 45 * f[I(-1 + i, j, k)] +                       \
      45 * f[I(1 + i, j, k)] - 9 * f[I(2 + i, j, k)] + f[I(3 + i, j, k)]) *                          \
     dx60)

#define D6y(f)                                                                                       \
    ((-f[I(i, -3 + j, k)] + 9 * f[I(i, -2 + j, k)] - 45 * f[I(i, -1 + j, k)] +                       \
      45 * f[I(i, 1 + j, k)] - 9 * f[I(i, 2 + j, k)] + f[I(i, 3 + j, k)]) *                          \
     dy60)

#define D6z(f)                                                                                       \
    ((-f[I(i, j, -3 + k)] + 9 * f[I(i, j, -2 + k)] - 45 * f[I(i, j, -1 + k)] +                       \
      45 * f[I(i, j, 1 + k)] - 9 * f[I(i, j, 2 + k)] + f[I(i, j, 3 + k)]) *                          \
     dz60)

/**************************
 * FD order (accuracy): 8 *
 * Derivative Order: 2    *
 **************************/

#define D8xx(f)                                                                                      \
    ((-9 * f[I(-4 + i, j, k)] + 128 * f[I(-3 + i, j, k)] - 1008 * f[I(-2 + i, j, k)] +               \
      8064 * f[I(-1 + i, j, k)] - 14350 * f[I(i, j, k)] + 8064 * f[I(1 + i, j, k)] -                 \
      1008 * f[I(2 + i, j, k)] + 128 * f[I(3 + i, j, k)] - 9 * f[I(4 + i, j, k)]) *                  \
     dxsq5040)

#define D8yy(f)                                                                                      \
    ((-9 * f[I(i, -4 + j, k)] + 128 * f[I(i, -3 + j, k)] - 1008 * f[I(i, -2 + j, k)] +               \
      8064 * f[I(i, -1 + j, k)] - 14350 * f[I(i, j, k)] + 8064 * f[I(i, 1 + j, k)] -                 \
      1008 * f[I(i, 2 + j, k)] + 128 * f[I(i, 3 + j, k)] - 9 * f[I(i, 4 + j, k)]) *                  \
     dysq5040)

#define D8zz(f)                                                                                      \
    ((-9 * f[I(i, j, -4 + k)] + 128 * f[I(i, j, -3 + k)] - 1008 * f[I(i, j, -2 + k)] +               \
      8064 * f[I(i, j, -1 + k)] - 14350 * f[I(i, j, k)] + 8064 * f[I(i, j, 1 + k)] -                 \
      1008 * f[I(i, j, 2 + k)] + 128 * f[I(i, j, 3 + k)] - 9 * f[I(i, j, 4 + k)]) *                  \
     dzsq5040)

#define D8xy(f)                                                                                      \
    ((9 * f[I(-4 + i, -4 + j, k)] - 96 * f[I(-4 + i, -3 + j, k)] + 504 * f[I(-4 + i, -2 + j, k)] -   \
      2016 * f[I(-4 + i, -1 + j, k)] + 2016 * f[I(-4 + i, 1 + j, k)] -                               \
      504 * f[I(-4 + i, 2 + j, k)] + 96 * f[I(-4 + i, 3 + j, k)] - 9 * f[I(-4 + i, 4 + j, k)] -      \
      96 * f[I(-3 + i, -4 + j, k)] + 1024 * f[I(-3 + i, -3 + j, k)] -                                \
      5376 * f[I(-3 + i, -2 + j, k)] + 21504 * f[I(-3 + i, -1 + j, k)] -                             \
      21504 * f[I(-3 + i, 1 + j, k)] + 5376 * f[I(-3 + i, 2 + j, k)] -                               \
      1024 * f[I(-3 + i, 3 + j, k)] + 96 * f[I(-3 + i, 4 + j, k)] + 504 * f[I(-2 + i, -4 + j, k)] -  \
      5376 * f[I(-2 + i, -3 + j, k)] + 28224 * f[I(-2 + i, -2 + j, k)] -                             \
      112896 * f[I(-2 + i, -1 + j, k)] + 112896 * f[I(-2 + i, 1 + j, k)] -                           \
      28224 * f[I(-2 + i, 2 + j, k)] + 5376 * f[I(-2 + i, 3 + j, k)] -                               \
      504 * f[I(-2 + i, 4 + j, k)] - 2016 * f[I(-1 + i, -4 + j, k)] +                                \
      21504 * f[I(-1 + i, -3 + j, k)] - 112896 * f[I(-1 + i, -2 + j, k)] +                           \
      451584 * f[I(-1 + i, -1 + j, k)] - 451584 * f[I(-1 + i, 1 + j, k)] +                           \
      112896 * f[I(-1 + i, 2 + j, k)] - 21504 * f[I(-1 + i, 3 + j, k)] +                             \
      2016 * f[I(-1 + i, 4 + j, k)] + 2016 * f[I(1 + i, -4 + j, k)] -                                \
      21504 * f[I(1 + i, -3 + j, k)] + 112896 * f[I(1 + i, -2 + j, k)] -                             \
      451584 * f[I(1 + i, -1 + j, k)] + 451584 * f[I(1 + i, 1 + j, k)] -                             \
      112896 * f[I(1 + i, 2 + j, k)] + 21504 * f[I(1 + i, 3 + j, k)] -                               \
      2016 * f[I(1 + i, 4 + j, k)] - 504 * f[I(2 + i, -4 + j, k)] + 5376 * f[I(2 + i, -3 + j, k)] -  \
      28224 * f[I(2 + i, -2 + j, k)] + 112896 * f[I(2 + i, -1 + j, k)] -                             \
      112896 * f[I(2 + i, 1 + j, k)] + 28224 * f[I(2 + i, 2 + j, k)] -                               \
      5376 * f[I(2 + i, 3 + j, k)] + 504 * f[I(2 + i, 4 + j, k)] + 96 * f[I(3 + i, -4 + j, k)] -     \
      1024 * f[I(3 + i, -3 + j, k)] + 5376 * f[I(3 + i, -2 + j, k)] -                                \
      21504 * f[I(3 + i, -1 + j, k)] + 21504 * f[I(3 + i, 1 + j, k)] -                               \
      5376 * f[I(3 + i, 2 + j, k)] + 1024 * f[I(3 + i, 3 + j, k)] - 96 * f[I(3 + i, 4 + j, k)] -     \
      9 * f[I(4 + i, -4 + j, k)] + 96 * f[I(4 + i, -3 + j, k)] - 504 * f[I(4 + i, -2 + j, k)] +      \
      2016 * f[I(4 + i, -1 + j, k)] - 2016 * f[I(4 + i, 1 + j, k)] + 504 * f[I(4 + i, 2 + j, k)] -   \
      96 * f[I(4 + i, 3 + j, k)] + 9 * f[I(4 + i, 4 + j, k)]) *                                      \
     dxdy705600)

#define D8xz(f)                                                                                      \
    ((9 * f[I(-4 + i, j, -4 + k)] - 96 * f[I(-4 + i, j, -3 + k)] + 504 * f[I(-4 + i, j, -2 + k)] -   \
      2016 * f[I(-4 + i, j, -1 + k)] + 2016 * f[I(-4 + i, j, 1 + k)] -                               \
      504 * f[I(-4 + i, j, 2 + k)] + 96 * f[I(-4 + i, j, 3 + k)] - 9 * f[I(-4 + i, j, 4 + k)] -      \
      96 * f[I(-3 + i, j, -4 + k)] + 1024 * f[I(-3 + i, j, -3 + k)] -                                \
      5376 * f[I(-3 + i, j, -2 + k)] + 21504 * f[I(-3 + i, j, -1 + k)] -                             \
      21504 * f[I(-3 + i, j, 1 + k)] + 5376 * f[I(-3 + i, j, 2 + k)] -                               \
      1024 * f[I(-3 + i, j, 3 + k)] + 96 * f[I(-3 + i, j, 4 + k)] + 504 * f[I(-2 + i, j, -4 + k)] -  \
      5376 * f[I(-2 + i, j, -3 + k)] + 28224 * f[I(-2 + i, j, -2 + k)] -                             \
      112896 * f[I(-2 + i, j, -1 + k)] + 112896 * f[I(-2 + i, j, 1 + k)] -                           \
      28224 * f[I(-2 + i, j, 2 + k)] + 5376 * f[I(-2 + i, j, 3 + k)] -                               \
      504 * f[I(-2 + i, j, 4 + k)] - 2016 * f[I(-1 + i, j, -4 + k)] +                                \
      21504 * f[I(-1 + i, j, -3 + k)] - 112896 * f[I(-1 + i, j, -2 + k)] +                           \
      451584 * f[I(-1 + i, j, -1 + k)] - 451584 * f[I(-1 + i, j, 1 + k)] +                           \
      112896 * f[I(-1 + i, j, 2 + k)] - 21504 * f[I(-1 + i, j, 3 + k)] +                             \
      2016 * f[I(-1 + i, j, 4 + k)] + 2016 * f[I(1 + i, j, -4 + k)] -                                \
      21504 * f[I(1 + i, j, -3 + k)] + 112896 * f[I(1 + i, j, -2 + k)] -                             \
      451584 * f[I(1 + i, j, -1 + k)] + 451584 * f[I(1 + i, j, 1 + k)] -                             \
      112896 * f[I(1 + i, j, 2 + k)] + 21504 * f[I(1 + i, j, 3 + k)] -                               \
      2016 * f[I(1 + i, j, 4 + k)] - 504 * f[I(2 + i, j, -4 + k)] + 5376 * f[I(2 + i, j, -3 + k)] -  \
      28224 * f[I(2 + i, j, -2 + k)] + 112896 * f[I(2 + i, j, -1 + k)] -                             \
      112896 * f[I(2 + i, j, 1 + k)] + 28224 * f[I(2 + i, j, 2 + k)] -                               \
      5376 * f[I(2 + i, j, 3 + k)] + 504 * f[I(2 + i, j, 4 + k)] + 96 * f[I(3 + i, j, -4 + k)] -     \
      1024 * f[I(3 + i, j, -3 + k)] + 5376 * f[I(3 + i, j, -2 + k)] -                                \
      21504 * f[I(3 + i, j, -1 + k)] + 21504 * f[I(3 + i, j, 1 + k)] -                               \
      5376 * f[I(3 + i, j, 2 + k)] + 1024 * f[I(3 + i, j, 3 + k)] - 96 * f[I(3 + i, j, 4 + k)] -     \
      9 * f[I(4 + i, j, -4 + k)] + 96 * f[I(4 + i, j, -3 + k)] - 504 * f[I(4 + i, j, -2 + k)] +      \
      2016 * f[I(4 + i, j, -1 + k)] - 2016 * f[I(4 + i, j, 1 + k)] + 504 * f[I(4 + i, j, 2 + k)] -   \
      96 * f[I(4 + i, j, 3 + k)] + 9 * f[I(4 + i, j, 4 + k)]) *                                      \
     dxdz705600)

#define D8yz(f)                                                                                      \
    ((9 * f[I(i, -4 + j, -4 + k)] - 96 * f[I(i, -4 + j, -3 + k)] + 504 * f[I(i, -4 + j, -2 + k)] -   \
      2016 * f[I(i, -4 + j, -1 + k)] + 2016 * f[I(i, -4 + j, 1 + k)] -                               \
      504 * f[I(i, -4 + j, 2 + k)] + 96 * f[I(i, -4 + j, 3 + k)] - 9 * f[I(i, -4 + j, 4 + k)] -      \
      96 * f[I(i, -3 + j, -4 + k)] + 1024 * f[I(i, -3 + j, -3 + k)] -                                \
      5376 * f[I(i, -3 + j, -2 + k)] + 21504 * f[I(i, -3 + j, -1 + k)] -                             \
      21504 * f[I(i, -3 + j, 1 + k)] + 5376 * f[I(i, -3 + j, 2 + k)] -                               \
      1024 * f[I(i, -3 + j, 3 + k)] + 96 * f[I(i, -3 + j, 4 + k)] + 504 * f[I(i, -2 + j, -4 + k)] -  \
      5376 * f[I(i, -2 + j, -3 + k)] + 28224 * f[I(i, -2 + j, -2 + k)] -                             \
      112896 * f[I(i, -2 + j, -1 + k)] + 112896 * f[I(i, -2 + j, 1 + k)] -                           \
      28224 * f[I(i, -2 + j, 2 + k)] + 5376 * f[I(i, -2 + j, 3 + k)] -                               \
      504 * f[I(i, -2 + j, 4 + k)] - 2016 * f[I(i, -1 + j, -4 + k)] +                                \
      21504 * f[I(i, -1 + j, -3 + k)] - 112896 * f[I(i, -1 + j, -2 + k)] +                           \
      451584 * f[I(i, -1 + j, -1 + k)] - 451584 * f[I(i, -1 + j, 1 + k)] +                           \
      112896 * f[I(i, -1 + j, 2 + k)] - 21504 * f[I(i, -1 + j, 3 + k)] +                             \
      2016 * f[I(i, -1 + j, 4 + k)] + 2016 * f[I(i, 1 + j, -4 + k)] -                                \
      21504 * f[I(i, 1 + j, -3 + k)] + 112896 * f[I(i, 1 + j, -2 + k)] -                             \
      451584 * f[I(i, 1 + j, -1 + k)] + 451584 * f[I(i, 1 + j, 1 + k)] -                             \
      112896 * f[I(i, 1 + j, 2 + k)] + 21504 * f[I(i, 1 + j, 3 + k)] -                               \
      2016 * f[I(i, 1 + j, 4 + k)] - 504 * f[I(i, 2 + j, -4 + k)] + 5376 * f[I(i, 2 + j, -3 + k)] -  \
      28224 * f[I(i, 2 + j, -2 + k)] + 112896 * f[I(i, 2 + j, -1 + k)] -                             \
      112896 * f[I(i, 2 + j, 1 + k)] + 28224 * f[I(i, 2 + j, 2 + k)] -                               \
      5376 * f[I(i, 2 + j, 3 + k)] + 504 * f[I(i, 2 + j, 4 + k)] + 96 * f[I(i, 3 + j, -4 + k)] -     \
      1024 * f[I(i, 3 + j, -3 + k)] + 5376 * f[I(i, 3 + j, -2 + k)] -                                \
      21504 * f[I(i, 3 + j, -1 + k)] + 21504 * f[I(i, 3 + j, 1 + k)] -                               \
      5376 * f[I(i, 3 + j, 2 + k)] + 1024 * f[I(i, 3 + j, 3 + k)] - 96 * f[I(i, 3 + j, 4 + k)] -     \
      9 * f[I(i, 4 + j, -4 + k)] + 96 * f[I(i, 4 + j, -3 + k)] - 504 * f[I(i, 4 + j, -2 + k)] +      \
      2016 * f[I(i, 4 + j, -1 + k)] - 2016 * f[I(i, 4 + j, 1 + k)] + 504 * f[I(i, 4 + j, 2 + k)] -   \
      96 * f[I(i, 4 + j, 3 + k)] + 9 * f[I(i, 4 + j, 4 + k)]) *                                      \
     dydz705600)

/**************************
 * FD order (accuracy): 8 *
 * Derivative Order: 1    *
 **************************/

#define D8x(f)                                                                                       \
    ((3 * f[I(-4 + i, j, k)] - 32 * f[I(-3 + i, j, k)] + 168 * f[I(-2 + i, j, k)] -                  \
      672 * f[I(-1 + i, j, k)] + 672 * f[I(1 + i, j, k)] - 168 * f[I(2 + i, j, k)] +                 \
      32 * f[I(3 + i, j, k)] - 3 * f[I(4 + i, j, k)]) *                                              \
     dx840)

#define D8y(f)                                                                                       \
    ((3 * f[I(i, -4 + j, k)] - 32 * f[I(i, -3 + j, k)] + 168 * f[I(i, -2 + j, k)] -                  \
      672 * f[I(i, -1 + j, k)] + 672 * f[I(i, 1 + j, k)] - 168 * f[I(i, 2 + j, k)] +                 \
      32 * f[I(i, 3 + j, k)] - 3 * f[I(i, 4 + j, k)]) *                                              \
     dy840)

#define D8z(f)                                                                                       \
    ((3 * f[I(i, j, -4 + k)] - 32 * f[I(i, j, -3 + k)] + 168 * f[I(i, j, -2 + k)] -                  \
      672 * f[I(i, j, -1 + k)] + 672 * f[I(i, j, 1 + k)] - 168 * f[I(i, j, 2 + k)] +                 \
      32 * f[I(i, j, 3 + k)] - 3 * f[I(i, j, 4 + k)]) *                                              \
     dz840)

#endif /* DERIVATIVES_H */
