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
 *  KleinGordon.h
 *  Common includes and prototypes for the KleinGordon thorn
 */

#ifndef KLEINGORDON_H
#define KLEINGORDON_H

/*******************
 * Cactus includes *
 *******************/
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/**************************************************
 * KleinGordon_Startup(void)                      *
 *                                                *
 * This functions registers a banner in Cactus's  *
 * startup.                                       *
 *                                                *
 * Input: Nothing                                 *
 *                                                *
 * Output: 0 on success                           *
 **************************************************/
int KleinGordon_Startup(void);

/**************************************************
 * KleinGordon_CheckParameters(CCTK_ARGUMENTS)    *
 *                                                *
 * This function checks the parsed parameters to  *
 * assert that they will produce valid results.   *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl)                                 *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void KleinGordon_CheckParameters(CCTK_ARGUMENTS);

/**************************************************
 * KleinGordon_Initialize(CCTK_ARGUMENTS)         *
 *                                                *
 * This function provies the scalar fields with   *
 * initial data to begin it's time evolution.     *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void KleinGordon_Initialize(CCTK_ARGUMENTS);

/**************************************************
 * KleinGordon_MoLRegister(CCTK_ARGUMENTS)        *
 *                                                *
 * This function registers variables within the   *
 * MoL thorn. Specifically it tells it which      *
 * variables are to be evolved in time and        *
 * which variables contain the RHS to be computed *
 * during the evolution                           *
 *                                                *
 * Input: CCTK_ARGUMENTS (the grid functions from *
 * interface.ccl                                  *
 *                                                *
 * Output: Nothing                                *
 **************************************************/
void KleinGordon_MoLRegister(CCTK_ARGUMENTS);

/**********************************************
 * KleinGordon_ZeroRHS(CCTK_ARGUMENTS)        *
 *                                            *
 * This function zeros the RHS variables in   *
 * order to prevent sporious NaNs.            *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_ZeroRHS(CCTK_ARGUMENTS);

/**********************************************
 * KleinGordon_RHS_4(CCTK_ARGUMENTS)          *
 *                                            *
 * This function computes the right hand side *
 * of the ADM scalar wave equation.           *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_RHS_4(CCTK_ARGUMENTS);

/**********************************************
 * KleinGordon_RHS_6(CCTK_ARGUMENTS)          *
 *                                            *
 * This function computes the right hand side *
 * of the ADM scalar wave equation.           *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_RHS_6(CCTK_ARGUMENTS);

/**********************************************
 * KleinGordon_RHS_8(CCTK_ARGUMENTS)          *
 *                                            *
 * This function computes the right hand side *
 * of the ADM scalar wave equation.           *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_RHS_8(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordon_RHSSync(CCTK_ARGUMENTS)                          *
 *                                                              *
 * Trigger a sync of the grid functions. This is achieved by    *
 * specifying a sync in the schedule and thus this function     *
 * does not need to anything.                                   *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
void KleinGordon_RHSSync(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_OuterBoundary(CCTK_ARGUMENTS)   *
 *                                             *
 * This function updates the grid functions    *
 * at the outer boundary points with boundary  *
 * conditions.                                 *
 *                                             *
 * Users can choose radiating or reflecting    *
 * boundary conditions                         *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_RHSBoundaries(CCTK_ARGUMENTS);

/****************************************************************
 * KleinGordon_Sync(CCTK_ARGUMENTS)                             *
 *                                                              *
 * Trigger a sync of the grid functions. This is achieved by    *
 * specifying a sync in the schedule and thus this function     *
 * does not need to anything.                                   *
 *                                                              *
 * Input: CCTK_ARGUMENTS (the grid functions from interface.ccl *
 *                                                              *
 * Output: Nothing                                              *
 ****************************************************************/
void KleinGordon_Sync(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_Boundary(void)                  *
 *                                             *
 * This function is a no-op (it does nothing). *
 * After taking a time step Carpet might apply *
 * wrong boundary conditions in the refinament *
 * boundaries.                                 *
 *                                             *
 * Scheduling this function in MoL_PostStep,   *
 * postrestrict and postregrid solves this.    *
 *                                             *
 * Fuerthermore, we register all BCs as 'none' *
 * to enforces the possible symmetry BCs       *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_Boundaries(CCTK_ARGUMENTS);

/**********************************************
 * KleinGordon_Energy(CCTK_ARGUMENTS)         *
 *                                            *
 * This function computes the PDE energy for  *
 * the wave equation.                         *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_Energy(CCTK_ARGUMENTS);

/****************************************************
 * KleinGordon_Energy(CCTK_ARGUMENTS)               *
 *                                                  *
 * This function computes the error of the solution *
 * by comparing it with the analytic gaussian pulse *
 *                                                  *
 * Input: CCTK_ARGUMENTS (the grid functions        *
 * from interface.ccl                               *
 *                                                  *
 * Output: Nothing                                  *
 ****************************************************/
void KleinGordon_Error(CCTK_ARGUMENTS);

/**********************************************************************
 * exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) *
 *                                                                    *
 * Computes the exect result of a gaussian fild in a Minkowski        *
 * background for a given time and position                           *
 *                                                                    *
 * Input: The 4-D point where the gaussian shoulde be computed        *
 *                                                                    *
 * Output: The result of the gaussian at the given point              *
 **********************************************************************/
CCTK_REAL exact_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL);

/*************************************************************************
 * dt_exact_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) *
 *                                                                       *
 * Computes the time derivative exect result of a gaussian field in a    *
 * Minkowski background for a given time and position                    *
 *                                                                       *
 * Input: The 4-D point where the time derivative if the gaussian should *
 * be computed                                                           *
 *                                                                       *
 * Output: The result of the time derivative of the gaussian at the      *
 * given point                                                           *
 *************************************************************************/
CCTK_REAL dt_exact_gaussian(CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL);

#endif /* DERIVATIVES_H */
