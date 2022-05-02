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
 * KleinGordon_ZeroError(CCTK_ARGUMENTS)      *
 *                                            *
 * This function zeros the Error variables in *
 * order to prevent sporious NaNs.            *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_ZeroError(CCTK_ARGUMENTS);

/**********************************************
 * KleinGordon_ZeroEnDen(CCTK_ARGUMENTS)      *
 *                                            *
 * This function zeros the rho_E variable     *
 * in order to prevent sporious NaNs.         *
 *                                            *
 * Input: CCTK_ARGUMENTS (the grid functions  *
 * from interface.ccl                         *
 *                                            *
 * Output: Nothing                            *
 **********************************************/
void KleinGordon_ZeroEnDen(CCTK_ARGUMENTS);

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

/***********************************************
 * KleinGordon_CalcTmunu_4(CCTK_ARGUMENTS)     *
 *                                             *
 * This function computes the energy momentum  *
 * tensor of the scalar field using 4-th order *
 * accurate finite differences.                *
 * of the ADM scalar wave equation.            *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_CalcTmunu_4(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_CalcTmunu_6(CCTK_ARGUMENTS)     *
 *                                             *
 * This function computes the energy momentum  *
 * tensor of the scalar field using 6-th order *
 * accurate finite differences.                *
 * of the ADM scalar wave equation.            *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_CalcTmunu_6(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_CalcTmunu_8(CCTK_ARGUMENTS)     *
 *                                             *
 * This function computes the energy momentum  *
 * tensor of the scalar field using 8-th order *
 * accurate finite differences.                *
 * of the ADM scalar wave equation.            *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_CalcTmunu_8(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_CalcEnDen_4(CCTK_ARGUMENTS)     *
 *                                             *
 * This function computes the energy density   *
 * of the scalar field using 4-th order        *
 * accurate finite differences.                *
 * of the ADM scalar wave equation.            *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_CalcEnDen_4(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_CalcEnDen_6(CCTK_ARGUMENTS)     *
 *                                             *
 * This function computes the energy density   *
 * of the scalar field using 6-th order        *
 * accurate finite differences.                *
 * of the ADM scalar wave equation.            *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_CalcEnDen_6(CCTK_ARGUMENTS);

/***********************************************
 * KleinGordon_CalcEnDen_8(CCTK_ARGUMENTS)     *
 *                                             *
 * This function computes the energy density   *
 * of the scalar field using 8-th order        *
 * accurate finite differences.                *
 * of the ADM scalar wave equation.            *
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl                          *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_CalcEnDen_8(CCTK_ARGUMENTS);

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

/***********************************************
 * KleinGordon_EnforceSymBound(CCTK_ARGUMENTS) *
 *                                             *
 * This function enforces the symmetry BC by   *
 * scheduling the RHS group to the "none" BC   *
 * type.
 *                                             *
 * Input: CCTK_ARGUMENTS (the grid functions   *
 * from interface.ccl)                         *
 *                                             *
 * Output: Nothing                             *
 ***********************************************/
void KleinGordon_EnforceSymBound(CCTK_ARGUMENTS);

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

/**
 * The general Gaussian solution of the wave equation in Minkowski
 * spacetime - cartesian version.
 *
 * @param t The time which the solution is to be evaluated. *
 * @param x The x cartesian coordinate
 * @param y The y cartesian coordinate
 * @param z The z cartesian coordinate
 * @param sigma The Gaussian width.
 * @return The general solution of the wave equation in Minkowski spacetime.
 */
CCTK_REAL cartesian_gaussian_solution(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                                      CCTK_REAL sigma);

/**
 * Time derivative of the general Gaussian solution of the wave equation in Minkowski
 * spacetime - cartesian version.
 *
 * @param t The time which the solution is to be evaluated. *
 * @param x The x cartesian coordinate
 * @param y The y cartesian coordinate
 * @param z The z cartesian coordinate
 * @param sigma The Gaussian width.
 * @return The time derivative of the general solution of the wave equation in Minkowski spacetime.
 */
CCTK_REAL cartesian_gaussian_solution_dt(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                                         CCTK_REAL sigma);

#endif /* DERIVATIVES_H */
