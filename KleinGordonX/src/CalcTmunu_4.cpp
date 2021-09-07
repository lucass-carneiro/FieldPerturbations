/*
 *  KleinGordonX - Thorn for scalar wave evolutions in arbitrary space-times
 *  Copyright (C) 2021  Lucas Timotheo Sanches
 *
 *  This file is part of KleinGordonX.
 *
 *  KleinGordonX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  KleinGordonX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with KleinGordonX.  If not, see <https://www.gnu.org/licenses/>.
 *
 * CalcTmunu.cpp
 * Compute the energy momentum tensor of the wave equation.
 */

#include "Derivative.hpp"
#include "KleinGordonX.hpp"

#include <iostream>

namespace KleinGordonX {

using Arith::vect;
using Loop::dim;
using Loop::GF3D2;
using Loop::GF3D2layout;
using Loop::loop_int;
using Loop::PointDesc;

extern "C" void KleinGordonX_CalcTmunu_4(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_KleinGordonX_CalcTmunu_4;
    DECLARE_CCTK_PARAMETERS;

    // Grid layout -----------------------------------------------
    const vect<int, dim> indextype = {0, 0, 0};
    const GF3D2layout layout(cctkGH, indextype);

    // Field grid functions --------------------------------------
    const GF3D2<const CCTK_REAL> gf_Phi(layout, Phi);
    const GF3D2<const CCTK_REAL> gf_K_Phi(layout, K_Phi);

    // ADM variables ---------------------------------------------
    const GF3D2<const CCTK_REAL> gf_alp(layout, alp);

    const GF3D2<const CCTK_REAL> gf_betax(layout, betax);
    const GF3D2<const CCTK_REAL> gf_betay(layout, betay);
    const GF3D2<const CCTK_REAL> gf_betaz(layout, betaz);

    const GF3D2<const CCTK_REAL> gf_gxx(layout, gxx);
    const GF3D2<const CCTK_REAL> gf_gxy(layout, gxy);
    const GF3D2<const CCTK_REAL> gf_gxz(layout, gxz);
    const GF3D2<const CCTK_REAL> gf_gyy(layout, gyy);
    const GF3D2<const CCTK_REAL> gf_gyz(layout, gyz);
    const GF3D2<const CCTK_REAL> gf_gzz(layout, gzz);

    // Tmunu variables -------------------------------------------
    const GF3D2<CCTK_REAL> gf_Ttt(layout, eTtt);
    const GF3D2<CCTK_REAL> gf_Ttx(layout, eTtx);
    const GF3D2<CCTK_REAL> gf_Tty(layout, eTty);
    const GF3D2<CCTK_REAL> gf_Ttz(layout, eTtz);
    const GF3D2<CCTK_REAL> gf_Txx(layout, eTxx);
    const GF3D2<CCTK_REAL> gf_Txy(layout, eTxy);
    const GF3D2<CCTK_REAL> gf_Txz(layout, eTxz);
    const GF3D2<CCTK_REAL> gf_Tyy(layout, eTyy);
    const GF3D2<CCTK_REAL> gf_Tyz(layout, eTyz);
    const GF3D2<CCTK_REAL> gf_Tzz(layout, eTzz);

    // Derivative objects -----------------------------------------
    const Derivative<Stencil::c5> d_Phi(gf_Phi);

    auto Tmunu_lambda = [&](const PointDesc &p) {
        // Local ADM variables
        const CCTK_REAL alpL = gf_alp(p.I);

        const CCTK_REAL betaxL = gf_betax(p.I);
        const CCTK_REAL betayL = gf_betay(p.I);
        const CCTK_REAL betazL = gf_betaz(p.I);

        const CCTK_REAL hxxL = gf_gxx(p.I);
        const CCTK_REAL hxyL = gf_gxy(p.I);
        const CCTK_REAL hxzL = gf_gxz(p.I);
        const CCTK_REAL hyyL = gf_gyy(p.I);
        const CCTK_REAL hyzL = gf_gyz(p.I);
        const CCTK_REAL hzzL = gf_gzz(p.I);

        // Local grid functions
        const CCTK_REAL PhiL = gf_Phi(p.I);
        const CCTK_REAL K_PhiL = gf_K_Phi(p.I);

        // Computing the inverse 3-metric (upper)
        const CCTK_REAL hdetL = -(hxzL * hxzL * hyyL) + 2 * hxyL * hxzL * hyzL - hxxL * hyzL * hyzL -
                                hxyL * hxyL * hzzL + hxxL * hyyL * hzzL;
        const CCTK_REAL ihxxL = (-hyzL * hyzL + hyyL * hzzL) / hdetL;
        const CCTK_REAL ihxyL = (hxzL * hyzL - hxyL * hzzL) / hdetL;
        const CCTK_REAL ihxzL = (-(hxzL * hyyL) + hxyL * hyzL) / hdetL;
        const CCTK_REAL ihyyL = (-hxzL * hxzL + hxxL * hzzL) / hdetL;
        const CCTK_REAL ihyzL = (hxyL * hxzL - hxxL * hyzL) / hdetL;
        const CCTK_REAL ihzzL = (-hxyL * hxyL + hxxL * hyyL) / hdetL;

        // Computing the covariant (lower) shift vector
        const CCTK_REAL ibetaxL = hxxL * betaxL + hxyL * betayL + hxzL * betazL;
        const CCTK_REAL ibetayL = hxyL * betaxL + hyyL * betayL + hyzL * betazL;
        const CCTK_REAL ibetazL = hxzL * betaxL + hyzL * betayL + hzzL * betazL;

        /* Reconstructing the 4-metric (lower).
         * It's only necessary to compute g_tt since the other componets are already
         * computed
         */
        const CCTK_REAL gttL =
            -power<2>(alpL) + ibetaxL * betaxL + ibetayL * betayL + ibetazL * betazL;

        // inverse 4-metric (upper)
        const CCTK_REAL igttL = -1.0 / power<2>(alpL);
        const CCTK_REAL igtxL = -1.0 * igttL * betaxL;
        const CCTK_REAL igtyL = -1.0 * igttL * betayL;
        const CCTK_REAL igtzL = -1.0 * igttL * betazL;
        const CCTK_REAL igxxL = ihxxL + igttL * betaxL * betaxL;
        const CCTK_REAL igxyL = ihxyL + igttL * betaxL * betayL;
        const CCTK_REAL igxzL = ihxzL + igttL * betaxL * betazL;
        const CCTK_REAL igyyL = ihyyL + igttL * betayL * betayL;
        const CCTK_REAL igyzL = ihyzL + igttL * betayL * betazL;
        const CCTK_REAL igzzL = ihzzL + igttL * betazL * betazL;

        // Derivatives of Phi
        auto grad_Phi = d_Phi.grad(p);

        const CCTK_REAL d_x_Phi = grad_Phi[0];
        const CCTK_REAL d_y_Phi = grad_Phi[1];
        const CCTK_REAL d_z_Phi = grad_Phi[2];
        const CCTK_REAL d_t_Phi =
            (betaxL * d_x_Phi + betayL * d_y_Phi + betazL * d_z_Phi) - 2.0 * alpL * K_PhiL;

        // The scalar quantity g^{ab} \nabla_{a} \phi \nabla_{b} \phi
        const CCTK_REAL nabladot = (igttL * d_t_Phi * d_t_Phi) + 2.0 * (igtxL * d_t_Phi * d_x_Phi) +
                                   2.0 * (igtyL * d_t_Phi * d_y_Phi) +
                                   2.0 * (igtzL * d_t_Phi * d_z_Phi) + (igxxL * d_x_Phi * d_x_Phi) +
                                   2.0 * (igxyL * d_x_Phi * d_y_Phi) +
                                   2.0 * (igxzL * d_x_Phi * d_z_Phi) + (igyyL * d_y_Phi * d_y_Phi) +
                                   2.0 * (igyzL * d_y_Phi * d_z_Phi) + (igzzL * d_z_Phi * d_z_Phi);

        // T_{\mu \nu} components
        gf_Ttt(p.I) += (d_t_Phi * d_t_Phi) + 0.5 * gttL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Ttx(p.I) += (d_t_Phi * d_x_Phi) + 0.5 * betaxL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Tty(p.I) += (d_t_Phi * d_y_Phi) + 0.5 * betayL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Ttz(p.I) += (d_t_Phi * d_z_Phi) + 0.5 * betazL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Txx(p.I) += (d_x_Phi * d_x_Phi) + 0.5 * hxxL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Txy(p.I) += (d_x_Phi * d_y_Phi) + 0.5 * hxyL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Txz(p.I) += (d_x_Phi * d_z_Phi) + 0.5 * hxzL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Tyy(p.I) += (d_y_Phi * d_y_Phi) + 0.5 * hyyL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Tyz(p.I) += (d_y_Phi * d_z_Phi) + 0.5 * hyzL * (power<2>(field_mass * PhiL) - nabladot);
        gf_Tzz(p.I) += (d_z_Phi * d_z_Phi) + 0.5 * hzzL * (power<2>(field_mass * PhiL) - nabladot);
    };

    loop_int<0, 0, 0>(cctkGH, Tmunu_lambda);
}
} // namespace KleinGordonX
