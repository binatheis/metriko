//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QEX_SANITIZATION_H
#define METRIKO_QEX_SANITIZATION_H
#include "common.h"

namespace metriko::qex {
    // truncate the coordinates of one corner to the local precision and snap them to the integer iso-lines
    inline complex sanitize_uv(const complex uv, const double max, const bool singular) {
        if (singular) return {std::round(uv.real()), std::round(uv.imag())};
        // rounding through a power of two clears the mantissa bits below the local precision
        double  delta = std::pow(2, log2(max));
        complex sign((0 < uv.real()) - (uv.real() < 0), (0 < uv.imag()) - (uv.imag() < 0));
        complex u = (uv + delta * sign) - delta * sign;
        // every coordinate within tol of an integer is that integer: a vertex near a grid point becomes the
        // point, otherwise the same point is found as an edge or face q-vertex in the neighbouring faces and
        // never pairs. a vertex near a single isoline lands exactly on it, which the tracer passes by stepping
        // around the vertex (pick_next_half)
        constexpr double tol = 1e-4;
        complex g = nearby_grid(u);
        if (std::abs(u.real() - g.real()) < tol) u.real(g.real());
        if (std::abs(u.imag() - g.imag()) < tol) u.imag(g.imag());
        return u;
    }

    inline void sanitization(
        Hmesh& mesh,
        const VecXi& matching,
        const VecXi& singular,
        const int rosyN
    ) {
        VecXc& cfn = mesh.cfn;
        VecXc heR;
        VecXc heT;
        compute_trs_matrix(mesh, matching, rosyN, heR, heT);

        auto corner_uv = [&](const Crnr& c) { return cfn(c.id); };

        // the slim orientation of every face. a snap is rejected when it would make an incident face
        // degenerate or flip it: the snap is meant to remove the numeric cracks between charts, not to
        // collapse the parametrization, and a collapsed face cannot be traced by qex (it produces dangling
        // q-ports and q-faces whose ring does not close)
        VecXd area0(mesh.nF);
        for (Face f: mesh.faces) {
            auto [c0, c1, c2] = f.crnrs();
            area0[f.id] = cross(corner_uv(c1) - corner_uv(c0), corner_uv(c2) - corner_uv(c0));
        }

        for (Vert v: mesh.verts) {
            vec<int>    cids; // the corners of v, whose value follows from the first one
            vec<complex> orig;
            for (Crnr cc: v.adjCrnrs()) {
                cids.push_back(cc.id);
                orig.push_back(cfn(cc.id));
            }

            double max = 0;
            for (Crnr cc: v.adjCrnrs()) {
                auto uv = cfn(cc.id);
                max = std::max(max, abs(uv.real()));
                max = std::max(max, abs(uv.imag()));
            }

            bool init = false;
            for (Half h: v.adjHalfs()) {
                if (!init) {
                    auto cc = h.next().crnr();
                    cfn(cc.id) = sanitize_uv(cfn(cc.id), max, singular[v.id] != 0);
                } else {
                    Crnr cp = h.twin().prev().crnr();
                    Crnr cc = h.next().crnr();
                    complex r = heR(h.twin().id);
                    complex t = heT(h.twin().id);
                    cfn(cc.id) = r * cfn(cp.id) + t;
                }
                init = true;
            }

            bool keep = true;
            for (Half h: v.adjHalfs()) {
                Face f = h.face();
                if (f.id == -1) continue; // boundary: the cut mesh is closed, but do not index a missing face
                auto [c0, c1, c2] = f.crnrs();
                double ar = cross(corner_uv(c1) - corner_uv(c0), corner_uv(c2) - corner_uv(c0));
                // a face already collapsed by slim has no sign to preserve: any non-degenerate orientation
                // is better than keeping it collapsed, and flipping it is worse than the degenerate it was
                if (ar == 0) { keep = false; break; }
                const double ref = area0[f.id] != 0 ? area0[f.id] : 1.;
                if ((ref > 0) != (ar > 0)) { keep = false; break; }
            }
            if (!keep) for (size_t k = 0; k < cids.size(); ++k) cfn(cids[k]) = orig[k];
        }
    }
}
#endif
