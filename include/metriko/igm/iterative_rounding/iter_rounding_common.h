//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_ITER_ROUNDING_COMMON_H
#define METRIKO_ITER_ROUNDING_COMMON_H
#include "metriko/common/typedef.h"
#include "metriko/solver/matrix_ops.h"
#include "injective_barrier.h"

namespace metriko {
    // diagnostics for the iterative rounding: how many rounding steps run, and how much of the parameterization
    // time they spend rebuilding the constant blocks of the constraint jacobian
    inline int&    rounding_steps()    { static thread_local int v = 0; return v; }
    // diagnostics for the LM factorization: how often the symbolic analysis runs (it must not, when the pattern
    // is stable) and how the time splits between the copy, the analysis and the numeric factorization
    inline int&    lm_analyze_calls() { static thread_local int v = 0; return v; }
    inline double& lm_copy_secs()     { static thread_local double v = 0; return v; }
    inline double& lm_analyze_secs()  { static thread_local double v = 0; return v; }
    inline double& lm_numeric_secs()  { static thread_local double v = 0; return v; }
    inline double& rounding_init_secs() { static thread_local double v = 0; return v; }
    class Integration {
    public:
        const SprsD &G2;
        const VecXd &F2;
        int xSize = 0;
        int ESize = 0;
        const int N;
        const int n;
        const double length; // grid scale
        const double wInteg;
        const double wConst; // a factor so that the fixed value do not move
        const double wClose;
        const double wBarrier;
        const bool localInjectivity;
        InjectiveBarrier *barrier;
        SprsD gInteg;
        SprsD gClose;
        SprsD gConst;
        SprsD X2F; // convert x to the field with parameter length

        Integration(
            const SprsD &intrinsicGrad_,
            const VecXd &intrinsicField_,
            const int N_,
            const int n_,
            const bool locInj_,
            const double length_,
            const double wInteg_,
            const double wConst_,
            const double wClose_,
            const double wBarrier_,
            InjectiveBarrier *barrier_
        ) : G2(intrinsicGrad_),
            F2(intrinsicField_),
            N(N_),
            n(n_),
            localInjectivity(locInj_),
            length(length_),
            wInteg(wInteg_),
            wConst(wConst_),
            wClose(wClose_),
            wBarrier(wBarrier_),
            barrier(barrier_)
        {
        }

        static void pre_iteration(const VecXd &prevx) { }

        void compute_jacobian(
            const VecXd &field,
            const VecXd& fInteg,
            const VecXd& fClose,
            const VecXd& fConst,
            const bool updateJ,
            VecXd &E,
            SprsD &J,
            const bool loop
        ) {
            if (localInjectivity) {
                SprsD gBarrier;
                VecXd fBarrier;
                barrier->update(field, updateJ, fBarrier, gBarrier);
                E.resize(fInteg.size() + fClose.size() + fConst.size() + fBarrier.size());
                E << fInteg * wInteg, fClose * wClose, fConst * wConst, fBarrier * wBarrier;
                if (updateJ) {
                    // the naive branch scales the barrier by wInteg here while the residual above uses
                    // wBarrier, so this J is not the exact derivative of E. it is load-bearing: the barrier
                    // then reaches the normal equations with weight wInteg * wBarrier = 1 in the gradient
                    // but wInteg^2 = 1e8 in the gauss-newton hessian, acting as curvature rather than as a
                    // driving force. using wBarrier here instead leaves it 1e8 times weaker and nefertiti
                    // at gridscale 0.02 hits the iteration cap with 10 inverted faces. the only consistent
                    // variant that still converges raises wBarrier to wInteg, at 24 instead of 13 iterations
                    gBarrier = loop?
                        gBarrier * barrier->gen_image_filed(field) * X2F * wBarrier :
                        gBarrier * barrier->gen_image_filed(field) * gClose * wInteg; // maybe wBarrier??
                    MatXi I(4, 1);
                    I << 0, 1, 2, 3;
                    sparse_block(I, {&gInteg, &gClose, &gConst, &gBarrier}, J);
                }
                return;
            }

            E.resize(fInteg.size() + fClose.size() + fConst.size());
            E << fInteg * wInteg, fClose * wClose, fConst * wConst;
            if (updateJ) {
                MatXi I(3, 1);
                I << 0, 1, 2;
                sparse_block(I, {&gInteg, &gClose, &gConst}, J);
            }
        }
    };
}
#endif
