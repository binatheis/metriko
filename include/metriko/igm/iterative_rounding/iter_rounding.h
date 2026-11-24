//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_ITER_ROUNDING_H
#define METRIKO_ITER_ROUNDING_H
#include <algorithm>
#include "metriko/solver/levenberg_marquardt.h"
#include "metriko/solver/sparse_solver.h"
#include "iter_rounding_common.h"
#include "injective_barrier.h"
#include "iter_rounding_init.h"
#include "iter_rounding_loop.h"

namespace metriko {
class CholeskyWrapper {
    vec<int> outer;
    vec<int> inner;

    template <class M>
    bool same_pattern(const M &A) const {
        return A.isCompressed()
            && std::ssize(outer) == A.outerSize() + 1
            && std::ssize(inner) == A.nonZeros()
            && std::equal(outer.begin(), outer.end(), A.outerIndexPtr())
            && std::equal(inner.begin(), inner.end(), A.innerIndexPtr());
    }


#ifdef METRIKO_HAVE_SUITESPARSE
    cholmod_common c;
    cholmod_factor* L = nullptr;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Jt;   // J row-major: the column-compressed storage of J'

    cholmod_sparse view() {   // J' as a cholmod matrix (unknowns x residuals), unsymmetric: cholmod factorizes A A'
        cholmod_sparse A{};
        A.nrow = Jt.cols(); A.ncol = Jt.rows(); A.nzmax = Jt.nonZeros();
        A.p = Jt.outerIndexPtr(); A.i = Jt.innerIndexPtr(); A.x = Jt.valuePtr();
        A.stype = 0;
        A.itype = CHOLMOD_INT;
        A.xtype = CHOLMOD_REAL;
        A.dtype = CHOLMOD_DOUBLE;
        A.sorted = 1;
        A.packed = 1;
        return A;
    }

public:
    CholeskyWrapper()  {
        cholmod_start(&c);
        c.supernodal  = CHOLMOD_SUPERNODAL;
        c.nthreads_max = detail::cholmod_deterministic() ? 1 : 0;
        // ordering: left to cholmod's trial. forcing METIS was measured and is NOT a win here (the trial picks
        // AMD: lnz 1.467e7 vs METIS 1.398e7, numeric 94.4 s vs 87.4 s, but METIS analysis 15.3 s vs 1.5 s, so
        // the total is worse): the matrix is a 2D-like lattice where AMD is already near optimal
    }
    ~CholeskyWrapper() { if (L) cholmod_free_factor(&L, &c); cholmod_finish(&c); }
    CholeskyWrapper(const CholeskyWrapper&) = delete;
    CholeskyWrapper& operator=(const CholeskyWrapper&) = delete;

    bool factorize(const SprsD &J) {
        const auto t_copy = std::chrono::steady_clock::now();
        Jt = J;
        Jt.makeCompressed();
        cholmod_sparse A = view();
        lm_copy_secs() += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_copy).count();
        if (!L || !same_pattern(Jt)) {
            ++lm_analyze_calls();
            const auto t_ana = std::chrono::steady_clock::now();
            if (L) cholmod_free_factor(&L, &c);
            L = cholmod_analyze(&A, &c);
            if (!L) { outer.clear(); inner.clear(); return false; }
            outer.assign(Jt.outerIndexPtr(), Jt.outerIndexPtr() + Jt.outerSize() + 1);
            inner.assign(Jt.innerIndexPtr(), Jt.innerIndexPtr() + Jt.nonZeros());
            lm_analyze_secs() += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ana).count();
            if (detail::timings_enabled()) {
                static const char* names[] = {"natural", "given", "amd", "metis", "nesdis", "colamd", "ccolamd", "camd"};
                const int ord = static_cast<int>(c.method[c.selected].ordering);
                std::fprintf(stderr, "[metriko] cholmod analyze: selected method %d, ordering %s, lnz %g\n",
                             c.selected, (ord >= 0 && ord < 8) ? names[ord] : "?", c.method[c.selected].lnz);
            }
        }
        const auto t_num = std::chrono::steady_clock::now();
        cholmod_factorize(&A, L, &c);
        lm_numeric_secs() += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_num).count();
        return c.status == CHOLMOD_OK && L->minor == L->n;   // minor < n: not positive definite
    }

    bool solve(const VecXd &rhs, VecXd &x) {
        cholmod_dense b{};
        b.nrow = rhs.size(); b.ncol = 1; b.nzmax = rhs.size(); b.d = rhs.size();
        b.x = const_cast<double*>(rhs.data()); b.xtype = CHOLMOD_REAL; b.dtype = CHOLMOD_DOUBLE;
        cholmod_dense* y = cholmod_solve(CHOLMOD_A, L, &b, &c);
        if (!y) return false;
        x = Eigen::Map<VecXd>((double*)y->x, y->nrow);
        cholmod_free_dense(&y, &c);
        return true;
    }
#else
    SparseLLT<SprsD> llt;

public:
    bool factorize(const SprsD &J) {
        SprsD A = J.transpose() * J;
        if (!same_pattern(A)) {
            llt.analyzePattern(A);
            if (llt.info() != Eigen::Success) { outer.clear(); inner.clear(); return false; }
            outer.assign(A.outerIndexPtr(), A.outerIndexPtr() + A.outerSize() + 1);
            inner.assign(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros());
        }
        llt.factorize(A);
        return llt.info() == Eigen::Success;
    }

    bool solve(const VecXd &rhs, VecXd &x) const { x = llt.solve(rhs); return true; }
#endif
};

    inline bool iterative_rounding(
        const VecXi &fixedIdcs,
        const VecXd &fixedVals,
        const SprsD &weightMatrix,
        const VecXi &singularIdcs,
        const VecXi &integerIdcs,
        const double length,
        const SprsD &Cfull,
        const SprsD &G2,
        const int N,
        const int n,
        const int nF,
        const bool seamless,
        const bool roundSeams,
        const bool localInjectivity,
        const bool verbose,
        const VecXd& intrinsicField,
        const std::function<bool(const VecXd& x)> &iter_cb,
        VecXd& fullx
    ) {
        using NI = NaiveIntegration;
        using SI = SeamlessIntegration<NaiveIntegration>;

        CholeskyWrapper llt_ni, llt_si;

        DiagonalDamping<NI> dd_ni(localInjectivity ? 0.01 : 0);
        LMSolver<CholeskyWrapper, NI, DiagonalDamping<NI> > isLM;

        InjectiveBarrier barrier(N, nF, intrinsicField);

        NI ni(
            G2,
            Cfull,
            intrinsicField,
            fixedIdcs,
            fixedVals,
            weightMatrix,
            length,
            localInjectivity,
            N,
            n,
            iter_cb,
            &barrier
        );

        //----- Naive Solution -----//
        if (!seamless) {
            // Either it doesn't need loc_inj or it is already loc_inj (because post_iteration has loc_inj check)
            if (!localInjectivity || ni.post_iteration(ni.XF_Small)) { fullx = ni.x0; return true; }
        }

        // the initial solve is the most expensive stage on some meshes (about 4 s per iteration at 60k vertices;
        // 40 of 100 iterations on gatubelo). truncating it is not safe: the seamless phase consumes the
        // converged field, and an early stop cascades into a crash later in the pipeline
        isLM.init(&llt_ni, &ni, &dd_ni);
        const auto t_lm = std::chrono::steady_clock::now();
        isLM.solve(verbose);
        if (detail::timings_enabled())
            std::fprintf(stderr, "[metriko] cholmod: %d analyses, copy %.1f s, analyze %.1f s, numeric %.1f s\n",
                         lm_analyze_calls(), lm_copy_secs(), lm_analyze_secs(), lm_numeric_secs());
        if (detail::timings_enabled())
            std::fprintf(stderr, "[metriko] naive LM: %d of %d iterations, %.3f s\n", isLM.curIter, isLM.maxIter,
                         std::chrono::duration<double>(std::chrono::steady_clock::now() - t_lm).count());

        // read-only check: F2 is stored expanded as N = 4 sign-symmetry copies (n = N/2 = 2 independent ones).
        // if the optimised copies stayed in the [I; -I] pattern, the expansion is a pure reparameterisation and the
        // LM could be posed on the reduced block; if they diverged, the extra freedom is being used and it cannot
        if (detail::timings_enabled()) {
            const VecXd F2chk = (ni.UExt * isLM.x).tail(2 * N * nF);
            if (F2chk.size() >= 4 * nF) {
                const VecXd a = F2chk.head(2 * nF);
                const VecXd b = F2chk.segment(2 * nF, 2 * nF);
                const double nrm = std::max(F2chk.norm(), 1e-300);
                std::fprintf(stderr, "[metriko] F2 sign symmetry: |F2| %.6g, |lower+upper|/|F2| %.3e, |second pair+first|/|F2| %.3e\n",
                             nrm, (a + b).norm() / nrm, (a - b).norm() / nrm);
            }
        }
        // todo need to assign minix!
        fullx = (ni.UExt * isLM.x).head(ni.UFull.rows());
        if (!seamless) return true;

        //----- Seamless Solution -----//
        VecXd X_  = isLM.x.head(ni.UFull.cols());
        VecXd F2_ = (ni.UExt * isLM.x).tail(2 * N * nF); // x is changed from F2, but how is this important??
        SI si(ni, X_, F2_, integerIdcs, singularIdcs, roundSeams);

        bool success = true;
        bool rounded = false;
        DiagonalDamping<SI> dd_si(localInjectivity ? 0.01 : 0);
        LMSolver<CholeskyWrapper, SI, DiagonalDamping<SI> > irLM;

        while (!si.leftIdcs.empty()) {
            if (!si.initFixedIndices()) continue;
            rounded = true;
            dd_si.currLambda = localInjectivity ? 0.01 : 0.;
            irLM.init(&llt_si, &si, &dd_si, 100, 1e-7, 1e-7);
            irLM.solve(verbose);
            if (!si.post_checking(irLM.x)) {
                success = false;
                break;
            }
        }

        fullx = rounded ? ni.UFull * irLM.x : ni.x0;
        return success;
    }
}
#endif
