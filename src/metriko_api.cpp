#include "metriko/api.h"
#include "metriko/lib.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#include <io.h>
#define METRIKO_DUP    _dup
#define METRIKO_DUP2   _dup2
#define METRIKO_CLOSE  _close
#define METRIKO_FILENO _fileno
#define METRIKO_NULLDEV "NUL"
#else
#include <fcntl.h>
#include <unistd.h>
#define METRIKO_DUP    dup
#define METRIKO_DUP2   dup2
#define METRIKO_CLOSE  close
#define METRIKO_FILENO fileno
#define METRIKO_NULLDEV "/dev/null"
#endif

namespace metriko {
namespace {
class StdStreamMute {
    int saved_out = -1, saved_err = -1;
    static void flush_all() {
        std::cout.flush();
        std::cerr.flush();
        std::fflush(stdout);
        std::fflush(stderr);
    }
public:
    explicit StdStreamMute(bool on) {
        if (!on) return;
        flush_all();
        std::FILE* nul = std::fopen(METRIKO_NULLDEV, "w");
        if (!nul) return;
        saved_out = METRIKO_DUP(METRIKO_FILENO(stdout));
        saved_err = METRIKO_DUP(METRIKO_FILENO(stderr));
        METRIKO_DUP2(METRIKO_FILENO(nul), METRIKO_FILENO(stdout));
        METRIKO_DUP2(METRIKO_FILENO(nul), METRIKO_FILENO(stderr));
        std::fclose(nul);
    }
    ~StdStreamMute() {
        if (saved_out < 0 && saved_err < 0) return;
        flush_all();
        if (saved_out >= 0) { METRIKO_DUP2(saved_out, METRIKO_FILENO(stdout)); METRIKO_CLOSE(saved_out); }
        if (saved_err >= 0) { METRIKO_DUP2(saved_err, METRIKO_FILENO(stderr)); METRIKO_CLOSE(saved_err); }
    }
    StdStreamMute(const StdStreamMute&) = delete;
    StdStreamMute& operator=(const StdStreamMute&) = delete;
};

class ThreadScope {
#ifdef _OPENMP
    int prev = 0;
    bool active = false;
public:
    explicit ThreadScope(int n) {
        if (n <= 0) return;
        prev = omp_get_max_threads();
        omp_set_num_threads(n);
        active = true;
    }
    ~ThreadScope() { if (active) omp_set_num_threads(prev); }
#else
public:
    explicit ThreadScope(int) {}
#endif
    ThreadScope(const ThreadScope&) = delete;
    ThreadScope& operator=(const ThreadScope&) = delete;
};

int find_root(vec<int>& parent, int i) {
    while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
}

void validate_input(const MatXd& V, const MatXi& F, const RemeshOptions& o) {
    if (!(o.scale > 0) || !std::isfinite(o.scale)) throw std::invalid_argument(std::format("scale must be a positive finite number (got {})", o.scale));
    if (o.slim_max_iter < 0) throw std::invalid_argument(std::format("slim_max_iter must be non-negative (got {})", o.slim_max_iter));
    if (!(o.slim_rel_tol >= 0) || !std::isfinite(o.slim_rel_tol)) throw std::invalid_argument(std::format("slim_rel_tol must be a non-negative finite number (got {})", o.slim_rel_tol));
    if (o.threads < 0) throw std::invalid_argument(std::format("threads must be non-negative (got {})", o.threads));
    if (V.cols() != 3) throw std::invalid_argument(std::format("V must have 3 columns (got {})", V.cols()));
    if (F.cols() != 3) throw std::invalid_argument(std::format("F must have 3 columns, i.e. a triangle mesh (got {})", F.cols()));
    if (V.rows() == 0 || F.rows() == 0) throw std::invalid_argument("the mesh is empty");
    if (!V.allFinite()) throw std::invalid_argument("V contains NaN or infinite coordinates");

    const int nV = static_cast<int>(V.rows());
    const int nF = static_cast<int>(F.rows());
    for (int i = 0; i < nF; ++i)
    for (int j = 0; j < 3; ++j) {
        const int v = F(i, j);
        if (v < 0 || v >= nV) throw std::invalid_argument(std::format("F[{}, {}] = {} is out of range [0, {})", i, j, v, nV));
    }

    std::unordered_map<std::uint64_t, int> edge_faces;
    edge_faces.reserve(static_cast<size_t>(nF) * 2);
    for (int i = 0; i < nF; ++i)
    for (int j = 0; j < 3; ++j) {
        const auto a = static_cast<std::uint64_t>(F(i, j));
        const auto b = static_cast<std::uint64_t>(F(i, (j + 1) % 3));
        if (a == b) throw std::invalid_argument(std::format("face {} is degenerate (repeated vertex {})", i, a));
        ++edge_faces[(std::min(a, b) << 32) | std::max(a, b)];
    }
    long long boundary = 0, nonmanifold = 0;
    for (const auto& [e, n] : edge_faces) {
        if (n == 1) ++boundary;
        else if (n > 2) ++nonmanifold;
    }
    if (boundary > 0) throw std::invalid_argument(std::format("the mesh has {} boundary edge(s); a closed (watertight) mesh is required", boundary));
    if (nonmanifold > 0) throw std::invalid_argument(std::format("the mesh has {} non-manifold edge(s) shared by more than 2 faces", nonmanifold));

    vec<int> parent(nV);
    std::iota(parent.begin(), parent.end(), 0);
    vec<char> used(nV, 0);
    for (int i = 0; i < nF; ++i)
    for (int j = 0; j < 3; ++j) {
        used[F(i, j)] = 1;
        parent[find_root(parent, F(i, j))] = find_root(parent, F(i, (j + 1) % 3));
    }
    int components = 0, unused = 0;
    for (int v = 0; v < nV; ++v) {
        if (!used[v]) { ++unused; continue; }
        if (find_root(parent, v) == v) ++components;
    }
    if (unused > 0) throw std::invalid_argument(std::format("{} vertex(es) are not referenced by any face", unused));
    if (components != 1) throw std::invalid_argument(std::format("the mesh has {} connected components; exactly 1 is required", components));
}
}

RemeshResult remesh(const MatXd& V, const MatXi& F, const RemeshOptions& opts) {
    validate_input(V, F, opts);
    const ThreadScope threads(opts.threads);
    const StdStreamMute mute(opts.quiet);
    detail::timings_enabled() = opts.report_timings;
    detail::failed_stage()    = RemeshStage::input;
    detail::stage_limit()     = opts.stage_timeout_seconds;
    detail::cholmod_deterministic() = opts.deterministic;
    detail::stage_start()     = std::chrono::steady_clock::now();
    const auto t_start = std::chrono::steady_clock::now();
    detail::deadline() = opts.time_budget_seconds > 0
        ? t_start + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(1.1 * opts.time_budget_seconds))
        : std::chrono::steady_clock::time_point::max();

    const bool   calign  = opts.field == RemeshField::CurvatureAligned;
    const int    attempts = std::max(1, opts.max_attempts);
    double       scale    = opts.scale;
    std::string  message;
    RemeshStage  stage = RemeshStage::input;
    const bool   dec_on = opts.decimate.enabled
                       && (opts.decimate.faces_threshold <= 0 || F.rows() >= opts.decimate.faces_threshold);

    // one attempt: decimate to this attempt's target, then run the pipeline. a retry decimates coarser as well
    // (0.6x per attempt): the decimated size drives the field / parameterization / SLIM cost quadratically, and
    // it is the strongest lever for fitting a time budget on large scans
    auto attempt_pipeline = [&](int attempt, double scale_i, MatXd& Vw, MatXi& Fw) {
        Vw = V;
        Fw = F;
        if (!dec_on) return;
        DecimateOptions dopt = opts.decimate;
        if (attempt > 0 && dopt.vertices_percentage == 0)
            dopt.max_verts = std::max(2000, int(std::lround(dopt.max_verts * std::pow(0.6, attempt))));
        detail::StageClock clock;
        clock.mark(RemeshStage::decimation);
        auto [V2, F2] = decimate_mesh(V, F, dopt);
        if (F2.rows() < 4) throw std::runtime_error(std::format("decimation left too few faces ({})", F2.rows()));
        if (F2.rows() != F.rows()) {
            std::cout << std::format("[metriko] decimated V={} F={} -> V={} F={}\n", V.rows(), F.rows(), V2.rows(), F2.rows());
            std::cout.flush();
            Vw = std::move(V2);
            Fw = std::move(F2);
            validate_input(Vw, Fw, opts);
        }
    };

    for (int attempt = 0; attempt < attempts; ++attempt) {
        const auto attempt_start = std::chrono::steady_clock::now();
        bool timed_out = false;
        MatXd Vw;
        MatXi Fw;
        try {
            attempt_pipeline(attempt, scale, Vw, Fw);
        } catch (const std::invalid_argument&) {
            throw;
        } catch (const std::exception& ex) {
            throw make_failure(RemeshStage::decimation, ex.what());
        }
        try {
            auto q = compute_quadrangulation_impl(Vw, Fw, scale, calign, opts.refine, opts.slim_max_iter, opts.slim_rel_tol);
            return RemeshResult{std::move(q.q_pos), std::move(q.q_idx), std::move(q.q_val),
                                std::move(Vw), std::move(Fw), scale, q.holes, std::move(q.q_fill)};
        } catch (const RemeshFailure& ex) {
            stage   = ex.stage;
            message = ex.what();
        } catch (const StageTimeout& ex) {
            stage   = ex.stage;
            message = ex.what();
            timed_out = true;
        } catch (const std::exception& ex) {
            stage   = detail::failed_stage();
            message = ex.what();
        }

        const double attempt_secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - attempt_start).count();
        const double spent = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
        // a coarser attempt is cheaper but repeats the field and the parameterization, so it is estimated at
        // 0.7x the cost of the attempt that just failed. a failed attempt that was cheap relative to the budget
        // gets one retry even when the estimate does not fit: it failed early, and a coarser + smaller attempt
        // is usually much cheaper than the estimate (an expensive failure has no room and is reported instead)
        const bool fits = opts.time_budget_seconds <= 0
                       || spent + 0.7 * attempt_secs <= opts.time_budget_seconds;
        // a stage that ran into its own timeout failed by construction: it gets a coarser retry while there is
        // any budget left, whatever the estimate says
        const bool early_failure = timed_out
                                || (opts.time_budget_seconds > 0 && attempt < 2 && attempt_secs <= 0.5 * opts.time_budget_seconds);
        const bool over_budget = (!fits && !early_failure)
                              || (opts.time_budget_seconds > 0 && spent >= opts.time_budget_seconds);
        if (!opts.retry_coarser || attempt + 1 >= attempts || over_budget) break;

        const double next = scale * opts.coarsen_factor;
        std::cout << std::format("[metriko] stage {} failed after {:.1f} s ({}); retrying with scale {} and a coarser decimation\n",
                                 stage_name(stage), spent, message, next);
        std::cout.flush();
        scale = next;
    }
    throw make_failure(stage, message);
}
}
