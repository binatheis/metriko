#ifndef METRIKO_API_H
#define METRIKO_API_H
#include <chrono>
#include <cstdio>
#include <stdexcept>
#include <string>
#include "metriko/common/typedef.h"

namespace metriko {
enum class RemeshField { Smoothest, CurvatureAligned };

// the pipeline stage that reported a failure, or that is being timed
enum class RemeshStage {
    input,
    decimation,
    field,
    parameterization,
    tmesh,
    quantization,
    collapse,
    embedding,
    tutte,
    slim,
    qex
};

inline const char* stage_name(RemeshStage s) {
    switch (s) {
        case RemeshStage::input:            return "input";
        case RemeshStage::decimation:       return "decimation";
        case RemeshStage::field:            return "field";
        case RemeshStage::parameterization: return "parameterization";
        case RemeshStage::tmesh:            return "tmesh";
        case RemeshStage::quantization:     return "quantization";
        case RemeshStage::collapse:         return "collapse";
        case RemeshStage::embedding:        return "embedding";
        case RemeshStage::tutte:            return "tutte";
        case RemeshStage::slim:             return "slim";
        case RemeshStage::qex:              return "qex";
    }
    return "unknown";
}

// the algorithm failed inside a stage: what() is "<stage>: <underlying message>"
struct RemeshFailure : std::runtime_error {
    RemeshStage stage;
    RemeshFailure(RemeshStage s, const std::string& msg) : std::runtime_error(msg), stage(s) {}
};

inline RemeshFailure make_failure(RemeshStage s, const std::string& msg) {
    return RemeshFailure(s, std::string(stage_name(s)) + ": " + msg);
}

namespace detail {
inline bool& timings_enabled() {
    static thread_local bool v = false;
    return v;
}

// the stage an exception escaped from. not restored on unwinding, so the caller still sees it
inline RemeshStage& failed_stage() {
    static thread_local RemeshStage v = RemeshStage::input;
    return v;
}

// wall-clock deadline of the whole call (time_point::max() when there is no budget). checked between stages and
// inside the long loops (quantization pricing, SLIM): a single stage can exceed it by the time one iteration of
// its hot loop takes
inline std::chrono::steady_clock::time_point& deadline() {
    static thread_local std::chrono::steady_clock::time_point v = std::chrono::steady_clock::time_point::max();
    return v;
}

// CHOLMOD's supernodal numeric factorization parallelizes internally, and its floating-point order (hence the
// last bits of the solve) then depends on the scheduling. deterministic (the default) keeps it single-threaded.
inline bool& cholmod_deterministic() {
    static thread_local bool v = true;
    return v;
}

inline bool budget_exceeded() {
    return std::chrono::steady_clock::now() > deadline();
}

// per-stage wall limit: a stage that runs past it is abandoned (the attempt is retried with a coarser and
// smaller mesh). 0 disables it. the long-running loops check it, so it interrupts them within an iteration.
inline std::chrono::steady_clock::time_point& stage_start() {
    static thread_local auto v = std::chrono::steady_clock::now();
    return v;
}
inline double& stage_limit() {
    static thread_local double v = 0;
    return v;
}
inline bool stage_overrun() {
    if (stage_limit() <= 0) return false;
    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - stage_start()).count();
    return dt > stage_limit();
}

// thrown when a stage runs past its limit. the caller retries with a coarser mesh instead of reporting it as
// an algorithmic failure
struct StageTimeout : std::runtime_error {
    RemeshStage stage;
    StageTimeout(RemeshStage s, const std::string& msg) : std::runtime_error(msg), stage(s) {}
};

// marks the pipeline stages: remembers which one was running when an exception escapes it, and prints the
// wall time of each stage when timings are enabled. one instance lives for the whole pipeline.
class StageClock {
    RemeshStage cur_ = RemeshStage::input;
    bool marked_ = false;
    std::chrono::steady_clock::time_point t0_ = std::chrono::steady_clock::now();
    bool printed_ = false;

    void print(RemeshStage s, bool partial) const {
        if (!timings_enabled() || (!marked_ && !partial)) return;
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
        std::fprintf(stderr, "[metriko] stage %-16s %8.3f s%s\n", stage_name(s), dt, partial ? " (failed)" : "");
        std::fflush(stderr);
    }
public:
    StageClock() = default;
    ~StageClock() {
        if (std::uncaught_exceptions() > 0) {
            failed_stage() = cur_;
            print(cur_, true);
        } else if (!printed_) {
            print(cur_, false);
            printed_ = true;
        }
    }
    void mark(RemeshStage s) {
        if (budget_exceeded())
            throw std::runtime_error("the time budget was exceeded during stage " + std::string(stage_name(cur_)));
        if (s == cur_) return;
        print(cur_, false);
        t0_ = std::chrono::steady_clock::now();
        stage_start() = std::chrono::steady_clock::now();
        cur_ = s;
        marked_ = true;
    }
    RemeshStage stage() const { return cur_; }
    StageClock(const StageClock&) = delete;
    StageClock& operator=(const StageClock&) = delete;
};
}
using StageTimeout = detail::StageTimeout; // a stage ran past its limit: retry with a coarser mesh

struct DecimateOptions {
    bool   enabled            = true;     // on (auto when faces_threshold hit); false disables entirely
    int    faces_threshold    = 300000;   // auto-enable at or above this face count; 0 = always
    double vertices_percentage = 0;       // % of input verts to keep; 0 = auto cap (max_verts)
    int    max_verts          = 90000;    // auto-target cap when vertices_percentage == 0 (the parameterization and SLIM cost grow superlinearly with it, and dense decimations hit the t-mesh/slim edge cases more often)
    double normal_deviation   = 0;        // degrees; 0 = unconstrained
    double aspect_ratio       = 0;        // 0 = unconstrained
    double min_edge_length    = 0;        // 0 = unconstrained
    int    max_valence        = 0;        // 0 = unconstrained
    double hausdorff_error    = 0;        // 0 = unconstrained
    double seam_threshold     = -1;       // <0 = disabled
    double seam_deviation     = 0;        // 0 = unconstrained
};

struct RemeshOptions {
    double         scale         = 0.02;
    RemeshField    field         = RemeshField::Smoothest;
    bool           refine        = true;
    int            slim_max_iter = 30;   // SLIM hits the cap on hard meshes without converging; 30 keeps the cost bounded
    double         slim_rel_tol  = 1e-4;
    int            threads       = 0;     // 0: leave the OpenMP default
    bool           quiet         = false; // silence stdout / stderr of the library during the call
    // a stage failure is retried with a coarser scale (scale *= coarsen_factor): dense parametrizations are
    // the main source of degenerate configurations, and they get rarer as the grid gets coarser
    bool           retry_coarser = true;
    int            max_attempts  = 3;     // including the first one
    double         coarsen_factor = 1.6;
    double         time_budget_seconds = 0; // > 0: do not start another attempt after this much wall time
    double         stage_timeout_seconds = 0; // > 0: abandon a stage that runs longer than this and retry coarser
    bool           deterministic = true;      // false: let CHOLMOD factorize in parallel (faster, not reproducible)
    bool           report_timings = false;  // print the wall time of every stage to stderr
    DecimateOptions decimate;
};

struct RemeshResult {
    MatXd pos;      // quad vertices
    MatXi idx;      // quad corner indices, one quad per row
    MatXi val;      // one row per quad, column 0: patch label
    MatXd used_pos; // the triangle mesh actually remeshed (decimated if decimation ran, else the input)
    MatXi used_idx;
    double used_scale = 0; // the scale the returned mesh was produced with (differs from the option on a retry)
    int    holes = 0;      // q-face rings the extractor had to drop; the output is closed regardless
    MatXi  fill_idx;       // triangles that close the loops left by those rings (quad-dominant output)
};

// throws std::invalid_argument on a bad input mesh or options, RemeshFailure when the algorithm fails
RemeshResult remesh(const MatXd& V, const MatXi& F, const RemeshOptions& opts = {});

// decimate a closed triangle mesh with pmp::decimate; throws std::invalid_argument / std::runtime_error
std::pair<MatXd, MatXi> decimate_mesh(const MatXd& V, const MatXi& F, const DecimateOptions& o);
}
#endif
