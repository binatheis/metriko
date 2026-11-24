#include <CLI/CLI.hpp>
#include <igl/readOBJ.h>
#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include "metriko/api.h"

#ifndef METRIKO_VERSION
#define METRIKO_VERSION "0.1.0"
#endif

namespace {
constexpr int EXIT_ALGORITHM = 1;
constexpr int EXIT_USAGE_IO  = 2;

bool write_quads(const std::string& path, const metriko::RemeshResult& r) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    for (int i = 0; i < r.pos.rows(); ++i) std::fprintf(f, "v %.17g %.17g %.17g\n", r.pos(i, 0), r.pos(i, 1), r.pos(i, 2));
    for (int i = 0; i < r.idx.rows(); ++i) std::fprintf(f, "f %d %d %d %d\n", r.idx(i, 0) + 1, r.idx(i, 1) + 1, r.idx(i, 2) + 1, r.idx(i, 3) + 1);
    for (int i = 0; i < r.fill_idx.rows(); ++i) std::fprintf(f, "f %d %d %d\n", r.fill_idx(i, 0) + 1, r.fill_idx(i, 1) + 1, r.fill_idx(i, 2) + 1);
    return std::fclose(f) == 0;
}

bool write_patches(const std::string& path, const metriko::RemeshResult& r) {
    std::ofstream f(path);
    if (!f) return false;
    for (int i = 0; i < r.val.rows(); ++i) f << r.val(i, 0) << '\n';
    return static_cast<bool>(f);
}
}

int run(int argc, char** argv);

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return EXIT_USAGE_IO;
    }
}

int run(int argc, char** argv) {
    CLI::App app{"metriko-cli: quad remeshing of closed triangle meshes (OBJ in, OBJ out)"};
    app.set_version_flag("--version", METRIKO_VERSION);

    std::string input, output, patches, field = "smoothest";
    bool no_refine = false, quiet = false, timings = false, no_retry = false, decimate_only = false, fast = false;
    metriko::RemeshOptions opts;
    opts.time_budget_seconds = 300; // production guard: a run stops at the next stage boundary after 5 minutes
    opts.stage_timeout_seconds = 150; // a single stage may not take more than half of it: it is retried coarser

    app.add_option("-i,--input", input, "input triangle mesh (OBJ)")->required()->check(CLI::ExistingFile);
    app.add_option("-o,--output", output, "output quad mesh (OBJ)")->required();
    app.add_option("--scale", opts.scale, "quad size relative to the bounding box")->capture_default_str()->check(CLI::PositiveNumber);
    app.add_option("--field", field, "direction field")->capture_default_str()->check(CLI::IsMember({"smoothest", "curvature_aligned"}));
    app.add_flag("--no-refine", no_refine, "disable the refinement of the extracted quad mesh");
    app.add_option("--slim-max-iter", opts.slim_max_iter, "maximum number of SLIM iterations (landing on the cap means it did not converge)")->capture_default_str()->check(CLI::NonNegativeNumber);
    app.add_option("--slim-rel-tol", opts.slim_rel_tol, "relative tolerance of the SLIM iterations")->capture_default_str()->check(CLI::NonNegativeNumber);
    app.add_option("--threads", opts.threads, "number of OpenMP threads (0: automatic)")->capture_default_str()->check(CLI::NonNegativeNumber);
    app.add_option("--patches", patches, "write the patch label of every quad to this text file");
    app.add_flag("--quiet", quiet, "silence the diagnostics of the library");
    app.add_flag("--fast", fast, "let CHOLMOD factorize in parallel: faster, but the output is no longer reproducible");
    app.add_flag("--timings", timings, "print the wall time of every pipeline stage to stderr");
    app.add_flag("--no-retry", no_retry, "do not retry with a coarser scale when a stage fails");
    app.add_option("--coarsen-factor", opts.coarsen_factor, "scale multiplier of each retry")->capture_default_str()->check(CLI::PositiveNumber);
    app.add_option("--max-attempts", opts.max_attempts, "including the first attempt")->capture_default_str()->check(CLI::PositiveNumber);
    app.add_option("--time-budget", opts.time_budget_seconds, "seconds after which no new retry is started (0: unlimited)")->capture_default_str();
    app.add_option("--stage-timeout", opts.stage_timeout_seconds, "seconds a single stage may run before it is abandoned and retried coarser (0: no limit)")->capture_default_str();
    app.add_flag("--decimate-only", decimate_only, "only decimate (with the --decimate-* options) and write the triangle mesh to the output");

    std::string decimated_out;
    bool dec_force = false, dec_off = false;
    app.add_flag("--no-decimate", dec_off, "disable decimation entirely");
    app.add_flag("--decimate", dec_force, "force decimation even below the face threshold")
        ->excludes("--no-decimate");
    auto& d = opts.decimate;
    app.add_option("--decimate-verts-pct", d.vertices_percentage, "decimate: % of input vertices to keep (0: automatic cap)")->capture_default_str();
    app.add_option("--decimate-max-verts", d.max_verts, "decimate: automatic target vertex cap")->capture_default_str();
    app.add_option("--decimate-faces-threshold", d.faces_threshold, "decimate: input face count that triggers automatic decimation")->capture_default_str();
    app.add_option("--decimate-normal-dev", d.normal_deviation, "decimate: max normal deviation in degrees (0: unconstrained)")->capture_default_str();
    app.add_option("--decimate-aspect-ratio", d.aspect_ratio, "decimate: max triangle aspect ratio (0: unconstrained)")->capture_default_str();
    app.add_option("--decimate-min-edge-len", d.min_edge_length, "decimate: minimum target edge length in mesh units (0: unconstrained)")->capture_default_str();
    app.add_option("--decimate-max-valence", d.max_valence, "decimate: maximum vertex valence (0: unconstrained)")->capture_default_str();
    app.add_option("--decimate-hausdorff", d.hausdorff_error, "decimate: maximum deviation from the original surface in mesh units (0: unconstrained)")->capture_default_str();
    app.add_option("--decimate-seam-thresh", d.seam_threshold, "decimate: texture seam detection threshold (<0: disabled; needs h:tex, not loaded for OBJ positions)")->capture_default_str();
    app.add_option("--decimate-seam-dev", d.seam_deviation, "decimate: max texture seam deviation in degrees (0: unconstrained)")->capture_default_str();
    app.add_option("--decimated-out", decimated_out, "also write the decimated triangle mesh to this OBJ file");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e) == 0 ? 0 : EXIT_USAGE_IO;
    }

    opts.field  = field == "curvature_aligned" ? metriko::RemeshField::CurvatureAligned : metriko::RemeshField::Smoothest;
    opts.refine = !no_refine;
    opts.quiet  = quiet;
    opts.report_timings = timings;
    opts.retry_coarser  = !no_retry;
    opts.deterministic  = !fast;
    if (dec_off)   opts.decimate.enabled = false;
    if (dec_force) opts.decimate.faces_threshold = 0;

    metriko::MatXd V;
    metriko::MatXi F;
    if (!igl::readOBJ(input, V, F) || V.rows() == 0 || F.rows() == 0) {
        std::cerr << "error: cannot read a triangle mesh from '" << input << "'\n";
        return EXIT_USAGE_IO;
    }

    if (decimate_only) { // no remeshing: useful to inspect or reuse the mesh the pipeline would work on
        try {
            opts.decimate.enabled = true;
            auto [V2, F2] = metriko::decimate_mesh(V, F, opts.decimate);
            std::FILE* f = std::fopen(output.c_str(), "wb");
            if (!f) { std::cerr << "error: cannot write '" << output << "'\n"; return EXIT_USAGE_IO; }
            for (int i = 0; i < V2.rows(); ++i) std::fprintf(f, "v %.17g %.17g %.17g\n", V2(i, 0), V2(i, 1), V2(i, 2));
            for (int i = 0; i < F2.rows(); ++i) std::fprintf(f, "f %d %d %d\n", F2(i, 0) + 1, F2(i, 1) + 1, F2(i, 2) + 1);
            if (std::fclose(f) != 0) { std::cerr << "error: cannot write '" << output << "'\n"; return EXIT_USAGE_IO; }
            std::cout << "decimated: vertices: " << V2.rows() << "  faces: " << F2.rows() << '\n';
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "error: decimation failed: " << e.what() << '\n';
            return EXIT_ALGORITHM;
        }
    }

    const auto t0 = std::chrono::steady_clock::now();
    metriko::RemeshResult r;
    try {
        r = metriko::remesh(V, F, opts);
    } catch (const std::invalid_argument& e) {
        std::cerr << "error: invalid input: " << e.what() << '\n';
        return EXIT_USAGE_IO;
    } catch (const std::exception& e) {
        std::cerr << "error: remeshing failed: " << e.what() << '\n';
        return EXIT_ALGORITHM;
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (!write_quads(output, r)) {
        std::cerr << "error: cannot write '" << output << "'\n";
        return EXIT_USAGE_IO;
    }
    if (!patches.empty() && !write_patches(patches, r)) {
        std::cerr << "error: cannot write '" << patches << "'\n";
        return EXIT_USAGE_IO;
    }
    if (!decimated_out.empty()) {
        std::FILE* f = std::fopen(decimated_out.c_str(), "wb");
        if (!f) { std::cerr << "error: cannot write '" << decimated_out << "'\n"; return EXIT_USAGE_IO; }
        for (int i = 0; i < r.used_pos.rows(); ++i)
            std::fprintf(f, "v %.17g %.17g %.17g\n", r.used_pos(i, 0), r.used_pos(i, 1), r.used_pos(i, 2));
        for (int i = 0; i < r.used_idx.rows(); ++i)
            std::fprintf(f, "f %d %d %d\n", r.used_idx(i, 0) + 1, r.used_idx(i, 1) + 1, r.used_idx(i, 2) + 1);
        if (std::fclose(f) != 0) { std::cerr << "error: cannot write '" << decimated_out << "'\n"; return EXIT_USAGE_IO; }
    }
    std::cout << "quads: " << r.idx.rows() << "  vertices: " << r.pos.rows()
              << "  scale: " << r.used_scale << "  elapsed: " << secs << " s\n";
    if (r.holes > 0)
        std::cerr << "warning: " << r.holes << " quad(s) could not be extracted (holes): the parametrization has\n"
                  << "         folded or collapsed uv triangles there; the rest of the mesh is valid\n";
    return 0;
}
