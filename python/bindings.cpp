#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>
#include <cstdint>
#include <stdexcept>
#include <string>
#include "metriko/api.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {
using MatX3dRow = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using MatX4iRow = Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor>;
using MatX3iRow = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;

using DoubleMat = nb::ndarray<const double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using IntMat    = nb::ndarray<const std::int32_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

metriko::RemeshField parse_field(const std::string& s) {
    if (s == "smoothest")          return metriko::RemeshField::Smoothest;
    if (s == "curvature_aligned")  return metriko::RemeshField::CurvatureAligned;
    throw std::invalid_argument("field must be 'smoothest' or 'curvature_aligned' (got '" + s + "')");
}

template <class T>
void set_if(T& dst, const std::optional<T>& v) { if (v) dst = *v; }

nb::tuple remesh(DoubleMat V, IntMat F, double scale, const std::string& field, bool refine,
                 int slim_max_iter, double slim_rel_tol, int threads, bool quiet,
                 const std::optional<bool>& decimate,
                 std::optional<double> decimate_verts_pct,
                 std::optional<int>    decimate_max_verts,
                 std::optional<int>    decimate_faces_threshold,
                 std::optional<double> decimate_normal_dev,
                 std::optional<double> decimate_aspect_ratio,
                 std::optional<double> decimate_min_edge_len,
                 std::optional<int>    decimate_max_valence,
                 std::optional<double> decimate_hausdorff,
                 std::optional<double> decimate_seam_thresh,
                 std::optional<double> decimate_seam_dev,
                 bool return_decimated,
                 bool retry, int max_attempts, double coarsen_factor, bool timings, double time_budget, double stage_timeout, bool return_stats) {
    metriko::RemeshOptions opts;
    opts.scale         = scale;
    opts.field         = parse_field(field);
    opts.refine        = refine;
    opts.slim_max_iter = slim_max_iter;
    opts.slim_rel_tol  = slim_rel_tol;
    opts.threads       = threads;
    opts.quiet         = quiet;
    opts.retry_coarser = retry;
    opts.max_attempts  = max_attempts;
    opts.coarsen_factor = coarsen_factor;
    opts.report_timings = timings;
    opts.time_budget_seconds = time_budget;
    opts.stage_timeout_seconds = stage_timeout;

    auto& d = opts.decimate;
    if (decimate && !*decimate) d.enabled = false;
    if (decimate && *decimate)  d.faces_threshold = 0;
    set_if(d.vertices_percentage, decimate_verts_pct);
    set_if(d.max_verts,           decimate_max_verts);
    set_if(d.faces_threshold,     decimate_faces_threshold);
    set_if(d.normal_deviation,    decimate_normal_dev);
    set_if(d.aspect_ratio,        decimate_aspect_ratio);
    set_if(d.min_edge_length,     decimate_min_edge_len);
    set_if(d.max_valence,         decimate_max_valence);
    set_if(d.hausdorff_error,     decimate_hausdorff);
    set_if(d.seam_threshold,      decimate_seam_thresh);
    set_if(d.seam_deviation,      decimate_seam_dev);

    const metriko::MatXd Vm = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(V.data(), V.shape(0), V.shape(1));
    const metriko::MatXi Fm = Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(F.data(), F.shape(0), F.shape(1));

    metriko::RemeshResult r;
    {
        nb::gil_scoped_release release;
        r = metriko::remesh(Vm, Fm, opts);
    }
    if (r.idx.cols() != 4) throw std::runtime_error("internal error: the output is not a quad mesh");
    MatX3dRow pos = r.pos;
    MatX4iRow idx = r.idx;
    Eigen::VectorXi patch = r.val.col(0);
    nb::dict stats;
    if (return_stats) {
        stats["used_scale"] = r.used_scale;
        stats["holes"]      = r.holes;
        stats["fill_idx"]   = r.fill_idx; // triangles that close the loops left by the dropped q-face rings
    }
    if (!return_decimated)
        return return_stats ? nb::make_tuple(std::move(pos), std::move(idx), std::move(patch), std::move(stats))
                            : nb::make_tuple(std::move(pos), std::move(idx), std::move(patch));
    MatX3dRow used_pos = r.used_pos;
    MatX3iRow used_idx = r.used_idx;
    if (return_stats)
        return nb::make_tuple(std::move(pos), std::move(idx), std::move(patch), std::move(used_pos), std::move(used_idx), std::move(stats));
    return nb::make_tuple(std::move(pos), std::move(idx), std::move(patch), std::move(used_pos), std::move(used_idx));
}
}

NB_MODULE(_core, m) {
    m.doc() = "metriko quad remeshing core";
    nb::exception<metriko::RemeshFailure>(m, "RemeshError", PyExc_RuntimeError);
    m.def("remesh", &remesh,
          "V"_a, "F"_a, "scale"_a = 0.02, "field"_a = "smoothest", "refine"_a = true,
          "slim_max_iter"_a = 30, "slim_rel_tol"_a = 1e-4, "threads"_a = 0, "quiet"_a = false,
          "decimate"_a = nb::none(), "decimate_verts_pct"_a = nb::none(),
          "decimate_max_verts"_a = nb::none(), "decimate_faces_threshold"_a = nb::none(),
          "decimate_normal_dev"_a = nb::none(), "decimate_aspect_ratio"_a = nb::none(),
          "decimate_min_edge_len"_a = nb::none(), "decimate_max_valence"_a = nb::none(),
          "decimate_hausdorff"_a = nb::none(), "decimate_seam_thresh"_a = nb::none(),
          "decimate_seam_dev"_a = nb::none(), "return_decimated"_a = false,
          "retry"_a = true, "max_attempts"_a = 3, "coarsen_factor"_a = 1.6, "timings"_a = false,
          "time_budget"_a = 300.0, "stage_timeout"_a = 150.0,
          "return_stats"_a = false);
}
