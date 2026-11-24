#ifndef METRIKO_LIB_H
#define METRIKO_LIB_H
#include <expected>
#include "metriko/api.h"
#include "hmesh/hmesh.h"
#include "nvec/face_rosy_field.h"
#include "igm/parameterization.h"
#include "quantization/quantization.h"
#include "tmesh/emesh.h"
#include "tmesh/emesh_validate.h"
#include "metriko/tutte/tutte_cutting.h"
#include "metriko/tutte/tutte_params.h"
#include "metriko/tutte/tutte_slim.h"
#include "metriko/qex/sanitization.h"
#include "metriko/qex/gen_q_vert.h"
#include "metriko/qex/gen_q_port.h"
#include "metriko/qex/gen_q_edge.h"
#include "metriko/qex/gen_q_face.h"
#include "patching.h"

namespace metriko {
struct QuadrangulationParts {
    std::unique_ptr<Hmesh> hmesh;
    std::unique_ptr<Mgrph> mgrph;
    std::unique_ptr<Emesh> emesh;
    std::unique_ptr<Hmesh> hmesh_cut; // the mesh with the t-mesh embedded, cut along the seam
    MatXd                  uv;        // per vertex of hmesh_cut: the slim result
    vec<qex::Qport>        q_ports;
    vec<qex::Qedge>        q_edges;
    vec<qex::Qface>        q_faces;
    MatXd                  q_pos;
    MatXi                  q_idx;
    MatXi                  q_val;
    MatXi                  q_fill; // triangles that close the boundary loops left by the dropped q-face rings
    int                    holes = 0;
};

inline QuadrangulationParts compute_quadrangulation_impl(
    const MatXd& V,
    const MatXi& F,
    const double scale,
    const bool calign = false,
    const bool refine = true,
    const int slim_iter = RemeshOptions{}.slim_max_iter,
    const double slim_tol = RemeshOptions{}.slim_rel_tol
) {
    detail::StageClock clock;
    clock.mark(RemeshStage::field);
    auto unit = grid_unit(V, scale);
    auto hm   = std::make_unique<Hmesh>(MatXd(V / unit), F);
    auto rawf = FaceRosyField(*hm, 4, calign? FieldType::CurvatureAligned : FieldType::Smoothest);
    auto seam = compute_seam(rawf);
    auto cutm = compute_cut_mesh(*hm, seam);
    auto cmbf = compute_combbed_field(rawf, seam);
    auto extf = compute_extrinsic_field(*cmbf, 4);

    clock.mark(RemeshStage::parameterization);
    auto rp = RosyParameterization(*hm, *cutm, extf, cmbf->singular, cmbf->matching, seam, 4, scale);
    rp.localInjectivity = true;
    rp.setup();
    rp.integ();

    hm->cfn = VecXc(hm->nF * 3);
    for (const Face f: hm->faces) {
        hm->cfn(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        hm->cfn(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        hm->cfn(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }

    clock.mark(RemeshStage::tmesh);
    auto mg = std::make_unique<Mgrph>(*hm, cmbf->matching, cmbf->singular);
    auto tm = std::make_unique<Emesh>(*mg);

    clock.mark(RemeshStage::quantization);
    auto X  = compute_quantization(*tm, *mg);
    validate_quantization(*tm, X);
    tm->set_x(X);

    clock.mark(RemeshStage::collapse);
    for (int i = 0; i < 100; ++i) {
        for (Ehalf th0 : tm->thalfs) {
            if (th0.id == -1 || th0.x != 0) continue;
            const auto& th1 = th0.twin();
            if (th0.tquad().thids(th0.side()).size() == 1) continue;
            if (th1.tquad().thids(th1.side()).size() == 1) continue;
            tm->collapse_thalf(th0.id);
        }

        for (const auto& [tqid, _] : tm->live_tquads())
            if (Tqaux a; tm->collapse_tquad_prepare(tqid, a)) tm->collapse_tquad_execute(a);

        if (rg::none_of(tm->thalfs, [](const Ehalf& th) { return th.id != -1 && th.x == 0; })) break;
    }
    validate_collapse_done(*tm);

    tm->collapse_tedge_snap(false);
    tm->collapse_tedge_snap(true);
    for (const auto& [teid, _] : tm->live_tedges()) { tm->collapse_tedge_snap_dedup(teid); }

    repair_crossing_tedges(*tm);

    clock.mark(RemeshStage::embedding);
    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(*hm, *tm, seam, cmbf->matching, cmbf->singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);

    clock.mark(RemeshStage::tutte);
    MatXd uv_ebd = compute_tutte_parameterization(*hm_emb, *tm, seam1, hdata);
    MatXd uv_cut(hm_cut->nV, 2);
    for (auto v: hm_cut->verts) uv_cut.row(v.id) = uv_ebd.row(v.half().next().crnr().id);
    clock.mark(RemeshStage::slim);
    MatXd uv = seamless_slim(*hm_emb, *hm_cut, seam1, singular1, uv_cut, slim_iter, slim_tol);

    // ------ qex on the slim result ------
    hm_emb->cfn = VecXc(hm_emb->nF * 3);
    for (int i = 0; i < hm_cut->nF; ++i) {
    for (int j = 0; j < 3; ++j) {
        int k = hm_cut->idx(i, j);
        hm_emb->cfn(i * 3 + j) = complex(uv(k, 0), uv(k, 1));
    }}

    clock.mark(RemeshStage::qex);
    qex::sanitization(*hm_emb, matching1, singular1, 4);

    vec<qex::Qport> q_ports;
    vec<qex::Qvert> vqvs, eqvs, fqvs;
    qex::generate_q_vert(*hm_emb, vqvs, eqvs, fqvs);
    qex::generate_vqvert_qport(*hm_emb, vqvs, q_ports);
    qex::generate_eqvert_qport(*hm_emb, eqvs, q_ports);
    qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);
    auto q_edges = qex::generate_q_edge(*hm_emb, matching1, q_ports);
    int n_open = 0, n_incomplete = 0;
    auto q_faces = qex::generate_q_faces(q_ports, q_edges, &n_open, &n_incomplete);
    const int n_dropped = n_open + n_incomplete;
    const double hole_ratio = n_dropped == 0 ? 0.
                             : double(n_dropped) / double(n_dropped + static_cast<int>(q_faces.size()));
    METRIKO_CHECK(!q_faces.empty() && hole_ratio <= 0.02,
        "{} of {} q-face rings do not close ({} dangling q-ports, {:.2f}% of the quads lost): the uv "
        "parametrization here has folded or collapsed triangles the extractor cannot trace through",
        n_dropped, n_dropped + static_cast<int>(q_faces.size()), q_ports.size() - static_cast<int>(q_edges.size()) * 2,
        100.0 * hole_ratio);
    MatXi q_fill;
    auto [q_pos, q_idx] = extract_quad_mesh(*hm, q_faces, refine, &q_fill);
    // the patch labels are replayed through the quad mesh's edge graph, which needs a closed mesh: with holes
    // there is no well-defined patch there, so label every quad with -1
    auto q_val = n_dropped == 0 ? label_quad_patches(*tm, cmbf->singular, q_faces, q_idx)
                                : MatXi::Constant(q_idx.rows(), 1, -1);
    q_pos *= unit;

    return {
        std::move(hm),
        std::move(mg),
        std::move(tm),
        std::move(hm_cut),
        std::move(uv),
        std::move(q_ports),
        std::move(q_edges),
        std::move(q_faces),
        std::move(q_pos),
        std::move(q_idx),
        std::move(q_val),
        std::move(q_fill),
        n_dropped
    };
}

inline std::expected<RemeshResult, std::string> compute_quadrangulation(
    const MatXd& V,
    const MatXi& F,
    const double s
) {
    try {
        auto q = compute_quadrangulation_impl(V, F, s);
        return RemeshResult{
            std::move(q.q_pos),
            std::move(q.q_idx),
            std::move(q.q_val),
            {},
            {},
            s,
            q.holes,
            std::move(q.q_fill)
        };
    }
    catch (const std::exception& ex) {
        return std::unexpected(std::string(stage_name(detail::failed_stage())) + ": " + ex.what());
    }
    catch (...) { return std::unexpected("unknown error"); }
}
}
#endif
