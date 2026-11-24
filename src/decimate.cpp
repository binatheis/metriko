#include "metriko/api.h"

#include <pmp/surface_mesh.h>
#include <pmp/algorithms/decimation.h>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace metriko {
std::pair<MatXd, MatXi> decimate_mesh(const MatXd& V, const MatXi& F, const DecimateOptions& o) {
    pmp::SurfaceMesh mesh;
    std::vector<pmp::Vertex> verts(static_cast<size_t>(V.rows()));
    for (int i = 0; i < V.rows(); ++i)
        verts[i] = mesh.add_vertex(pmp::Point(V(i, 0), V(i, 1), V(i, 2)));
    for (int i = 0; i < F.rows(); ++i)
        mesh.add_triangle(verts[F(i, 0)], verts[F(i, 1)], verts[F(i, 2)]);
    mesh.garbage_collection();

    const auto nv = static_cast<unsigned>(V.rows());
    unsigned target;
    if (o.vertices_percentage > 0) {
        target = std::clamp<unsigned>(
            static_cast<unsigned>(std::lround(nv * o.vertices_percentage / 100.0)),
            1, nv);
    } else {
        target = o.max_verts > 0 ? std::min<unsigned>(static_cast<unsigned>(o.max_verts), nv) : nv;
    }
    if (target >= nv) return {V, F};

    try {
        pmp::decimate(mesh, target,
                      static_cast<pmp::Scalar>(o.aspect_ratio),
                      static_cast<pmp::Scalar>(o.min_edge_length),
                      static_cast<unsigned>(std::max(0, o.max_valence)),
                      static_cast<pmp::Scalar>(o.normal_deviation),
                      static_cast<pmp::Scalar>(o.hausdorff_error),
                      // <0: "disabled"; pmp only consults it when the mesh carries h:tex, which we never load
                      static_cast<pmp::Scalar>(o.seam_threshold >= 0 ? o.seam_threshold : 1e-2),
                      static_cast<pmp::Scalar>(o.seam_deviation));
    }
    catch (const std::exception& ex) {
        throw std::runtime_error(std::string("decimation failed: ") + ex.what());
    }

    MatXd V2(static_cast<long>(mesh.n_vertices()), 3);
    MatXi F2(static_cast<long>(mesh.n_faces()), 3);
    for (auto v : mesh.vertices()) {
        const pmp::Point& p = mesh.position(v);
        V2.row(static_cast<long>(v.idx())) = Eigen::RowVector3d(
            static_cast<double>(p[0]), static_cast<double>(p[1]), static_cast<double>(p[2]));
    }
    for (auto f : mesh.faces()) {
        int j = 0;
        for (auto v : mesh.vertices(f)) F2(static_cast<long>(f.idx()), j++) = static_cast<int>(v.idx());
    }
    return {std::move(V2), std::move(F2)};
}
}
