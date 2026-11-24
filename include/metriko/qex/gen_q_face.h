//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QEX_GEN_Q_FACE_H
#define METRIKO_QEX_GEN_Q_FACE_H
#include "common.h"

namespace metriko::qex {
    // walks the q-edge pairing around every q-vertex. a ring that cannot be walked entirely (a dangling q-port
    // in the way) or whose four edges do not close on the port they started from is *dropped*, not fatal: the
    // reference implementation (libQEx, Ebke et al. 2013, §4.6) treats exactly this case as a warning and leaves
    // a hole, because it is a local defect of the parametrization, not a broken extractor. the caller counts the
    // drops and decides whether they are few enough to return a mesh (lib.h)
    inline vec<Qface> generate_q_faces(
        const vec<Qport> &qps,
        const vec<Qedge> &qes,
        int* n_open       = nullptr, // rings whose four q-edges do not close
        int* n_incomplete = nullptr  // rings the walk could not finish (a dangling q-port)
    ) {
        vec<Qface> qfs;
        vec<Qhalf> qhs;
        for (int i = 0; i < qes.size(); i++) {
            auto &qe = qes[i];
            qhs.emplace_back(qe, i * 2 + 0, true);
            qhs.emplace_back(qe, i * 2 + 1, false);
        }

        vec po2qh(qps.size(), -1);
        vec visit(qhs.size(), false);

        for (const auto &qh: qhs) po2qh[qh.port1().idx] = qh.idx;

        for (size_t i = 0; i < qhs.size(); ++i) {
            if (visit[i]) continue;
            vec qfhs{qhs[i]};
            visit[i] = true;
            auto *qp = &qhs[i].port2();
            for (int i = 0; i < 3; i++) {
                int h = po2qh[qps[qp->prev_id].idx]; if (h < 0) break;
                qfhs.emplace_back(qhs[h]);
                visit[h] = true;
                qp = &qhs[h].port2();
            }
            if (qfhs.size() != 4) { if (n_incomplete) ++*n_incomplete; continue; }
            if (qps[qp->prev_id].idx != qhs[i].port1().idx) { if (n_open) ++*n_open; continue; }
            qfs.emplace_back(qfhs);
        }
        return qfs;
    }
}
#endif
