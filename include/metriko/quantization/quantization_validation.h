//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef METRIKO_QUANTIZATION_VALIDATION_H
#define METRIKO_QUANTIZATION_VALIDATION_H
#include <queue>
#include "metriko/tmesh/emesh.h"

namespace metriko {
// the zero-width t-edges must not form a path between two fixed junction nodes (a loop would collapse onto
// nothing). only the nodes a zero-width tedge touches are considered: the check runs once per step of the
// pricing loop, so it has to be O(#zero tedges), not O(#nodes)
inline bool compute_validation(
    const Mgrph &mg,
    const Emesh &tm,
    const VecXd &X
) {
    if ((X.array() < 0).any()) return false;

    umap<int, int> idx;   // tnode -> compact index
    vec<int>       nids;  // compact index -> tnode
    vec<vec<int>>  adjs;
    auto intern = [&](int nid) {
        auto it = idx.find(nid);
        if (it != idx.end()) return it->second;
        const int i = static_cast<int>(nids.size());
        idx.emplace(nid, i);
        nids.push_back(nid);
        adjs.emplace_back();
        return i;
    };
    for (const auto& [teid, ns]: tm.live_tedges()) {
        if (X[teid] != 0) continue;
        const int a = intern(ns.front());
        const int b = intern(ns.back());
        adjs[a].push_back(b);
        adjs[b].push_back(a);
    }

    vec<char> visited(adjs.size(), 0);
    for (int s = 0; s < static_cast<int>(adjs.size()); ++s) {
        if (visited[s]) continue;
        int n_fixed = 0;
        std::queue<int> q;
        q.push(s);
        visited[s] = 1;
        while (!q.empty()) {
            const int c = q.front();
            q.pop();
            if (mg.mnodes[nids[c]].jt == JunctionType::F) ++n_fixed;
            for (int n: adjs[c]) if (!visited[n]) { visited[n] = 1; q.push(n); }
        }
        if (n_fixed > 1) return false;
    }
    return true;
}
}
#endif
