# The output

## What comes back

### C++ — `metriko::RemeshResult`

| field | type | meaning |
|---|---|---|
| `pos` | `MatXd (k, 3)` | quad vertices |
| `idx` | `MatXi (q, 4)` | four vertex indices per quad, counter-clockwise seen from outside |
| `val` | `MatXi (q, 1)` | patch label of every quad (`-1` when rings were dropped) |
| `fill_idx` | `MatXi (t, 3)` | triangles that close the loops left by the dropped rings |
| `holes` | `int` | number of q-face rings the extractor dropped (informative) |
| `used_scale` | `double` | scale the returned mesh was built with (differs on a retry) |
| `used_pos` / `used_idx` | | the triangle mesh that was actually remeshed (post-decimation) |

### Python — `metriko.remesh(...)`

Returns `(Q_pos, Q_idx, patch)`. With `return_decimated=True` it appends `(V_used, F_used)`; with
`return_stats=True` it appends a dict with `used_scale`, `holes` and `fill_idx`.

Note that the helper `metriko.save_obj` writes **quads only**. If `stats["holes"] > 0`, the mesh is closed by the
triangles in `stats["fill_idx"]`; write those too (or use the CLI, which writes both in one file):

```python
Q_pos, Q_idx, patch, stats = metriko.remesh(V, F, return_stats=True)
fill = stats["fill_idx"]
with open("quads.obj", "w") as fh:
    for x, y, z in Q_pos:
        fh.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
    for a, b, c, d in Q_idx:
        fh.write(f"f {a+1} {b+1} {c+1} {d+1}\n")
    for a, b, c in fill:
        fh.write(f"f {a+1} {b+1} {c+1}\n")
```

## Closedness and manifoldness

The output is verified to be **watertight and manifold** on every mesh in our test set: each undirected edge is
shared by exactly two faces, so there are no boundary edges and no non-manifold edges, and the Euler
characteristic matches the input's topology (for a genus-0 input, χ = 2).

That property comes from two mechanisms:

1. the extraction produces quads on a closed integer grid, so the quad mesh is closed — except for the rings that
   cannot be traced through (dropped, see below);
2. `fill_holes` closes the boundary loops those drops leave: a 4-edge loop becomes a quad, longer loops are
   filled with a triangle fan from their centroid. Pinched boundaries (a vertex shared by two holes) are split
   into simple cycles first, so the fill stays manifold.

## Quad-dominance

* `holes == 0`: the output is **pure quads**.
* `holes > 0`: a few cells are filled with triangles. On our worst mesh (drag_30k, 25 dropped rings) the output
  was 12,609 quads with no fill at all; on Azure full it was 9,206 quads + 22 triangles (99.8 % quads). In other
  words: quad-dominant, and the triangle count stays proportional to the number of dropped rings.

## Patch labels

`patch` is the chart (patch) of every quad, useful for texturing: each patch is a disk in the parametrization, and
the seams between patches are where the integer grid jumps.

When rings are dropped, the labels cannot be transferred reliably and the whole array is set to `-1` (the
alternative would be silently wrong labels).

## File format

The CLI writes a single OBJ with `v` lines for the vertices, then `f a b c d` for the quads and `f a b c` for the
fill triangles (if any). Both are 1-based. No normals, no uv, no materials.

The Python front-end returns arrays; use `metriko.save_obj` for the simple case or write your own OBJ/PLY as above.

## Checking a result yourself

A 20-line check with numpy:

```python
import numpy as np
from collections import Counter

def edge_report(V, F):
    ec = Counter()
    for f in F:
        n = len(f)
        for i in range(n):
            a, b = int(f[i]), int(f[(i + 1) % n])
            ec[(min(a, b), max(a, b))] += 1
    return {
        "faces": len(F),
        "boundary_edges": sum(1 for c in ec.values() if c == 1),
        "nonmanifold_edges": sum(1 for c in ec.values() if c > 2),
        "euler": len(np.unique(F)) - len(ec) + len(F),
    }
```

`boundary_edges` and `nonmanifold_edges` must both be 0.
