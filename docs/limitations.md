# Limitations — explained from the ground up

This page is written for someone who does not work with mesh processing every day. It first explains the terms
metriko uses (a small glossary), then states what metriko can and cannot do, in the same terms.

If you already know what a manifold mesh, a T-mesh and a seamless parametrization are, jump to
[The limitations](#the-limitations).

---

## Glossary

### Mesh, triangle, vertex, edge, face

A **mesh** is a surface represented by small polygons. In metriko it is always a **triangle mesh**: the surface is
covered by flat triangles that share their corners and edges.

* A **vertex** is a corner point.
* An **edge** is a segment joining two vertices.
* A **face** (here: a triangle) is bounded by three edges.

"Triangle mesh" is just the raw input. metriko's output is a **quad mesh**: the faces are quadrilaterals instead.

### Boundary, watertight / closed

Look at every edge and count how many triangles use it.

* Any edge used by **exactly 2** triangles is an interior edge: the surface continues on both sides.
* Any edge used by **exactly 1** triangle is a **boundary edge**: the surface stops there, like the rim of a sheet
  of paper.

A mesh with no boundary edges is called **closed** or **watertight** (think of a balloon: no opening anywhere).
A mesh with boundary edges is **open** (a sheet, a scan with a hole).

metriko **only accepts closed meshes**, because its whole pipeline (the parametrization must cover the entire
surface, the grid must wrap around) assumes there is no rim. If your mesh has a hole, cap it first (most DCC tools
and MeshLab have a "fill holes" operation).

### Manifold

A surface is **manifold** when, zoomed in enough, every point looks like a small flat disk. Two common violations:

* an edge used by **3 or more** triangles (three sheets meeting along the same seam);
* a **pinch vertex** ("bowtie"): two otherwise separate pieces of surface touching at a single vertex. Near that
  vertex the surface looks like two cones glued at their tips, not like a disk.

metriko rejects non-manifold meshes because the integer grid and the quad layout are only well defined on a
manifold surface.

### Connected components — "a single-component mesh"

Two triangles **belong to the same component** if you can walk from one to the other through edges (face →
neighbour face). A model of two separate objects (say, a cup and a saucer exported together) has **two
components**.

metriko processes **exactly one component**. If your file has several objects, split it first (or remesh them one
by one). This is what "single component" means in the error messages; it is not about the object having "parts" in
the design sense — a cup with a handle is one component.

### Genus (number of handles)

For a closed surface, the **genus** is the number of "handles" you could pass a finger through: a sphere has
genus 0, a donut (torus) has genus 1, a two-handled mug has genus 2.

It matters because a closed surface of genus *g* satisfies the **Euler formula**

```
vertices − edges + faces = 2 − 2·g
```

which is how we verify that an output mesh has the same topology as the input (for a genus-0 input: 2). metriko
handles genus > 0, but the more handles, the more singularities the field needs.

### Quad mesh, quad-dominant

A **quad mesh** is a mesh whose faces have four sides. Quad meshes are preferred for subdivision surfaces,
texturing and animation because the edge flow follows the shape.

**Quad-dominant** means "mostly quads, with a few triangles". metriko produces pure quads when it can, and a
small number of fill triangles when the extraction has to drop a few cells (see
[Quad-dominance](#quad-dominance-not-always-pure-quads)).

### Remeshing, scale

**Remeshing** means rebuilding the connectivity of a mesh while keeping its shape: same surface, different
triangles/quads.

**`scale`** is the target size of the quads, expressed as a fraction of the model's **bounding-box diagonal** (the
diagonal of the smallest box containing the model). `scale = 0.02` means "quads about 1/50 of the model's
diagonal", i.e. roughly 50 quads across the model. It is relative, so it does not depend on whether your model is
in millimetres or metres.

### Cross field (direction field)

A **cross field** assigns, to every point of the surface, a set of four perpendicular directions ("a cross"), like
the grain of wood but four-way. The quads are aligned with this field, which is why the quad edges follow the
shape's natural directions instead of an arbitrary grid.

### Singularity

At most points the four cross directions can be followed smoothly from a triangle to its neighbours. At a few
points they cannot: going around that point, the cross comes back rotated by 90°, 180° or 270°. Such a point is a
**singularity** (a "cone" point). Singularities are unavoidable on a closed surface (a combed hairy ball always
has a cowlick), and they are where the quad grid has an irregular vertex (valence ≠ 4).

metriko's `holes`/patch labels and the number of singularities are related: more singularities mean more
constraints on the grid.

### Seamless parametrization, cut, uv

To flatten a closed surface you must **cut** it open somewhere (like peeling an orange). The cut curve is the
**seam**: along it, the flattened layout can be shifted/rotated without changing the 3D result.

A **parametrization** assigns a 2D coordinate (u, v) to every point. **Seamless** means the 2D layout is consistent
across the cut up to an integer grid translation and a 90° rotation — this is what allows the integer grid to be
traced across seams without the quads becoming skewed.

### Integer grid, quantization

The quads come from overlaying an **integer grid** (like squared paper) on the flattened surface: every cell of
the grid becomes a quad.

The parametrization is continuous, so the grid positions are *real* numbers. **Quantization** is the step that
moves them to *integers* so the grid cells become actual quads instead of "almost-quads". This is a combinatorial
optimization and it is the stage that most often dominates the time on hard meshes.

### T-mesh (the concept you asked about)

The **T-mesh** is the intermediate structure between the continuous parametrization and the final quad mesh. Think
of it as *a quad mesh with T-junctions*:

* its edges are **grid arcs**: pieces of the integer grid lines, running between singularities, seams and
  crossings;
* where two arcs cross, or where an arc ends on another arc, there is a **T-junction** (the "T" of T-mesh: an arc
  can terminate on the side of another arc instead of continuing to the next crossing);
* the regions bounded by arcs are **t-quads**: they are *supposed* to become quads, but in the T-mesh they are not
  yet square — their side lengths are still real numbers.

So the T-mesh is the **combinatorial skeleton** of the result: it already knows how many quads there will be and
how they connect, but not yet their exact sizes. Quantization then assigns an integer length to every arc, and a
t-quad is only valid if its four sides add up consistently (it must become a unit square). After quantization, a
**collapse** stage removes the t-quads that turned out to have zero width, and the remaining grid is the quad mesh
you get.

Why you should care: **almost every "hard mesh" failure happens in the T-mesh**, because when the integer grid
becomes very dense or the field is very irregular, the T-mesh can become *degenerate*:

* a t-quad whose two opposite arcs end up at the same place (**zero-width t-quad**);
* an arc that crosses another one in a way that makes a t-quad **self-intersect**;
* a segment whose two ends are no longer connected through a face of the mesh (the error
  `boundary segment ... lies in no face`).

Those are the cases metriko reports as `collapse` / `embedding` errors.

### Quantization, pricing (the cost you see in the logs)

Quantization works by repeatedly improving the integer assignment. The **pricing** loop looks for a *strip loop*
(a closed chain of arcs along the grid) whose cost would decrease if we add +1 or −1 along it, applies that
change, and repeats. Each step is cheap; the number of steps is **data dependent**, which is why a small mesh can
sometimes spend far longer in `quantization` than a larger one.

---

## The limitations

### Input domain

metriko is built for **closed, manifold, single-component triangle meshes**. Each of those words means what the
glossary says:

| requirement | rejected when | how to fix your mesh |
|---|---|---|
| closed / watertight | any edge used by exactly 1 triangle | fill the holes |
| manifold | any edge used by 3+ triangles, or a pinch vertex | delete/repair the offending faces |
| single component | two face-disconnected pieces | split the file, process each part |
| triangles | quads/ngons in the C++ API | triangulate (the Python loader does it for you) |

Rejection happens at the `input` stage, with `std::invalid_argument` (C++) / `ValueError` (Python), before any
expensive work. metriko never silently repairs your mesh: a wrong input would produce a wrong output.

Not required: UVs, normals, a specific scale, a specific orientation, or a specific face order.

### Meshes that fail structurally (the "degenerate T-mesh" class)

Some meshes produce a degenerate T-mesh (see the glossary) at every resolution tried. On our reference set this is
one mesh out of nine: `dragonita` (498k faces). It failed at three *different* stages across three retry attempts:

| attempt | configuration | outcome |
|---|---|---|
| 1 | scale 0.02, 90k vertices | `embedding: no half between vertices` |
| 2 | scale 0.032, 54k vertices | `collapse: left 2 zero-width tquads` |
| 3 | scale 0.0512, 39k vertices | `collapse: tquad 869 is self-intersected by one of its thalfs` |

This is not a single bug to patch: it is a robustness limit of the T-mesh construction. Such meshes are reported
with a stage-tagged error instead of producing a broken quad mesh. If you meet one, the practical options are to
decimate more aggressively (`--decimate-max-verts`), to repair the mesh (welding near-duplicate vertices and
removing sliver triangles often removes the degeneracy), or to accept a different tool for that model.

### Quad-dominance, not always pure quads

The extractor walks the integer grid to build quads. Occasionally a ring does not close — the flattened (u, v)
layout is locally folded, so the walk cannot trace a closed cell. Those rings are dropped (the behaviour of the
reference QEx implementation), counted in `holes`, and the holes they leave are **closed with fill triangles**.

Consequences:

* `holes == 0` → the output is **pure quads**;
* `holes > 0` → the output is **quad-dominant**, with a number of triangles proportional to the dropped rings
  (on our worst case: 22 triangles for 9,206 quads, 99.8 % quads);
* if more than 2 % of the rings fail, the run **fails** instead: that means the parametrization is broken, not a
  few local cells.

### Reproducibility

With the defaults (`deterministic=True`), metriko removes the scheduling-dependent floating-point order that it
controls (the SLIM assembly, the CHOLMOD factorization threading, Eigen's internal parallelism). Two runs of the
same mesh, on the same machine, with the same thread count, give the same quad count, vertex count and hole
count. That is verified, not a mathematical guarantee:

* a different thread count, a different CPU (different SIMD paths in the linear algebra library) or `--fast` can
  change the last bits, and on a *marginal* mesh (one that sits close to a degenerate configuration) that can
  change the outcome;
* wall-clock time is not reproducible at all: it depends on the machine, its load, its temperature, and on the
  time budget (which is wall-clock based).

`--threads 1` is the strongest reproducibility setting. If you need to reproduce a specific result, record the
version, the options, the thread count and the machine.

### Time

The defaults are **300 s total** and **150 s per stage**. Those are *caps*, not a service-level guarantee:

* on an idle machine, every mesh in our 8-mesh reference set finishes under 300 s (19 s to 261 s);
* on a loaded machine the same runs take 2-2.5× longer, and the meshes in the 150-280 s range cross the cap.

The budget is checked at stage boundaries, so a stage that starts before the deadline is allowed to finish (with a
10 % grace): a run can exceed the budget by up to about one stage. Use `--timings` to see where the time goes, and
raise `--time-budget`/`--stage-timeout` to match your hardware.

### Memory

The parameterization dominates: it assembles and factorizes a large sparse system whose size grows with the vertex
count (about 18 unknowns per vertex in the current formulation). The peak grows with it and with the fill of the
factorization; we have not published a memory benchmark. `decimate.max_verts` is the knob that lowers it. There is
no streaming or out-of-core mode.

### Platform

The distributed binaries are:

* Linux x86_64 (`manylinux_2_28`), Windows x86_64, macOS arm64;
* CPython ≥ 3.12, `abi3` (one wheel for 3.12, 3.13, …);
* the x86_64 builds assume **AVX2/FMA** (the bundled OpenBLAS targets `HASWELL`). On an older CPU, import fails
  with an illegal-instruction error; build from source with a different OpenBLAS target in that case.

### Not implemented

* Other quadrangulation algorithms (only the QGP pipeline with QEx extraction is exposed).
* Anisotropic or feature-preserving quad sizing: a single scalar `scale` for the whole mesh.
* Texture atlases / UV output: the patch labels are exposed, the uv coordinates are not.
* Partial or incremental remeshing of a region: the whole mesh is always processed.
* Open surfaces in any form.
* Multi-component inputs in a single call.

---

## What to do when something fails

1. **Read the stage.** The error message starts with the stage name; [errors.md](errors.md) has a table of what
   each stage means and the first thing to try.
2. **Look at `--timings`.** If the failing stage was over its timeout, the mesh is simply too big for the budget:
   lower `--decimate-max-verts` or raise the budget.
3. **Check `used_scale`.** If the retry ladder produced the result, the density is coarser than you asked for.
   Raise `--max-attempts` or `--coarsen-factor` if you need to know.
4. **For `collapse`/`embedding` errors**, try a coarser `--scale` or a different decimation. This class is
   mesh-specific and resolution-dependent; a 20 % coarser grid often escapes it.
5. **If the mesh is precious and keeps failing**, repair it outside metriko (weld duplicates, remove degenerate
   triangles) and retry — those repairs address the root cause of most degenerate T-meshes.
