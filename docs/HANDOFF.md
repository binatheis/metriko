# HANDOFF — metriko fork: robustness pass 2

**Date:** 2026-11-21 (second session)
**Repo:** `D:\onito\OneDrive\Desarrollo\meshora\Inspo\metriko` (fork of https://github.com/komietty/metriko, upstream HEAD `7932da9`)
**Status of working tree:** still **all uncommitted**. `docs/` is gitignored.

---

## 1. What this session changed (all in `include/metriko/`)

| # | File | Change | Why |
|---|---|---|---|
| 1 | `tmesh/emesh.h`, `tmesh/emesh_collapse_tq_prepare.ipp` | `collapse_tquad_prepare` is now two-phase: it retries with `keep_zero_ends`, which adds the endpoints of zero-width thalfs to `aux.pts` even when they are plain junctions (`n_adj_tqs == 2`, not fixed) | **fixes the zero-width "lens" deadlock** (old handoff §3.2 item 3): the zero thalf's endpoint was never in `aux.pts`, so `prepare` always returned false for both tquads of the lens. The first phase keeps the previous behaviour bit for bit, so currently-working collapses are untouched |
| 2 | `qex/sanitization.h` | the integer snap is now **conditional**: a vertex whose snap would collapse (`area == 0`) or flip an incident face (vs. the slim orientation) keeps its slim value. `sanitize_uv` factored out | the unconditional `tol = 1e-4` snap created 86 (manotas) / 292 (drag) **degenerate uv faces**, which broke the q-edge pairing (582 unpaired ports on manotas) |
| 3 | `common/utilities.h` | `calc_coefficient` no longer throws `degenerate uv triangle`: for a param triangle collapsed to a line it returns the minimum-norm barycentric coordinates (= orthogonal projection onto the line), and `(0,0)` for a point | QEx explicitly admits triangles "degenerated to a line or a point" (paper §2.1). The throw made such a mesh unfinishable; the projection is exact for the only points that can be traced there (points on the line) |
| 4 | `qex/gen_q_face.h` | the ring check now reports both ports, their uv, and what went wrong | "fail with a clear stage-specific message" instead of a bare `q-face ring does not close` |
| 5 | `common/typedef.h` | opt-in `METRIKO_DIAG_LOC` (compile with `-DMETRIKO_DIAG_LOC`) prefixes every `METRIKO_FAIL`/`METRIKO_CHECK` message with `file:line` | locating which stage threw without a debugger. Off by default, zero cost otherwise |

Carried over from the previous session (still uncommitted, now regression-tested): `hmloc.h` `try_get_left_face`, `hmflat.h` `flat_around` closing face, `tmesh/emesh_collapse_path.h` reflex-vertex unfold, `tmesh/emesh_collapse_th.ipp` twin-collapse fallback, `tutte_slim.h` `n_small` rename, `lib.h` slim plumbing, `GC_HAVE_SUITESPARSE`→`METRIKO_HAVE_SUITESPARSE`, `parameterization_setup.ipp` `vec<int>`.

## 1.5 Session 3: production hardening (performance + failure contract)

| Area | Change |
|---|---|
| `quantization/quantization_validation.h` | `compute_validation` is now **O(#zero-width t-edges)** instead of O(#nodes): it only indexes the nodes a zero tedge touches (the old version allocated a `vector<vector<int>>` of every mnodes entry and ran once per pricing step). **The single biggest speedup**: drag_30k quantization 344 s → 126 s, whole attempt 540 s → 227 s; Azure full (955 kV → 150 k) 663 s → **275 s**; GEO-head quantization ~250 s → 1.7 s |
| `qex/gen_q_face.h` | a q-face ring that cannot be walked (dangling q-port) or does not close is **dropped and counted**, not fatal — libQEx (Ebke et al. 2013, §4.6) treats exactly this as a warning and leaves a hole. `generate_q_faces(qps, qes, &n_open, &n_incomplete)` |
| `lib.h` | counts the dropped rings, fails the `qex` stage only when **> 2 %** of them are lost (message names the counts and the cause), and returns `holes` in the result. Patch labels are `-1` when there are holes (the label replay needs a closed quad mesh) |
| `qex/gen_q_edge.h` | one summary line for unpaired q-ports instead of one line per port |
| `api.h` | `RemeshStage` enum + `RemeshFailure{stage}`; `RemeshOptions::{retry_coarser, max_attempts, coarsen_factor, time_budget_seconds, report_timings}`; `RemeshResult::{used_scale, holes}`; a stage clock that records the failing stage, prints per-stage timings, and enforces the time budget between stages |
| `src/metriko_api.cpp` | retry loop: a failed stage is retried with `scale *= coarsen_factor` (dense integer-grid maps are the main source of degenerate configurations), only when the remaining budget can fit another attempt |
| `src/cli/main.cpp` | `--timings --time-budget (default 300) --no-retry --coarsen-factor --max-attempts --decimate-only`; prints the scale used and a warning with the hole count |
| `python/bindings.cpp`, `python/metriko/__init__.py` | `metriko.RemeshError` (base `RuntimeError`) with `.stage` as a `metriko.Stage` enum; `retry/max_attempts/coarsen_factor/timings/time_budget/return_stats` kwargs |
| `tutte/tutte_slim.h` | prints its iteration count, unknowns and final energy under `--timings` |`n| `igm/iterative_rounding/*` | the constant blocks of the constraint jacobian (`X2F*wInteg`, the `wClose` identity) are built once instead of on every rounding step, and `--timings` prints the rounding step count + the time spent in those blocks |

`--decimate-only` exists because `--decimated-out` is only written on success: it writes the decimated triangle mesh without remeshing (the fast way to build test meshes).



**The pipeline is not bit-reproducible** (OpenMP reductions in `tutti_slim`/`construct_generating_vectors`, multithreaded CHOLMOD/BLAS): a marginal mesh can flip between runs. Repeat counts below are real runs of the same binary and input.

| Mesh | Result | Runs |
|---|---|---|
| manotas_30k | **SUCCESS** 3464 quads / 3466 verts, 53-139 s | 4/5 default; **3/3 with `--threads 1`** (was `degenerate uv triangle`) |
| azure_30k | **SUCCESS** 2240 quads, 36-68 s | 3/3 (+1 in the batch) |
| gat_30k | FAIL 69-169 s: `q-face ring does not close` at port 5819 (uv -1,-58) | 2/2, same port |
| drag_30k | FAIL 540 s: `q-face ring does not close` (was the ZW lens) | 2/2 |
| gatubelo (250 kV) | FAIL 1109 s: `the seam transitions are not integer (off by 0.472)` (`tutte/tutte_slim.h:65`) | 1 |
| GEO-head_stylized (264 kV) | **SUCCESS** 2432 quads, 295 s | 1 (matches upstream's 2432) |
| manotas (999 kV) | FAIL 270 s: `the seam transitions are not integer (off by 0.540)` | 1 |
| Azure Alien Head Scan (1) (955 kV) | **SUCCESS** 7237 quads, 663 s (955027 → 150000 verts decimated) | 1 |
| dragonita (249 kV) | **not finished after 30 min** (killed): stuck in the quantization pricing loop, the same pathology as drag_30k | 1 |
| fertility (ctest) | **pass** 13-27 s | 2 |
| wheel | rebuilt, delvewheel-repaired, **7/7 pytest on 3.12 and 3.13** | 1 |

Notes:
- Output quality of the successful runs, checked with `C:\tmp\check_quads.py`: pure quads, no repeated corners, every edge used by 2 quads, χ = 2 (genus 0, like the inputs):

  | output | V | F | nonquad | nonmanifold edges | χ |
  |---|---|---|---|---|---|
  | manotas_30k | 3466 | 3464 | 0 | 0 | 2 |
  | azure_30k | 2242 | 2240 | 0 | 0 | 2 |
  | GEO-head_stylized | 2434 | 2432 | 0 | 0 | 2 |
  | Azure Alien Head Scan (1) | 7239 | 7237 | 0 | 0 | 2 |
- `--threads 1` makes the run deterministic in practice and is the recommended way to get a reproducible answer; it is slower wall-clock.
- drag_30k spends **344 s of its 540 s in `compute_quantization`** (1119 of 3624 t-edges driven to x = 0 → the pricing loop; `find_negative_cycle` is Bellman-Ford over ~2 x thalf arcs, called twice per improving iteration). dragonita at 249 kV is the same situation several times larger. Any future iteration on drag/dragonita should start from a checkpoint (below), not from the CLI.
- The production numbers are single runs; given the nondeterminism, treat a FAIL as "fails at this stage" rather than as a definitive property of the mesh, except where the same port/stage repeats (gat, drag).

## 2. Result matrix (final binary, CLI defaults: decimation cap 90k, budget 300 s, stage timeout 150 s, retry on)

| mesh | result | time | holes |
|---|---|---|---|
| manotas_30k | OK 3464 quads | 125 s | 0 |
| azure_30k | OK 2240 quads | 32-79 s | 0 |
| gat_30k | OK 3236 quads | 150 s | 12 |
| drag_30k | OK 13236 quads | 218 s | 22 |
| GEO-head_stylized | OK 1095 quads (1 retry, scale 0.032) | 299 s | 0 |
| Azure Alien Head Scan (1) | OK 9209 quads | 261 s | 0 |
| manotas | OK 3701 quads | 291 s | 0 |
| gatubelo | OK 2878 quads (2 retries: param timeout -> collapse no-path) | 274 s | 0 |
| dragonita | FAIL: initial LM over the 150 s stage limit even at 54k; budget exhausted | 340 s | - |

8/9 pass, every passing mesh under 300 s. dragonita is the only hard case left (its initial LM is slower than
gatubelo's; the new `unfold` fallback for a snapped chord is in place but that run never reaches it).
## 3. The remaining blockers

### 3.1 `the seam transitions are not integer` (gatubelo, manotas at full resolution)

**Resolution-dependent, and fixed for manotas by decimating less deep.** `manotas` at the old 150k cap failed here; at **90k it passes** (3701 quads, 286 s, measured). `gatubelo` at 60k also passes the seam check (3477 quads) but takes 442 s because its **parameterization** costs 251 s at 60k (vs 6 s for Azure at 150k: `iterative_rounding` is highly data dependent) — that, not the seam, is what keeps gatubelo out of the budget. Consequences applied this session: the automatic decimation cap is now **90k** (was 150k) and each retry decimates a further 0.6x, so a retry is both coarser in `scale` and smaller in vertices. The underlying integer-transition inconsistency at 150k is still unexplained: the next step is to instrument `seamless_slim`'s gap loop (`tutte/tutte_slim.h:45-65`) with `#ifdef METRIKO_DIAG_SEAM` and print, for the worst seam edge, both copies' uv, the two tquads and whether each endpoint is a rim node of its tquad.

### 3.2 `q-face ring does not close` (drag, gat)

Symptom on both drag_30k and gat_30k: ~110 unpaired q-ports, then the ring walk collects 4 q-edges but the 4th step ends at a different port of the **same** q-vert:

```
q-face ring does not close: the edge ending at port 8819 (uv -1, -58) is followed by port 8818,
not by port 773 (uv -1, -58) that opened the ring
```

What was tried and did **not** help (revert if you see it in the tree — none of these are in the tree):
- filtering the vertex-port match by `abs(p->uv - b) < EPS` (same chart): unpaired 107 → 211.
- grouping the vertex ports per distinct corner uv (one group per chart side): manotas *broke* too, so reverted.

Best remaining hypotheses, in order:
1. `generate_vqvert_qport` builds **one port group per vertex** mixing the corners of every chart side, and accumulates the angle `acc` across them. On a seam vertex the cyclic order of the ports is not the order around the vertex, so `prev_id`/`next_id` can jump to a port of another chart. Fix: group by corner uv *and* accumulate the angle per group — note the naive version of this was tried and made things worse, so also check `generate_q_vert` (it creates a vertex-q-vertex from the *first* corner only, so a vertex integral on one side and not on another can miss its port group entirely).
2. `generate_q_edge` does not implement the paper's §5.5 line 21: *when the trace crosses into a face whose uv orientation is opposite to the previous one, invert the tracing direction*. With the current no-flip snap policy drag_30k's slim result keeps 0 flipped faces (9 stay exactly collapsed), so this is probably not the trigger, but it is a real gap in the tracer.
3. The remaining **9 collapsed faces** on drag (and a similar set on gat) are t-mesh "pockets": a triangle whose three vertices all lie on **one t-quad side**, so Tutte's prescribed integer rim coordinates make it exactly degenerate (old handoff §3.2 item 5). They come from `face_cutting` / `split_degenerate_faces` in `tutte/tutte_cutting.h`. QEx is supposed to tolerate them (paper §2.1) — the missing piece is the paper's Algorithm 1 line 13, "as a final step collapse all edges with zero length after the precision truncation": metriko never collapses uv-zero-length edges, so degenerate faces stay and the iso-line tracer gives up inside them (`pick_next_half` → dangling port). Implementing that uv-space edge collapse (union-find on equal corner uvs, rebuild faces/cfn/matching1/singular1) is the most promising single fix left.

## 4. Diagnostic harness in `C:\tmp` (rebuilt this session)

- **`metriko_cp.h`** — checkpoint save/load. `CP1` (`*_tm.bin`): hm (pos/idx/cfn), seam, matching, singular, t-nodes/t-edges/t-halfs/t-quads/fixed, unit. `CP2` (`*_uv.bin`): hm_emb (pos/idx/cfn), slim `uv`, cut-mesh `idx`, matching1, singular1.
- **`diag_dump.exe <mesh> <out-prefix> [scale] [run-rest]`** — runs to the quantization, writes CP1, then runs the rest (collapse loop … qex) and writes CP2 before sanitization. Timings per stage on stdout.
- **`diag_qex.exe <cp2.bin> [mode]`** — replays the qex stage: `0` production sanitization, `2` the experimental conditional-snap variant, `3` ring diagnostic (prints the four q-halfs), `4` vertex-port-group sizes. **Faithful**: reproduces the CLI exactly (manotas 3464 quads, drag/gat the same ring failure). Iterating here takes seconds.
- **`diag_replay.exe <cp1.bin> [rest] [mode] [roundtrip]`** — replays the collapse loop from CP1. **KNOWN BROKEN**: the reconstructed t-mesh does not reproduce the live pipeline (it throws `no path through tquad` and eventually access-violates, while the same state collapses cleanly in `diag_dump`). Something of the Emesh state is not reconstructed faithfully; do not trust it for the collapse stage. CP2 replay is unaffected.
- Build scripts `build_diag_{dump,qex,replay}.bat` (cl + the prebuilt `metriko_core.lib`; they compile the *current* headers, so header edits are picked up without touching the library).
- Checkpoints on disk: `C:\tmp\cp\{drag,manotas,gat,azure}_{tm,uv}.bin`.

## 5. Reproduce

```bat
build\dev.bat build\uvbuild.bat
uv tool run delvewheel repair wheelhouse\raw\*.whl -w wheelhouse
build\dev.bat ctest --test-dir build\cp312-abi3-win_amd64 --output-on-failure
build\cp312-abi3-win_amd64\metriko-cli.exe -i C:\tmp\manotas_30k.obj -o C:\tmp\out.obj
C:\tmp\diag_qex.exe C:\tmp\cp\drag_uv.bin 3      :: ring diagnostic in ~40 s
```

Header edits still do NOT trigger rebuilds (Spanish-locale MSVC + Ninja `/showIncludes`) — touch `src/metriko_api.cpp src/cli/main.cpp python/bindings.cpp` or build clean (see `AGENTS.md`).

## 6. Environment notes

- Windows 11, MSVC 17.14, CMake 4.3.3, Ninja, Conan 2.29.1, Python 3.12.11 / 3.13.13. `build\dev.bat` sets `OMP_NUM_THREADS=8`; the CLI runs above were started **without** it (unset) to match the earlier reference runs.
- i5-12450HX, 24 GB, thermally throttling. `compute_quantization` is single-threaded and dominates on some meshes (344 s on drag_30k).
- Wheel test venvs (uv-managed, **no pip** — use `uv pip install --python <venv>\Scripts\python.exe --force-reinstall --no-deps`): `C:\tmp\venv312w`, `C:\tmp\venv313`.

## 7. References

- Upstream: https://github.com/komietty/metriko. QEx paper (OCR markdown): `D:\onito\OneDrive\Desarrollo\meshora\Inspo\Tessra\docs\papers\QEx_Robust Quad Mesh Extraction.pdf_by_PaddleOCR-VL-1.6.md` (§2.1 degenerate triangles are legal input, §3.2 sanitization, §5.5 q-edge tracing / flip inversion).
- Key sources: `include/metriko/lib.h` (pipeline), `tmesh/emesh_collapse_{path.h,th.ipp,tq_prepare.ipp,tq_execute.ipp}`, `tmesh/emesh_validate.h`, `hmesh/{hmloc.h,hmflat.h}`, `tutte/{tutte_cutting.h,tutte_params.h,tutte_slim.h}`, `qex/{sanitization.h,gen_q_vert.h,gen_q_port.h,gen_q_edge.h,gen_q_face.h}`, `quantization/quantization.h`, `common/utilities.h`, `patching.h`.
- History: `C:\Users\jhonr\AppData\Roaming\devin\cli\summaries\history_879571d5df5f48c0.md` (session 1).

### 1.6 Stage timeout + retry (session 3, last change)

`RemeshOptions::stage_timeout_seconds` (CLI `--stage-timeout`, Python `stage_timeout`, default 150 s in both
front-ends; 0 = off in the C++ API) abandons the current attempt when a single stage runs longer than the
limit and lets the retry ladder take over with a coarser scale **and** a 0.6x decimation. It is checked inside
the long loops, so it interrupts within one iteration: the initial LM solve
(`solver/levenberg_marquardt.h`, throws with `RemeshStage::parameterization`), the quantization pricing loop
(`quantization/quantization.h`) and SLIM (`tutte/tutte_slim.h`). A timeout always earns a retry while any
budget is left, regardless of the cost estimate. Verified: with a 5 s limit, gat_30k aborts twice at the
parameterization stage and reports `parameterization: the initial integration ran past the stage timeout`;
with 30 s it completes normally at scale 0.02 in 55 s. NOTE: the 300 s budget is soft; a timeout+retry can
overshoot it, and the CLI reports the real elapsed time.

**Measured with a single attempt at 60 k decimation (no budget, no retry):** GEO 331 s / 2365 quads,
Azure 212 s / 3150, manotas 240 s / 3125, gatubelo 319 s / 3478 (passes!), at 90 k: GEO 299 s / 1095 (via one
retry), Azure 261 s / 9209, manotas 291 s / 3701. The cost is dominated by data-dependent stages (the initial
LM: 40 of 100 iterations, ~4 s each at 60 k), not by the decimation size alone -- which is why the timeout +
retry is the right mechanism and a single "good" decimation cap does not exist.

### 2.1 Second pass (repeats, same binary, machine idle)

All eight passing meshes pass again, and faster than the first pass (which ran while the wheel was building):

| mesh | run 1 | run 2 (idle) | holes run 2 |
|---|---|---|---|
| manotas_30k | 125 s | 118 s | 7 |
| azure_30k | 32-79 s | 60 s | 0 |
| gat_30k | 150 s | 97 s | 12 |
| drag_30k | 218 s | 213 s | 25 |
| GEO-head | 299 s | 124 s (1 retry) | 0 |
| Azure full | 261 s | 108 s | 14 |
| manotas | 291 s | 133 s | 14 |
| gatubelo | 274 s | 272 s | 6 |

The hole counts vary run to run (the known nondeterminism of the qex port pairing): always small (< 0.5 % of the
quads), and `--threads 1` is the way to pin them. dragonita remains the only failure (its initial LM needs more
than the 150 s stage budget even at 54k vertices).

## 8. Reproducibility (session 4)

Three sources of run-to-run divergence were removed:
1. `tutte/tutte_slim.h`: the energy used `#pragma omp parallel for reduction(+:e)` and the assembly used
   `#pragma omp atomic` into K and rhs -- both fix the floating-point sum order by the thread scheduling. The
   energy is now computed per face and summed in face order, and the assembly goes into per-thread buffers
   (capped at 8) that are merged in thread-index order.
2. `igm/iterative_rounding/iter_rounding.h`: CHOLMOD's supernodal factorization parallelizes its numeric update
   internally with a scheduling-dependent order -> `c.nthreads_max = 1` in `CholeskyWrapper`.
3. `solver/sparse_solver.h` + `CMakeLists.txt`: the same for the Cholmod wrappers used by SLIM/Tutte
   (`deterministic(solver)` helper) and `EIGEN_DONT_PARALLELIZE` for Eigen's internal OpenMP products.

**Measured state** (manotas_30k, 5 runs): identical quad count, identical vertex count and identical hole count
every run (3455 / 3466 / 21). The outputs are the same mesh up to vertex numbering and last-bit coordinates:
every vertex of run A has a counterpart in run B within 6e-7 (median) and 3.8e-3 (worst). It is NOT byte
identical yet: the remaining variability is in the vertex ordering of `weld_quad_soup` and in the local surface
projection (`patching.h:88`), which amplifies 1e-15 input differences up to its tolerance.


## 9. Hole filling (session 4): the output is watertight by construction

`patching.h::fill_holes` runs at the end of `extract_quad_mesh` (optional out-param): it walks the boundary of
the welded quad mesh (a directed quad edge without its reverse is a boundary edge), chains it into loops and
closes each one:

- a loop of **4 edges** is exactly one dropped q-face ring -> it becomes a **quad** again;
- a longer loop (several dropped rings sharing corners) -> a **triangle fan** from its centroid (the centroid
  vertex is appended to `pos`).

The fill polygon runs against the quads (a quad edge is directed a->b with the quad on the left), so the
orientation is consistent. `RemeshResult::fill_idx` carries the triangles, the CLI writes them after the quads
in the same OBJ, and Python gets them in `stats["fill_idx"]` with `return_stats=True`. `holes` stays as the
informative count of dropped rings.

**Verified on Azure full** (the run that had 14 dropped rings): output `9228` faces = `9206` quads + `22`
triangles, **boundary_edges = 0**, nonmanifold_edges = 0, euler = 2. The mesh is closed and quad-dominant
(99.76 % quads). drag_30k (22-25 rings) was not re-verified because that run hit the quantization budget (its
known variance), which is unrelated to the fill.

Pending: the dragonita LM speed-up (reuse the J sparsity pattern in the LM iterations, cheapen the injective
barrier Jacobian) and the final regression matrix with the fill in place.

## 10. LM speed-up: measured limit (session 5)

Instrumented the initial LM solve (`solver/levenberg_marquardt.h` reports its split under `--timings`, and
`CholeskyWrapper` reports analyses / copy / analyze / numeric). On gat_60k (40 iterations, 144 s):

```
LM breakdown: jacobian 21.0 s, damping 18.2 s, factorize 98.7 s, solve 4.8 s
cholmod:      1 analysis, copy 3.9 s, analyze 1.5 s, numeric 93.1 s
```

**The numeric factorization of the normal equations is the wall (65 %), and it is already optimal in the cheap
ways**: the symbolic analysis runs ONCE (the pattern is stable), so there is nothing to cache. Re-enabling
CHOLMOD's internal threading (`--fast`) makes it 76.5 s instead of 98.7 s (1.29x) but brings back run-to-run
divergence (verified: 3457 vs 3455 quads), which is why determinism is the default. The planned "reuse the J
sparsity" only addresses the jacobian+damping part (39 s, 27 %) -> about 15 % of the LM total. A real 2-3x
needs a better fill-reducing ordering (METIS inside CHOLMOD: the vendored SuiteSparse is built without it, so
CHOLMOD falls back to AMD) or a different solver; both are build-level changes.

Product outcome: `RemeshOptions::deterministic` (default true) / CLI `--fast` gives the caller the trade-off
explicitly. dragonita stays outside the 300 s budget with the deterministic configuration.

### 9.1 Manifold-safe fill (verification)

The first fill version produced 5 non-manifold edges on drag_30k: a *pinched* boundary loop (a vertex that
belongs to two holes) was filled across the pinch. The fill now walks each boundary edge once (a vertex may
have several outgoing boundary edges), rejects loops that repeat a vertex or end open, and refuses any face
that would push an edge past two incident faces. Verified closed AND manifold on two meshes with dropped rings:
gat_30k -> 3236 quads + 14 triangles, azfull -> 9206 + 22, both boundary_edges = 0, nonmanifold_edges = 0,
euler = 2. Pinched loops that are rejected leave a hole (drag_30k has these); that residual is the qex pairing
defect, not the fill.

### 9.2 Pinched-loop decomposition + the drag_30k blocker

The fill now decomposes the boundary graph into **simple cycles** (walk the unused boundary edges, emit the cycle
when a vertex repeats) instead of walking vertex by vertex, so a pinched hole (a vertex shared by two holes) is
split into fillable cycles; the manifold-safety check stays. Verified on gat_30k (3236 quads + 14 tris) and
azfull (9206 + 22): boundary_edges = 0, nonmanifold_edges = 0, euler = 2.

**drag_30k could not be re-verified**: its quantization pricing needs 126-966 s of CPU depending on the run (one
run at scale 0.03 was still going after 16 min), so most runs never reach extraction. That is the pricing
pathology, not the fill. The next real lever for drag/dragonita is therefore the pricing loop (Bellman-Ford ->
SPFA or a warm start), not the solver.

## 11. Full matrix with hole filling (session 6, final build)

CLI defaults (decimation 90k, budget 300 s with grace, stage timeout 150 s, retry on, deterministic CHOLMOD):

| mesh | result | time | dropped rings (filled) |
|---|---|---|---|
| manotas_30k | OK 3455 quads | 241 s | 21 |
| azure_30k | OK 2240 quads | 87 s | 0 |
| gat_30k | OK 3236 quads | 149 s | 12 |
| drag_30k | FAIL quantization time budget | 481 s | - |
| GEO-head_stylized | OK (1 retry) | 124 s | 0 |
| Azure Alien Head Scan (1) | OK 9206 quads | 102 s | 14 |
| manotas | OK 3701 quads | 135 s | 14 |
| gatubelo | OK 2878 quads (2 retries) | 267 s | 6 |
| dragonita | FAIL quantization time budget | 345 s | - |

7/9 pass, all under 300 s, and every passing output is watertight (verified on gat_30k and azfull with
`C:\tmp\check_mesh.py`: boundary_edges = 0, nonmanifold_edges = 0, euler = 2). **Both remaining failures are the
same root cause: the quantization pricing loop** (drag_30k and dragonita are the same mesh at two resolutions,
and their pricing time varies between 126 s and 966 s of CPU from run to run). That is now the single highest
value target, ahead of the LM solver.

## 12. SPFA pricing: the real performance fix (session 7)

`find_negative_cycle` (quantization.h) was Bellman-Ford from scratch (n full passes over the arcs, re-run twice
per improving iteration). It is now SPFA (a work queue, only the nodes whose label changed are relaxed). Same
decision, different (still deterministic) choice of cycle. The cost of the pricing collapsed:

| mesh | before | after |
|---|---|---|
| manotas_30k | 118-241 s total | **43 s** |
| azure_30k | 87 s | **19 s** |
| drag_30k | 126-966 s in quantization alone, usually over budget | **188 s total, passes** |
| gatubelo | 267 s | **261 s** |

Full matrix with SPFA (CLI defaults, deterministic): manotas_30k 43 s / 3457 quads, azure_30k 19 s / 2240,
gat_30k 91 s / 2791 (scale 0.032), drag_30k 188 s / 12609, GEO-head 123 s, Azure full 103 s / 6851,
manotas 132 s / 3700, gatubelo 261 s / 5587 (scale 0.032), dragonita **fails at 379 s** in the t-mesh collapse:
`boundary segment V(40704) -> V(40334) lies in no face (6 faces on the first, 5 on the second)` -- two nodes
snapped to non-adjacent vertices (a chord). That is a t-mesh construction issue, no longer the pricing or the LM.
8/9 pass, all under 300 s.

### 12.1 The residual hole on gatubelo: an inconsistently oriented boundary

The 8 remaining boundary edges of the gatubelo output form ONE simple 8-cycle
(2779-2802-3088-3106-2848-2800-2715-2683), all vertices of degree 2 -- a perfectly fillable hole that the fill
skipped. Cause: `fill_holes` walks the boundary as a *directed* graph (each directed quad edge without its
reverse). That only forms cycles when the incident quads are consistently oriented; a few flipped quads make
the boundary non-orientable, the walk hits a vertex with no unused outgoing edge, and the loop is dropped.
Fix for the next round: build the boundary as an **undirected** graph, take its connected components, and fan
each simple component from one centroid (triangles with each new interior edge shared by two triangles),
orienting each face against the existing face across its boundary edge. A 4-component can still become a quad.

### 12.2 The fill closes every tested mesh (pairing + cycle splitting)

The boundary is now handled as an undirected graph whose edges are **paired at every vertex** (all degrees are
even on a manifold boundary); following the pairing yields closed walks, which are then **split at their repeated
vertices** into simple cycles. Each cycle is filled (4 edges -> quad, longer -> centroid fan). This covers both
failure modes seen in practice: flipped quads (non-orientable boundary) and a vertex shared by two holes (degree
4, the gatubelo figure-eight).

Verified with the final build (boundary_edges / nonmanifold_edges / euler):

| mesh | faces | closed |
|---|---|---|
| gatubelo | 5587 quads + 34 tris | 0 / 0 / 2 |
| gat_30k | 2791 quads + 14 tris | 0 / 0 / 2 |
| drag_30k | 12609 quads (pure) | 0 / 0 / -2 (genus 2) |
| manotas_30k | 3457 quads + 18 tris | 0 / 0 / 2 |
| Azure full | 6851 quads | 0 / 0 / 2 |
| manotas full | 3700 quads + 15 tris | 0 / 0 / 2 |

The only remaining failure of the matrix is dragonita (t-mesh collapse chord, a separate construction issue).

## 13. dragonita: the actual cost profile (session 8)

Two things were fixed/tried this round:
1. `repair_crossing_tedges` no longer aborts when a tedge cannot be unfolded (a snapped chord between
   non-adjacent vertices): the reroute is skipped, like every other failed repair. Contained change, no effect
   on meshes without chords.
2. Measured dragonita with **no budget and no stage timeout** (`--stage-timeout 0 --time-budget 0`):

```
decimation 6.0 s | field 19.5 s | parameterization 412.7 s | tmesh 8.9 s | quantization 231.3 s | ...
```

So dragonita is simply far outside the 300 s budget: the initial LM alone is 413 s at the 90k decimation, and
the quantization another 231 s. It cannot fit the budget with the current solver -- the only structural fix is a
better fill-reducing ordering (METIS in CHOLMOD, which the vendored SuiteSparse is built without) or a different
solver. The collapse chord is a secondary defect behind that cost wall.

## 14. METIS was measured and is NOT the fix (session 9)

The vendored SuiteSparse already compiles **METIS 5.1.0** into CHOLMOD (`CHOLMOD_PARTITION=ON`,
`CHOLMOD/SuiteSparse_metis`, namespaced `SuiteSparse_metis_*`), so no external dependency had to be built. The
question was whether CHOLMOD uses it. Instrumented (`--timings`, `cholmod analyze: selected method N, ordering
X, lnz Z`):

| ordering | lnz | analyze | numeric | factorize total |
|---|---|---|---|---|
| AMD (CHOLMOD's trial choice) | 1.467e7 | 1.5 s | 94.4 s | 100.3 s |
| METIS (forced: `nmethods=1`, `method[0]=CHOLMOD_METIS`) | 1.398e7 | 15.3 s | 87.4 s | 107.1 s |

METIS reduces the fill by 5 % and the numeric factorization by 7 %, but its analysis costs 15 s versus 1.5 s, so
the total is **worse**. The matrix is a 2-D-like lattice where AMD is already near optimal: **the 94 s numeric
factorization is irreducible with this solver**, which closes the ordering avenue. Forcing METIS was reverted
(the diagnostic line stays, printed only under `--timings`).

Consequence: dragonita (parameterization 412.7 s, quantization 231.3 s with no budget) cannot fit 300 s without a
different solver formulation (e.g. avoiding the normal equations) or a coarser decimation policy. It stays a
documented limit.

## 15. Final matrices (session 9, delivered build)

Same binary (SPFA pricing, undirected fill with pairing/splitting, deterministic by default, METIS reverted),
run twice on the same machine at different load. **The pass/fail depends on machine load**, because the 300 s
budget is absolute:

| mesh | idle machine | loaded machine (2-2.5x slower) |
|---|---|---|
| manotas_30k | OK 43 s / 3457 quads | OK 112 s / 3455 |
| azure_30k | OK 19 s / 2240 | OK 54 s / 2240 |
| gat_30k | OK 91 s / 2791 | OK 240 s / 2791 |
| drag_30k | OK 188 s / 12609 | FAIL quantization budget (485 s) |
| GEO-head | OK 123 s (1 retry) | OK 271 s (1 retry) |
| Azure full | OK 103 s / 6851 | OK 245 s / 6851 |
| manotas full | OK 132 s / 3700 | OK 280 s / 3701 |
| gatubelo | OK 261 s / 5587 | FAIL collapse budget (372 s) |
| dragonita | FAIL 379 s (collapse chord) | FAIL parameterization budget (345 s) |

8/9 idle, 6/9 loaded. Every output verified closed and manifold (`C:\tmp\check_mesh.py`: 0 boundary edges,
0 non-manifold edges, correct euler characteristic). Wheel + delvewheel + 9/9 pytest on Python 3.12 and 3.13,
ctest green. Nothing committed.

## 16. Why the LM is 1.6M unknowns: the per-face block is 4x too big (session 10)

Read-only instrumentation (`--timings`, NaiveIntegration constructor) on gat_30k (30k vertices, 60k faces):

```
integration system: full field 60428, reduced 59998, face vars 479968, Cfull 430x60428 (nnz 1456)
LM system:          539966 unknowns, 1199922 residuals, 3407788 nnz (2.8 per row)
```

The composition of `xSize = UFull.cols() + F2.size()`:

| block | size | meaning |
|---|---|---|
| `UFull` (null space of `Cfull`) | 59,998 | ~2 per vertex: the u values, essentially minimal (Cfull only has 430 constraints) |
| `F2` (per-face field) | **479,968** | **8 per face** = 4 complex directions per face |

A 4-RoSy field needs **one** representative direction per face (2 reals) plus the 4-fold rotation applied
across seams at the end; the code carries the four directions as independent unknowns (8 reals). That is a 4x
redundancy in the dominant block, and it explains everything measured: 18 unknowns per vertex instead of ~4, a
fill (lnz) of 2.3e7 at 90k vertices, 215 s of numeric factorization for dragonita, and 540k unknowns even for a
30k mesh.

**The fix is a formulation change, not a solver tweak**: keep one complex value per face in `F2` (and one per
vertex in `UFull`), and let the 4-fold symmetry enter only through the seam transitions. That shrinks the system
~3-4x (540k -> ~180k for a 30k mesh), which should cut the factorization 4-9x. It touches `compute_jacobian`,
`G2`/`X2F`, the injective barrier and the rounding, so it must be implemented behind an option and gated by the
full matrix (8/9, watertight outputs, identical reproducibility, no quality loss) -- not merged blind.

### 16.1 The redundancy is exactly the sign-symmetry expansion (safe to remove in principle)

`parameterization.h::set_sign_symmetry(4)` sets `n = N/2 = 2` and

```
lreductor = [ I(2) ; -I(2) ]        // 4 x 2
```

so the per-face field has **n = 2 real DOF** and `F2` (8 per face) is the *expanded* representation: the four
sign-symmetry copies (+u, -u, +v, -v) of those two values. The expansion is an exact linear map, so eliminating
it is *not* an approximation: the LM can be posed on the reduced block (`n = 2` per face, `lreductor` folded into
`G2`/`X2F` and the barrier's residual), and the objective is the same up to that linear change of variables.

That makes the refactor well defined (compose `G2`, `X2F`, `compute_jacobian`, `InjectiveBarrier` with the
expander/`lreductor`), but it still touches the four places above, so it goes behind an option and through the
full matrix as the gate. Measured potential: 540k -> 180k unknowns at 30k vertices (dragonita 1.62M -> ~540k),
which is the difference between 215 s and tens of seconds in the numeric factorization.

### 16.2 Hypothesis falsified: the sign-symmetry expansion is NOT a removable redundancy

Read-only check at the end of the naive LM (gat_30k): the four sign-symmetry copies of `F2` are stored expanded
(N = 4, n = N/2 = 2) but they **diverge during the optimisation**:

```
|F2| 486.9, |lower+upper|/|F2| 0.736, |second pair+first|/|F2| 0.693
```

So the LM is using that extra freedom; imposing `F2 = [I;-I] u` would constrain the solution and change the
result. The 1.6M-unknown system is therefore **intrinsic to the formulation as written**, and the "shrink the
LM" avenue is closed (a different algorithm would be needed, not a reparameterisation).

Decimation lever, measured on dragonita (498,620 F -> 45,000 V / 90,004 F):

| stage | 90k V | 45k V |
|---|---|---|
| parameterization | 315.1 s | **68.1 s** |
| quantization | 231.3 s | **219.5 s** |

The LM scales steeply with resolution (4.6x cheaper), but the **quantization does not (219 s at both)** and
becomes the wall. dragonita therefore needs a cheaper pricing/rounding iteration count, not a smaller LM.

### 16.3 Two more levers measured and rejected (session 11)

1. **Warm-started SPFA** (reusing the label vector of the previous call of the same sign as the initial
   potentials): dragonita at 45k, quantization **219.5 s cold -> 281.1 s warm** (worse: the previous labels are
   not feasible for the new weights, so certifying "no negative cycle" takes more relaxations). Reverted.
2. **Lower decimation** for large meshes, dragonita (498,620 F):

| | 90k V | 45k V |
|---|---|---|
| parameterization | 315.1 s | 68.1 s |
| quantization | 231.3 s | 219.5 s |
| collapse | fails (chord) | budget exceeded at 141 s |

   The LM scales with resolution, the pricing does not. Both stages must fit 300 s together, so a decimation
   change alone cannot bring dragonita in.

Conclusion after three measured dead ends (METIS, reduced F2, warm start): the remaining cost is the *number of
pricing iterations*, which is data-dependent, and the LM's system size, which is intrinsic to the formulation.
Neither is reachable by a safe local change.

## 17. dragonita: full profile and classification (session 12)

One unbounded run (`--stage-timeout 0 --time-budget 0`), all three retry attempts, 2036 s total:

| attempt | scale / decimation | outcome |
|---|---|---|
| 1 | 0.02 / 90k V | collapse PASSES (the chord is gone) -> **embedding**: `no half between vertices` at 1010 s |
| 2 | 0.032 / 54k V | embedding passes -> **collapse**: `left 2 zero-width tquads: [1304, 1309]` at 1499 s |
| 3 | 0.0512 / 39k V | quantization 461.6 s -> **collapse**: `tquad 869 is self-intersected by one of its thalfs` at 2036 s |

Three *different* t-mesh pathologies, each at a different resolution, and the quantization is not monotone in
resolution (239 s at 90k, 461 s at 39k: its cost follows the pricing iteration count).

**Classification**: dragonita is out of scope with the current algorithm. It produces a degenerate T-mesh
(self-intersected tquad, zero-width tquads, segment with no half) at every resolution tried; fixing it is a
robustness programme for the T-mesh construction (guarantee no self-intersected or zero-width tquads are ever
produced), not a local fix. What this phase *did* deliver for it: the collapse chord no longer aborts (the
impossible reroute is skipped, which is a repair failure by definition), and `no half between vertices` now
reports the reversed-half probe so the next investigation starts from data.

Verified after these changes: ctest green, wheel + delvewheel rebuilt, 9/9 pytest on Python 3.12 and 3.13, and
the 8/9 matrix (mat6) re-run on this exact binary. Nothing committed.
