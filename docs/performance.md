# Performance

## The numbers

Reference set, CLI defaults (300 s budget, 150 s stage timeout, retry on, deterministic), **idle machine**
(12-core laptop, OpenBLAS `HASWELL`):

| mesh | input | result | total |
|---|---|---|---|
| manotas_30k | 30k V | 3,457 quads | 43 s |
| azure_30k | 30k V | 2,240 quads | 19 s |
| gat_30k | 30k V | 2,791 quads (scale 0.032) | 91 s |
| drag_30k | 30k V | 12,609 quads (pure quads) | 188 s |
| GEO-head_stylized | 90k V (decimated) | 1,095 quads (1 retry) | 123 s |
| Azure Alien Head Scan (1) | 955k F → 90k V | 6,851 quads | 103 s |
| manotas | 999k F → 90k V | 3,700 quads | 132 s |
| gatubelo | 66k V | 5,587 quads (2 retries) | 261 s |

Same binary on a **loaded machine** (about 2-2.5× slower): manotas_30k 112 s, azure_30k 54 s, gat_30k 240 s,
GEO-head 271 s, Azure full 245 s, manotas full 280 s, while drag_30k (445 s), gatubelo (378 s) and GEO-head's
retry path crossed the 300 s budget.

**Read this as a range, not a promise.** The 300 s budget is absolute; the wall time of a given mesh depends on the
machine, the load and the thermal state. If you run in batch or on shared hardware, raise the budget (e.g. 600 s)
or expect the marginal meshes to fail and be retried with different settings.

## What dominates

Two stages account for almost everything on large meshes.

### `parameterization` — the LM solve

Measured breakdown on a 60k-vertex mesh (40 iterations, 144 s):

```
LM breakdown: jacobian 21.0 s, damping 18.2 s, factorize 98.7 s, solve 4.8 s
cholmod:      1 analysis, copy 3.9 s, analyze 1.5 s, numeric 93.1 s
```

* The **numeric factorization** of the normal equations is ~65 % and it is irreducible here: the symbolic analysis
  already runs once, and the fill-reducing ordering is already AMD's (forcing METIS gives a smaller fill but a much
  more expensive analysis, so the total is worse).
* The system has ~18 unknowns per vertex in the current formulation. Making it smaller would require a different
  formulation; this was investigated and is not a safe local change.
* Practical consequence: the cost scales steeply with the vertex count. 30k → ~20 s, 60k → ~145 s,
  90k → ~300-400 s. **`decimate.max_verts` is therefore the most effective performance knob.**

### `quantization` — the pricing loop

* Ordinary meshes: 5-40 s.
* Meshes where hundreds of t-edges are driven to zero: 100-460 s. The cost follows the number of pricing
  iterations (a data-dependent quantity), not the mesh size: the same mesh at 45k and 90k vertices took 219 s and
  231 s.
* The pricing already uses SPFA (a work-queue Bellman-Ford); a warm start across iterations was tried and is
  *slower*, and Bellman-Ford from scratch is what it replaced (3-8× slower).

## How to profile

Run the CLI with `--timings`:

```bash
metriko-cli -i input.obj -o out.obj --timings
```

You get one line per stage plus the internals:

```
[metriko] stage parameterization  412.666 s
[metriko] naive LM: 52 of 100 iterations, 301.324 s
[metriko] cholmod: 1 analyses, copy 7.9 s, analyze 2.7 s, numeric 203.8 s
[metriko] LM system: 1620038 unknowns, 3600082 residuals, 10292984 nnz (2.9 per row)
[metriko] stage quantization      231.261 s
[metriko] slim: 4 of 30 iterations, 184657 unknowns, energy 2.1e+05
[metriko] stage collapse           46.760 s
```

Reading it:

* `naive LM: N of 100 iterations` — the parameterization stopped at `N` iterations. Hitting 100 means it did not
  converge; between 30 and 60 is normal.
* `cholmod: numeric` — if this is most of the parameterization, the mesh is simply large for the solver: lower
  `decimate.max_verts`.
* `LM system: X unknowns` — for a given vertex count this is fixed by the formulation (~18 per vertex).
* `slim: N of 30` — hitting 30 is normal (the cap), not a failure.
* A stage that is *skipped* by a retry shows as two lines for the same stage in the log.

## Practical tuning

| symptom | change |
|---|---|
| parameterization over budget | `--decimate-max-verts 60000` (or 45000) |
| quantization over budget | raise `--stage-timeout`, or `--decimate-max-verts` lower, or accept the retry |
| a retry produced a coarser mesh than wanted | raise `--max-attempts`, lower `--coarsen-factor`, or fix the underlying density |
| everything is uniformly slow | check `--timings`: if `decimation` dominates, the input is huge; if `field`, the mesh is irregular |
| too slow in CI | `--decimate-max-verts 30000 --slim-max-iter 10` gives a quick, lower-quality result |

## Parallelism

* OpenMP is used in the field, SLIM, the Tutte solve and parts of the quantization; `--threads N` limits it.
* The CHOLMOD factorizations are single-threaded by default for reproducibility (`--fast` re-enables their
  threading, ~1.3× on the initial solve).
''s Jacobian and its factorization dominate. A 90k-vertex mesh needs roughly 2-4 GB; 150k vertices
  can exceed 8 GB. Lower `decimate.max_verts` on memory-constrained machines.
