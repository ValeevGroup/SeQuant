# Cost analysis

A standalone utility that estimates the cost of evaluating a tensor equation
without running it. It reads a serialized equation, optimizes and binarizes it,
then writes a Markdown report on the largest intermediates, the most expensive
contractions, peak storage, total FLOPs, and (optionally) cache reuse. All
figures are symbolic `AsyCost` polynomials in the index-space sizes, evaluated
at the sizes given in the driver.

The optimize/binarize step mirrors the choices MPQC's `SeQuantEngine` makes
(symmetrizer stripping, per-summand factorization, volatile-leaf weighting,
cache gating), so the estimates track what a real MPQC evaluation would do.

## Build & run

Built as part of SeQuant (target `cost_analysis`). It takes a single argument, a
[JSON](https://www.json.org/) _driver_:

```bash
cost_analysis --driver examples/ccsd_r2.json
```

Paths inside the driver (`equation_file`, `output.path`) are resolved relative
to the driver's own directory.

## Driver

| Block | Key | Meaning | Default |
|---|---|---|---|
| `context` | `spbasis` | `spinor` or `spinfree` | `spinor` |
| | `field` | `real` or `complex` | `complex` |
| | `convention` | registry: `min_sr`/`sr`/`mr`/`f12` | `sr` |
| | `aux` | factorization spaces to register: `["df"]` (Κ) and/or `["thc"]` (L) | none |
| `sizes` | `<label>` | approximate size per index space, e.g. `{"i": 10, "a": 38}` | — |
| `optimize` | `objective` | `dense_flops` / `dense_size` / `dense_peak_size` | `dense_flops` |
| | `volatile_leaf` | amplitude label being solved (e.g. `t`, `R`); empty disables volatile weighting | `R` |
| | `reorder`, `cse_subnet`, `volatile_weight`, `machine_balance`, `fast_mem_elems` | forwarded to `sequant::OptimizeOptions` | see `OptSpec` |
| `cache` | `enabled` | run the cache-reuse simulation | `true` |
| | `min_repeats` | reuse threshold to cache an intermediate | `2` |
| `output` | `path` | report filename | `cost_analysis.md` |
| | `top_n` | rows in the largest/expensive tables | `20` |
| | `dump_tree` | also write each result's binarized tree to `<name>.tree.txt` | `false` |
| `results` | `[{name, equation_file}]` | equations to analyze | — |

Each `equation_file` holds one `<head> = <rhs>` in SeQuant serialization V1, e.g.

```
R{a1,a2;i1,i2} = g{a1,a2;i1,i2} + g{a3,a4;i1,i2} t{i3,i4;a3,a4} R{a1,a2;i3,i4}
```

## Workflow

Per result equation (`process()` in `cost_analysis.cpp`):

1. **Deserialize** the `<head> = <rhs>` string into a `ResultExpr`.
2. **Strip** a leading (anti)symmetrizer — it is re-expressed by the head, not
   evaluated.
3. **Flatten** nested sums/products (mirrors `SeQuantEngine`; a no-op for
   already-flat serialized input).
4. Per summand: **`optimize`** (choose a binary contraction order for the
   objective) then **`binarize`** into an evaluation tree — *once*; the tree is
   the single source of truth for every reported number.
5. **Catalog** the tree's internal nodes, keyed by structural hash/equality so
   equal intermediates across terms collapse into one entry (with a `uses`
   count). Record each node's memory, local FLOPs, and O/V/X space signature.

A final **cache simulation** (`cache_manager`) runs over all results' trees to
count intermediates that would be cached / persist across iterations.

## Report

For each result: a one-line summary (terms, distinct/reused intermediates,
largest, peak storage), total FLOPs (symbolic and evaluated), then **Largest
intermediates**, **Most expensive contractions**, and a **Shape census** grouped
by O/V/X signature. A trailing **Cache** table summarizes the simulation.

## Tests

Two reference tests (`ctest -R cost_analysis_`) run the tool on `examples/` and
diff the report against a frozen `*.md.expected`:

- `cost_analysis_ccsd_r2` — spin-orbital CCSD R2 (a `Sum` of terms).
- `cost_analysis_df_r1` — a single density-fitted product (the non-`Sum` path).
