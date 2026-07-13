# Cost analysis

A standalone utility that estimates the cost of evaluating a tensor equation
without running it. It reads a serialized equation, optimizes and binarizes it,
then writes a Markdown report on the largest intermediates, the most expensive
contractions, peak storage, total operation count, and (optionally) cache reuse.
Figures are symbolic `AsyCost` polynomials in the index-space sizes; memory and
storage are additionally evaluated to megabytes at the driver's sizes (element
width follows the field: 8 bytes real, 16 bytes complex). Operation counts are
reported symbolically only.

The optimize/binarize step mirrors the choices MPQC's `SeQuantEngine` makes
(symmetrizer stripping, per-summand factorization, volatile-leaf weighting,
cache gating), so the estimates track what a real MPQC evaluation would do.

## Build & run

Built as part of SeQuant (target `cost_analysis`). It takes a single argument, a
[JSON](https://www.json.org/) _driver_:

```bash
cost_analysis --driver examples/ccsd_r2.json
```

Paths inside the driver (`equation_file`, `output.path`) — and the `dump_tree`
`<name>.tree.txt` files — are resolved/written relative to the driver's own
directory, not the invocation directory. The final "Report written to …" line
prints the absolute path so the output is easy to locate.

## Driver

| Block | Key | Meaning | Default |
|---|---|---|---|
| `context` | `spbasis` | `spinor` or `spinfree` | `spinor` |
| | `field` | `real` or `complex` | `real` |
| | `convention` | registry: `min_sr`/`sr`/`mr`/`f12` | `sr` |
| | `aux` | factorization spaces to register: `["df"]` (Κ) and/or `["thc"]` (L) | none |
| `sizes` | `<label>` | approximate size per index space, e.g. `{"i": 10, "a": 38}` | — |
| `optimize` | `objective` | `dense_flops` / `dense_size` / `dense_peak_size` | `dense_flops` |
| | `volatile_leaf` | amplitude label being solved (e.g. `t`, `R`); empty disables volatile weighting | `R` |
| | `reorder`, `cse_subnet`, `volatile_weight`, `machine_balance`, `fast_mem_elems` | forwarded to `sequant::OptimizeOptions` | see `OptSpec` |
| `cache` | `enabled` | run the cache-reuse simulation | `true` |
| | `min_repeats` | reuse threshold to cache an intermediate | `2` |
| | `max_footprint` | cache budget, in **elements** (not MB); `0` = unlimited | `0` |
| `output` | `path` | report filename | `cost_analysis.md` |
| | `top_n` | rows in the largest/expensive tables | `20` |
| | `dump_tree` | also write each result's binarized tree to `<name>.tree.txt` | `false` |
| `equations` | `[{name, equation_file}]` | equations to analyze | — |

> [!IMPORTANT]
> Every index space that appears in an equation should be given a size. Any
> space omitted from `sizes` falls back to its standard-registry default, so
> the reported MB/operation figures stay well-defined but reflect that default
> rather than your problem — size all used spaces (including any `aux` space) to
> get meaningful numbers.

Each `equation_file` holds one `<head> = <rhs>` in SeQuant serialization V1, e.g.

```
R{a1,a2;i1,i2} = g{a1,a2;i1,i2} + g{a3,a4;i1,i2} t{i3,i4;a3,a4} R{a1,a2;i3,i4}
```

> [!NOTE]
> CSV/PNO expressions are not supported: `sequant::AsyCost` is not proto-index aware, so
> it would size such an index as a bare one and report wrong costs.

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
   count). Record each node's memory, local operation count, and O/V/X space
   signature.

A final **cache simulation** (`cache_manager`) runs over all results' trees to
count intermediates that would be cached / persist across iterations.

## Report

One `##` section per result — a one-line summary and the **Largest
intermediates**, **Most expensive contractions**, and **Shape census** tables —
followed by a single **Cache** section pooled across all results. Memory is the
tensor's element count (product of its index-space sizes) times the field's
element width (8 B real / 16 B complex), in MB; operation counts stay symbolic.

The quantities that aren't self-evident:
- **distinct intermediates** — structurally-unique internal contraction nodes;
  equal intermediates across terms collapse to one, scalars excluded.
- **peak storage** — the transient working set of the *heaviest single
  contraction* (its two operands plus result, via `min_storage`), maxed over the
  terms. **Not** the cache footprint — the **Cache** table covers that.
- **Shape census `Size`** — the largest instance of each O/V/X signature.
- **Cache** — a simulation over the volatile-gated `cache_manager` (as MPQC
  builds it), pooled across all results. **Cached** / **Persistent** count
  intermediates reused ≥ `min_repeats` (volatile-leaf gated) / those that persist
  across iterations; **footprints** are their summed element counts in MB — a
  cache-content total, distinct from the per-result transient **peak storage**.

## Tests

Two reference tests (`ctest -R cost_analysis_`) run the tool on `examples/` and
diff the report against a frozen `*.md.expected`:

- `cost_analysis_ccsd_r2` — spin-orbital CCSD R2 (a `Sum` of terms) whose two
  ladder terms share one `g*t` intermediate, exercising the reuse census.
- `cost_analysis_df_r1` — a single density-fitted product (the non-`Sum` path).
