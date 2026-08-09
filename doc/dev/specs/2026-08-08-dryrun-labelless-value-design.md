# DryRun label-free value design

_Date: 2026-08-08. Scope: `SeQuant/core/eval/backends/dryrun/` (cost model only). Status: design._

## Problem

A value in an eval DAG has **no intrinsic labels**. Only an *op* binds labels to
a value, and those labels are meaningful only within that op; the same value
consumed by two ops is bound to two different label sets (that is the whole point
of CSE with divergent slots). The op's einsum annotation (`lannot`/`rannot`/
`this_annot`) is the sole source of labels.

The DryRun cost-model `Result` (`ResultDryRun` / `ResultDryRunNested`) violates
this: it stores `indices_`, a *labeled* `svector<Index>`, as the value, and code
reads those labels as if they were stable across ops. But `indices_` only holds
whatever labels the value's **producer** used  --  which say nothing about how a
**consumer** binds it. Two bugs trace to this:

1. **(FIXED  --  `f5d05d3bc`)** The per-op FLOPs derived the contracted index set
   from the operands' stored `indices_` instead of the op's annotations. For any
   CSE-shared value used under a different binding  --  e.g. both legs of a
   `(g·C)(g·C)` contraction are the *same* cached value  --  the stored indices are
   the producer's labels, so `out ∪ contracted` unioned modes from different
   label contexts. On the C60 residual one op read `6.65e16`
   (`120⁴·4320·42³`) vs the correct `~1.3e13` (`120³·4320·42²`)  --  80% of the
   whole `dryrun_flops`, a spurious ~5×/80% "overcompute" that vanished once the
   op was costed from `lannot`/`rannot`/`this_annot`.

2. **(OPEN  --  this spec)** The slice widths `ExtentOverrides =
   map<Index, size_t>` are keyed by the value's stored (producer) labels. Under
   batching an op binds different labels, so the overrides do not match the modes
   the op sees. Worse, `overrides` do not only *size*  --  `mode_batches()` reads
   them to decide **how many blocks** to slice a mode into  --  so a targeted
   label-remap changes batch *structure*, not just sizes (verified: it flipped a
   split from 6 ops to a baseline's 13, backwards). This cannot be spot-fixed;
   the label-keyed value is load-bearing in several positional-but-label-keyed
   sites at once.

## Goal

Make the DryRun value **label-free**. It carries a positional **shape** (per
mode: space + composite proto structure) sufficient for `CostModel` sizing, and
carries no labels. Every per-op quantity  --  cost, slice overrides, batch counts,
result shape  --  is resolved **positionally** against the op's annotation, the only
source of labels.

## Design

1. **Positional overrides.** `ExtentOverrides`: `map<Index, size_t>` →
   `map<size_t, size_t>` (mode *position* → realized sliced extent).
   `slice_mode(mode, …)` and `mode_batches(mode, …)` already take a positional
   `mode`; they set/read overrides by position. `CostModel::memsize`/`flops`
   apply an override to the mode at that position of the index list they are
   handed  --  never by `Index` identity.

2. **Value shape, not labels.** `Result`'s `indices_` becomes a shape descriptor
   whose entries are used **only** for sizing (space extent + composite
   `inner_pow` via proto arity), never compared across ops. Start minimal: keep
   an `Index` list but treat it as shape-only (guarantee no code reads its labels
   for cross-op matching). Escalate to a dedicated label-free `Shape` type (per
   mode `{IndexSpace, svector<IndexSpace> proto}`) only if that removes remaining
   label misuse the compiler surfaces.

3. **Labels from the annotation.** `prod`/`sum` already receive
   `lannot`/`rannot`/`this_annot`. Cost and the result's shape are computed from
   these (already true for FLOPs after bug 1's fix). The result's positional
   overrides are built by mapping each operand's positional overrides through the
   annotation's mode correspondence  --  position `k` of an operand's shape and of
   its annotation are the same einsum slot.

4. **`CostModel::memsize`/`flops`** accept positional overrides (mode index →
   extent) applied to the index-list argument by position, replacing the
   `Index`-keyed `make_extent_fn` lookup.

## Compiler-driven migration

Change `ExtentOverrides` to `map<size_t, size_t>` **first**. The build errors
then enumerate every site that keyed by `Index` label  --  `make_extent_fn`,
`memsize`, `flops`, all `DryRunOps` bodies (`prod`/`sum`/`slice_mode`/
`mode_batches`/`write_into_slice`/`pre_sized_zeros_over_mode`), and the
`peak_profile.hpp` / `cost_profile.hpp` callers. Fix each to positional. Exposing
every "wrong thing" this way is the point of removing labels from the value.

Files: `backends/dryrun/cost_model_object.hpp` (`ExtentOverrides`,
`make_extent_fn`, `memsize`, `flops`), `backends/dryrun/result.hpp` (all ops),
`peak_profile.hpp` (`cell_footprint` builds an `ExtentOverrides`).

## Scope boundary

DryRun backend only. The wet backends (`tiledarray`/`tapp`/`btas`) never had
either bug: they contract via annotation-driven einsum
(`tapp_ops::contract(…)`) and store real tensors, so labels come from the op and
sizes from the data. **No numerics and no wet execution path are affected  --  this
is cost prediction only** (`eval:dry_run` / `cost_profile`).

## Testing

- Unbatched `(g·C)(g·C)` cost stays correct (C60 forest `none` avoidable = 0).
- The split test needs a **correctness** assertion, not a work-magnitude proxy.
  Measured after the flops fix, the split does **fewer** ops (6) than baseline
  (13): per-occurrence homing builds each leg cleanly, whereas the baseline
  shares one value and then re-slices it for the divergent legs. So *both* the
  old `dryrun_flops >` proxy and a naive `dryrun_n_ops >` are backwards -- they
  passed only because the pre-fix over-sizing inflated the split. Assert the
  real invariant instead: the divergent value is materialized at two distinct
  occ homes (or, minimally, the router changes the schedule vs. baseline).
  Design the exact check during implementation.
- Batched witness: aux-only and occ+aux `cost_profile` peak/flops sensible, and
  `mode_batches` counts **unchanged** for non-relabeled ops (positional ==
  prior behavior there)  --  a regression guard that the refactor only fixes the
  relabeled case.
- Full `[dryrun]` / `[eval]` / `[dryrun-costmodel]` suites green.

## Open questions

- Shape representation: `Index`-list-as-shape (minimal) vs a dedicated label-free
  `Shape`. Resolve during implementation from what the compiler-driven pass
  leaves misusing labels.
- Whether `AssembledCoverage` (already keyed by mode position) needs any change  -- 
  expected none, but confirm it composes with positional overrides.
