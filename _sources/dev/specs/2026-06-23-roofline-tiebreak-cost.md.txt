# Roofline secondary-cost model for the peak-objective tie-break

_Status:_ design of record (supersedes the additive `footprint_weight` prototype)
_Date:_ 2026-06-23
_Scope:_ SeQuant single-term optimizer cost model (`DensePeakSize`,
`DensePeakSizeBatched`) and its MPQC wiring (`SeQuantEngine`). Single-node
target; the distributed case is the same formula with the binding level shifted.

## 1. Problem

The peak objectives minimize peak memory (primary axis) and break ties among
schedules within `peak_flops_tolerance` of the minimum peak using a **secondary
cost**. Today that secondary cost is **flops only**, replay-weighted:

```
cost = Σ_ops  w_op · flops_op,    w_op = volatile_weight if op's subtree contains t, else 1
```

This under-ranks **bandwidth-bound** contractions: a memory-traffic-heavy op
(e.g. a single-PNO-index contraction, intensity `≈ 2d`) has modest flops but
large wall time, so the flop tie-break treats folding the amplitude in early
(making such an intermediate volatile and rebuilding it every iteration) as
nearly free. Consequence: the fold-t/defer-t crossover sits too high in
`volatile_weight` (≈26 for water_14), and below it the optimizer picks the
schedule that rebuilds an expensive bandwidth-bound intermediate every
iteration. The audit established this is **not** a cache bug (the eval cache's
NV/V-frontier rule caches every persistent-feeding-volatile node correctly,
0 frontier-misses); it is a cost-model deficiency.

A prototype added an additive `footprint_weight·footprint` term. It works on
the reproduction but is wrong in three ways: (a) additive double-charges
compute-bound ops that never pay the bandwidth (so it is not inert in the dense
case), (b) it charged the result footprint only, missing contraction-type ops
whose cost is operand reads, and (c) it dropped cluster/proto indices in the
hand-analysis, mis-sizing intermediates. This spec replaces it with a roofline
model.

## 2. Model

### 2.1 Per-op roofline cost

For one binary contraction at subset `n` with children `lp`, `rp`, work in
**elements** throughout (footprints are element counts; flops are FMA counts):

```
flops    = flops_of(idx[lp], idx[rp], idx[n])      // proto-aware (inner_pow), already computed
traffic  = S[lp] + S[rp] + S[n]                     // compulsory single-pass operand+result
Q        = max( traffic , kappa · flops / sqrt(M) ) // add finite-cache re-reads (Hong-Kung)
cost_op  = w · max( flops , beta · Q )              // roofline; w = volatile_weight or 1
```

- `flops`, `S[·]` are proto-explicit / block-sparse (cluster indices and their
  proto-indices are real extents; the existing `footprint_counter` / `inner_pow`
  machinery already sizes them -- only hand-analysis ever dropped them).
- `traffic = S[lp]+S[rp]+S[n]` is exactly the `both` quantity the peak DP
  already computes in `relax`.
- `w` (replay weight) wraps the whole `max`: a volatile op pays
  `max(flops, beta·Q)` on **every** replay. This is the per-iteration-rebuild
  physics and is unchanged from today.

### 2.2 The `sqrt(M)` term (finite-cache reuse)

`Q` is the data moved across the binding fast↔slow boundary. The compulsory
term `traffic` is the single-pass lower bound (each operand/result element once).
When the working set exceeds the fast level, a blocked contraction re-reads
operands; the I/O-optimal bound (Hong & Kung 1981; Loomis-Whitney) for a
matmul-shaped op is `Q = Θ(flops/√M)`, i.e. each loaded element is reused at
most `√M`-deep. Holding three tiles resident gives `b ≈ √(M/c₀)`, `c₀ ≈ 3`; the
exact constant is BLAS-blocking-dependent, so:

```
M_eff = M / c0                       // c0 ~ 3 nominal (A,B,C tiles); calibratable
Q_block = kappa · flops / sqrt(M_eff)
```

`kappa` (default 1) folds the FMA factor, packing overhead, and BLAS
inefficiency; `M_eff` folds the tile-count constant. Both are absorbed by
calibration (§4). Writing `Q = max(traffic, Q_block)` self-selects the regime:

- `I_arith = flops/traffic ≈ 2·min(index groups)`. If `I_arith < √M_eff`
  (low arithmetic reuse -- the single-PNO contraction), `Q_block < traffic` so
  `Q = traffic`: charged its true single-pass traffic.
- If `I_arith > √M_eff` (high reuse but cache-limited), `Q = Q_block`: charged
  the blocking-limited re-reads.

The effective reuse is `min(I_arith, √M_eff)` -- the small-dimension cap and the
finite-cache cap, as one `max`.

### 2.3 Regime behavior (why `max`, not additive)

`cost_op = max(flops, beta·Q)`:

- **Compute-bound** (`flops > beta·Q`, dense CC: all index groups large and
  `√M_eff ≥ beta`): `cost = flops`. `beta`, `M` are **inert** -- shape-independent,
  identical to today's flop tie-break. This is the dense-case invariant, now a
  property of the formula rather than a flag.
- **Bandwidth-bound** (`beta·Q > flops`, single-PNO contractions): `cost = beta·Q`,
  charged ×`w`. This is what tips fold-t → defer-t.
- **Capacity-bound** (large op, poor blocking, `√M_eff < I_arith`): the
  `Q_block` branch binds -- the case the single-pass roofline omits.

### 2.4 Parameters and the level abstraction

Two physical numbers per binding level:

- `beta` -- machine balance, **flops per element of traffic**: `beta = 8·F/B`
  (8 bytes/double). `O(100-500)` for CPU.
- `M` -- capacity of the binding fast level, **elements** (LLC for single node:
  `LLC_bytes/8`).

Plus two calibration constants `kappa` (≈1) and `c0` (≈3) on the `Q_block` term.

**Distributed = same formula, level shifted** (per agreement): instantiate with
the network level -- DRAM becomes the "fast" level (`M` = per-node resident
elements), the interconnect the "slow" level (`B → B_net`, `F → F_agg`,
`beta_net = 8·F_agg/B_net`, much larger). No separate code path; only the
`(beta, M, kappa, c0)` tuple changes.

## 3. Implementation

### 3.1 `SeQuant/core/optimize/cost_model.hpp`

`PeakModel` and `PeakBatchedModel` gain members (replacing the prototype's
`footprint_weight`):

```cpp
double machine_balance = 0.0;   // beta; 0 = roofline OFF (pure flop tie-break, today's behavior)
double fast_mem_elems  = 0.0;   // M (elements); 0 or beta==0 disables the sqrt(M) term
double block_tiles     = 3.0;   // c0
double block_prefactor = 1.0;   // kappa
```

Replace, in both models' `relax()`, the secondary-cost term

```cpp
double const w = (ctx.volatile_mask & n) ? volatile_weight : 1.0;
double const cflops = w * ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]);
```

with a helper `roofline_cost(flops, traffic)`:

```cpp
double const flops   = ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]);
double const traffic = ctx.S[lp] + ctx.S[rp] + ctx.S[n];   // PeakModel
//        (PeakBatchedModel: full/unsliced sizes, ctx.sz(lp,0)+ctx.sz(rp,0)+ctx.sz(n,0))
double Q = traffic;
if (machine_balance > 0.0 && fast_mem_elems > 0.0)
  Q = std::max(traffic, block_prefactor * flops /
                            std::sqrt(fast_mem_elems / block_tiles));
double const op = (machine_balance > 0.0)
                      ? std::max(flops, machine_balance * Q)
                      : flops;                              // beta==0 -> today's flop-only
double const cflops = w * op;
```

- `machine_balance == 0` reproduces the current flop tie-break exactly (safe
  default; keeps the oracle/peak unit tests untouched -- the **peak axis is
  unchanged**, only the secondary cost differs, and only when `beta>0`).
- `w` and the Pareto-frontier machinery (`pareto_insert`, `reconstruct`,
  `peak_flops_tolerance`) are unchanged; only the per-op secondary cost changes.
- `leaf()` stays 0 (a leaf read is counted in its consuming contraction's
  `traffic` via `S[lp]`/`S[rp]`; charging it again would double-count).

**Batched model -- full vs sliced.** The secondary cost is total work over a
replay, so use **full (unsliced)** `flops` and `traffic` (`ctx.sz(·,0)`), not the
per-slice sizes. Slicing reduces *peak* (the primary axis), not the total
traffic/flops of a replay. (Open refinement §5.)

### 3.2 `SeQuant/core/optimize/single_term.hpp`

`single_term_opt` already receives `footprint_weight`; rename/extend its peak
arms to forward `machine_balance`, `fast_mem_elems`, `block_tiles`,
`block_prefactor` into the two model constructors (drop the `(void)`s). The
`DenseFLOPs`/`DenseSize` arms keep the existing additive `footprint_weight`
(it remains meaningful there as a one-time footprint penalty; orthogonal).

### 3.3 `SeQuant/core/optimize/options.hpp`

`OptimizeOptions` gains `machine_balance`, `fast_mem_elems`, `block_tiles`,
`block_prefactor` (defaults `0,0,3,1` → roofline off). `opt_pure_product`
(`optimize.cpp`) forwards them to the two peak arms (it already forwards
`footprint_weight`, `peak_flops_tolerance`).

### 3.4 MPQC `SeQuantEngine` (`sequant_engine.{h,cpp}`)

- `SeQuantEngineDefaults` + members: `machine_balance_`, `fast_mem_elems_`
  (and optionally `block_*`). Read `sequant:optimize:machine_balance`,
  `sequant:optimize:fast_mem_elems` (non-negative).
- In `make_optimize_options`, set `opts.machine_balance`/`opts.fast_mem_elems`
  **only when `cache_` is on and a volatile leaf exists** -- same gating
  rationale as `volatile_weight` (the per-replay rebuild model is meaningful
  only when persistent intermediates are cached).
- Remove the prototype's `footprint_weight_` engine wiring.

### 3.5 Cleanup

Revert the `CacheClassify` diagnostic in `cache_manager.hpp` (it served its
audit purpose: proved 0 frontier-misses) before committing.

## 4. Calibration (optional, recommended)

`beta` and the `Q_block` prefactor are where the model is sensitive. Defaults
(`beta≈200`, `M≈8 MB/8 = 10⁶ elem`, `c0=3`, `kappa=1`) put the dense/PNO
crossover in a wide basin and should suffice. A one-time startup
micro-calibration -- time 2-3 representative contractions (one compute-bound,
one bandwidth-bound) and fit `(beta, kappa·/√(M/c0))` -- removes the guesswork
and adapts to the machine/BLAS. Out of scope for the first cut; defaults +
keyval overrides ship first.

## 5. Open questions / refinements

- **Batched full-vs-sliced.** Using full sizes ignores that slicing can make a
  per-slice working set fit cache (raising `√M_eff` effectively). First cut:
  full sizes. Refinement: compute `Q_block` against the per-slice working set.
- **`beta`/`kappa` degeneracy.** They are degenerate in the `Q_block` branch
  (`beta·kappa·flops/√M`) but not in the compulsory branch (`beta·traffic`),
  so both are kept; calibration fits the product where it matters.
- **Per-op-class `beta`.** Batched-Summa Hadamard ops lean more bandwidth-bound
  than GEMMs; a single scalar captures the dominant split, an op-class `beta` is
  a later refinement.
- **Interaction with `peak_flops_tolerance`.** The ε-band now selects min
  *roofline cost* within ε of min peak -- the intended secondary objective
  (min wall-time among near-peak-optimal). No change to the band mechanism.

## 6. Validation plan

1. `machine_balance=0` ⇒ bit-identical to current behavior; full SeQuant unit
   suite green, peak/oracle tests unaffected.
2. **water_8, vw=10** (`tpno=1e-6`): with a physical `beta` (not `fw=10`),
   the schedule reaches the defer-t `[6,3,3,3]` μ̃-chain pattern (matches vw=100);
   energy unchanged to ~1e-10.
3. **Dense-case invariant:** a dense (non-CSV) CCSD term's factorization is
   **unchanged** vs `machine_balance=0` (the `max` collapses to flops).
4. **14-mer** (cluster): vw=20 + physical `beta` drops the 34 GB-chain rebuild
   from every-iteration to once (the original payoff).
5. Sweep `beta` to confirm a wide insensitive basin (robust default).
