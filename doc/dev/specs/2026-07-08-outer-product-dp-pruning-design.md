# Outer-product (disconnected-subgraph) pruning for the single-term DP

**Goal:** Prune the single-term contraction-ordering DP so it never spends work
on *outer products* -- bipartitions whose two parts share no **contractible**
index (an index that is summed somewhere in the term). This collapses the
per-term search from ~3^n (n = number of tensors) to the number of *connected*
subnetworks, which is polynomial for the path/tree-shaped tensor networks CC
residual terms actually produce. Hadamard/partial contractions over a shared
index that is summed only later (a hyperedge) are **kept** -- they are legitimate
intermediates, not outer products. The optimal factorization of a connected
network is unchanged (loss-free); only the search space shrinks.

**Architecture:** Precompute, once per term, a pairwise **contractible-index
adjacency** over the tensors (edge `i-j` iff they share an index that is not an
output/target index), and from it the connected components (union-find) and a
per-subset connectivity table `connected[n]`. In the generic DP driver
`solve_single_term` (`SeQuant/core/optimize/cost_model.hpp`), skip disconnected
subsets and relax a bipartition `(lp, rp)` only when both parts are connected
subnetworks. In the CSE metadata precompute (`build_subnet_metadata`,
`single_term_detail.hpp`, additive objectives) skip the same disconnected
subsets, so no outer-product subnetwork is ever canonicalized. A term whose
adjacency graph has more than one component is a genuine **product** -- which
CC residual summands never are (linked-cluster theorem) -- and is recognized up
front from `connected[full]`; there the guards go inert and the term is
optimized by the original unpruned enumeration. The relax-loop prune is fully
model-independent; the only per-model touch is a uniform `connected[]` parameter
threaded through `build_context`.

**Tech stack:** SeQuant C++20, `SeQuant/core/optimize/`.

---

## Problem and motivation

The generic single-term DP driver (`cost_model.hpp:27` `solve_single_term`)
enumerates, for every subset `n` of the tensors, every bipartition `(lp, rp)`
and calls `m.relax(...)`:

```cpp
for (size_t n = 1; n < st.size(); ++n) {
  if (std::popcount(n) == 1) { st[n] = m.leaf(ctx, n); continue; }
  typename Model::State acc = m.init(ctx, n);
  for (auto&& [lp, rp] : bits::bipartitions(n))
    if (lp != 0 && rp != 0) m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);
  st[n] = std::move(acc);
  m.finalize(ctx, n, st);
}
```

For the batched objective (`PeakBatchedModel::relax`, `cost_model.hpp:694`) each
`relax` iterates over all `2^m` sliced-set contexts `B` (m = number of batchable
aux indices) and, per `B`, all subsets of `Acand` -- **before** it ever inspects
the child frontiers. So every outer-product bipartition still pays the full
`2^m * 2^|Acand|` inner cost even though it can never be part of an optimal
schedule. The number of `(subset, bipartition)` pairs is `~3^n`, and a large
fraction of them are outer products (splits across a disconnection).

### Evidence

Under the faithful C60 pVDZ-F12 regime (`test_eval_dryrun.cpp`), a **single**
13-tensor residual term optimized under the perf-first `DenseTimeSpaceBatched`
objective **did not complete in 38 minutes**. The full 55-term residual is
therefore un-sweepable, and -- more importantly -- MPQC's own end-to-end
perf-first path inherits the same wall-clock cost. Perf-first uses no CSE, so
the cost is entirely these wasted relaxes, not per-subset canonicalization.

## The pruning rule: connected-subgraph DP

A bipartition `(lp, rp)` of a subset `n` is an **outer product** iff `lp` and
`rp` share no index that is ever summed -- i.e. no *contractible* index. Merging
them is then a pure tensor product that only inflates rank and can never be part
of an optimal schedule of a connected term. Sharing an index that survives to
the result (an output/spectator index) does **not** save a merge: nothing is
combined, it is still an outer product.

Crucially, a merge that shares a contractible index which is summed only *later*
(a **hyperedge** -- an index on three-plus tensors, or on two tensors whose
result is consumed further up) is **not** an outer product. It is a legitimate
Hadamard/partial-contraction intermediate on the way to that later sum, and it
is kept. Requiring a merge to *share a contractible index* (rather than to
*sum one immediately*) is exactly what keeps these and removes the need for any
whole-network fallback (see below).

Model the tensors as a graph with an edge between two tensors iff they share a
contractible index (adjacency below). Then:

1. A subset `n` is **connected** iff the subgraph induced on its tensors is
   connected. A connected subset can always be assembled by allowed merges
   (grow it one tensor at a time; connectivity guarantees each addition shares a
   contractible index with the current group). A **disconnected** subset can
   only be formed by an outer product across its components -- never wanted.
2. The optimal contraction tree of a connected network is built entirely from
   connected subsets: its root split `(lp, rp)` has both parts connected (a
   disconnected part would itself need an internal outer product), recursively.

The rule follows:

```cpp
// precompute
connected = connected_subsets(adjacency, nt);   // vector<char>, size 2^nt
bool full_connected = connected[(1u << nt) - 1];

for (size_t n = 1; n < st.size(); ++n) {
  if (std::popcount(n) == 1) { st[n] = m.leaf(ctx, n); continue; }
  if (full_connected && !connected[n]) continue;   // never form a disconnected subset
  typename Model::State acc = m.init(ctx, n);
  for (auto&& [lp, rp] : bits::bipartitions(n))
    if (lp != 0 && rp != 0 &&
        (!full_connected || (connected[lp] && connected[rp])))
      m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);
  st[n] = std::move(acc);
  m.finalize(ctx, n, st);
}
```

When the full network is connected (`full_connected`, the only case that arises
for CC -- see "Product terms" below), disconnected subsets are skipped entirely
(no `init`/`relax`/`finalize` -- they pay nothing), and a connected subset is
relaxed only through bipartitions whose **both** children are connected, so its
skipped disconnected neighbours are never read as children. When the full
network is *not* connected, the guards are inert and the loop is the original
unpruned enumeration.

### Adjacency: share a contractible (non-target) index

`contractible_adjacency(network, tidxs)`: tensors `i` and `j` are adjacent iff
they share at least one index that is **not** an output/target index of the term
(`x` with `x not in tidxs`). Such an `x` is summed somewhere in the term -- now,
if it is private to `i,j`, or later, when the last of its carriers is contracted
in.

This is a genuine **pairwise, context-independent** graph. The target list is
consulted only to *classify each index* (contractible vs. output) -- a property
of the index, not a per-merge context.

- A **hyperedge** -- a contractible index on three-plus tensors -- makes all its
  carriers mutually adjacent (a clique). This is what lets the DP form the
  necessary Hadamard intermediates. Example: `A{i;;r} B{j;;r} C{k;;r}` with
  targets `{i,j,k}` and `r` contractible. `r` links `A,B,C` into one component;
  the DP is free to form `AB{i,j,r}` (a Hadamard over `r`) and sum `r` when `C`
  joins. No merge here is a disconnected outer product.
- A **spectator / output** index -- an external occupied index carried only as a
  composite protoindex, an aux (mu-tilde / K) index that survives to the output
  -- is in `tidxs`, so it is **not** contractible and creates **no** edge. Two
  branches that share *only* the external occ `i,j` are not adjacent: their
  occ-Hadamard product is never formed early. ("Sharing batched/spectator
  indices only is still an outer product.")

Adjacency compares indices at the granularity they contract, so two tensors
carrying different composites with the same protoindices (`a<i,j>` vs `b<i,j>`)
are **not** adjacent through the shared occupied protos -- those protos are
outputs, never summed.

### Precompute: index degree + target flag

One pass over all tensor slots yields, per index: its **degree** (number of
slots it occupies) and whether it **is a target** (`in tidxs`). From these:

- **adjacency / classification:** an index is contractible iff it is not a
  target (and degree >= 2, i.e. actually shared). `contractible_adjacency` and
  the component union-find are built directly from the per-index carrier lists.
- **open-set membership in O(1):** index `x` is *excluded from* `open[n]` (i.e.
  summed / internal to `n`) iff `x` is not a target **and** all `degree(x)` of
  its slots lie in `n`. The memorized degree is exactly what decides "is it
  internal yet?" for hyperedges without rescanning; targets are always open.

`connected_subsets(adjacency, nt)` then fills `connected[n]` for all `2^nt`
subsets of a component by a BFS/union-find over the induced subgraph. Cost is
`O(2^nt * nt^2)` (n=13: ~1.4M boolean ops), negligible against the DP it guards.
It is computed **once per component** and shared by every consumer (it depends
only on `network` and `tidxs`): built in the `run_single_term_opt` /
`run_single_term_opt_axes` entry points and handed to both `build_context` (CSE
metadata) and `solve_single_term` (relax loop).

### Product terms: detected up front, pruning simply disabled

Whether a term is a single connected network or a product of several components
is read off `connected[(1<<nt)-1]` -- i.e. up front, not by filling the `2^nt`
table and discovering a dead end at the top.

- **One component** (every CC residual summand -- the linked-cluster theorem
  guarantees CC residual diagrams are connected -- and every case the perf-first
  blowup comes from): `full_connected == true`, the pruned DP runs.
- **More than one component** (a genuine product; does not arise for CC
  residuals): `full_connected == false`, and the guards above are inert, so the
  driver and the CSE metadata both fall back to the original unpruned
  enumeration for that term. It is correct and, because such terms do not occur
  in the target workload, its cost is irrelevant.

Note this is *not* the over-pruning fallback we started from: the contractible-
index adjacency already keeps every Hadamard/hyperedge merge, so connected terms
(including hyperedge terms like `A{i;;r}B{j;;r}C{k;;r}`) are fully pruned with
**no** fallback. Only genuinely disconnected products -- which CC never produces
-- take the unpruned path, and they are recognized in `O(nt^2)` before any DP
work. Optimizing each component separately (and combining by outer products) is
a possible future refinement, deferred as unnecessary for the target workload.

### CSE metadata skips disconnected subsets too

The CSE-enabled additive objectives precompute, per subset of size >= 2, a
canonical-subnet id via `build_subnet_metadata` (`single_term_detail.hpp:445`),
which builds and bliss-canonicalizes a `TensorNetwork` for **every** such
subset. A disconnected subset is an outer product the pruned DP never forms, so
its canonical cost is **never looked up** -- canonicalizing it is pure waste
(a bliss-graph build is not cheap). The same pruning applies here:

```cpp
for (size_t n = 0; n < results.size(); ++n) {
  if (std::popcount(n) < 2) continue;
  if (!connected[n]) continue;            // outer-product subset: never an intermediate
  ... canonicalize, assign meta id ...
}
```

Disconnected subsets keep the `numeric_limits<size_t>::max()` sentinel `meta_id`
(exactly as singletons do today). `finalize` (`cost_model.hpp:204`) reads
`ctx.meta_ids[n]` only for subsets the driver actually DPs -- now connected-only
-- so no sentinel id is ever dereferenced.

## Correctness

**Loss-free for connected components (under the standard property).** Every
schedule the pruned DP removes contains a *disconnected* intermediate -- an
outer product `A(x)B` between two groups sharing no contractible index, which
materializes the full index union of both. For the dense flop / footprint / peak
objectives here, forming such an intermediate strictly increases both the
intermediate size and the subsequent contraction's flops relative to a connected
split whenever one exists. Under this *no-beneficial-outer-product* property, the
unpruned optimum of a connected component is built entirely from connected
subsets through connected-child splits (statement (2)), so
`pruned-optimal == unpruned-optimal` -- identical factorization *and* cost -- for
every objective. Hadamard/partial-contraction intermediates over hyperedges are
**not** removed (they share a contractible index), so nothing legitimate is lost.

This is the same heuristic `opt_einsum`/`netcon` apply by default. It is not
*unconditionally* identical to exhaustive search: there exist pathological
networks where an outer-product intermediate lowers total cost. Those do not
arise for the CC/DF-PNO terms this optimizer targets, and the **parity test
(below) is the empirical guarantee**: it asserts pruned == unpruned on the real
terms we rely on. A parity failure would flag a genuine outer-product-optimal
case for reconsideration rather than silently mis-optimizing.

**Product terms are handled by the unpruned path.** A multi-component term
(`full_connected == false`) disables the guards, so it is optimized exactly as
today -- no pruning, no change in result. Correctness there is trivially the
status quo.

**Scope of the change.** The prune is in the shared driver, so it applies to
every objective (`DenseFLOPs`, `DenseSize`, `DenseSpaceTime`,
`DenseSpaceTimeBatched`, `DenseTimeSpaceBatched`) and both the sequence-only
(`run_single_term_opt`) and axes-reporting (`run_single_term_opt_axes`) entry
points, which share `solve_single_term`.

## Testing

1. **Parity (the primary gate).** For a battery of existing multi-tensor test
   terms (the small optimize-test networks and a handful of real CC terms),
   assert the pruned DP returns the byte-identical optimal `EvalSequence` **and**
   modelled cost as the unpruned DP, for **each** objective. Implemented by
   running the DP twice behind a temporary "disable pruning" switch (test-only)
   and comparing. This is the correctness net for the whole change.
2. **Adjacency / connectivity unit tests.** `contractible_adjacency` /
   `connected_subsets` on hand-built networks: a path; a star; the
   `A{i;;r}B{j;;r}C{k;;r}` hyperedge (must be **one** component and the DP must
   form the `AB{i,j,r}` Hadamard intermediate); a term sharing only a spectator
   (output) index (must **not** create an edge); a genuine two-component product.
3. **Product-term (multi-component) safety net.** A deliberately disconnected
   two-component term optimizes to the same result as before the change
   (`full_connected == false` disables the guards). Confirms the unpruned path is
   intact and `connected[full]` correctly detects the product.
4. **Performance.** The 13-tensor C60 term that did not finish in 38 min under
   perf-first now optimizes in seconds (assert an upper bound on wall time, or
   simply that it completes within the test's normal budget).

## Non-goals / future work

- No change to any objective's cost formula, selection policy, or the batched
  sliced-set (`B`/`Acand`) machinery.
- **Per-component optimization of product terms.** Multi-component terms disable
  pruning (unpruned DP). Optimizing each component separately and combining by
  outer products is a possible future refinement, deferred because CC residual
  summands are always single-component (linked-cluster) so it would never run.
