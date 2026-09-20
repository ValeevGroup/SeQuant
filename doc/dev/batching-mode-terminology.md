# Batched tensor evaluation: terminology reference

Canonical vocabulary for **batched (memory-blocked) tensor-network evaluation**
in SeQuant's cost model (`core/optimize/`) and eval engine (`core/eval/`). New
specs and code use these terms. The last section maps them to the legacy
`aux`-era code identifiers so older specs, commits, and comments stay readable —
those are **not** retrofitted; this table is the translation key.

Grounded in the TAPP index-role taxonomy: Brandejs, Hornblad, Valeev, Heinecke,
Hammond, Matthews, Bientinesi, *Tensor Algebra Processing Primitives (TAPP):
Towards a Standard for Tensor Operations*, arXiv:2601.07827 (2026), Sec. 2.2.

---

## Modes and index roles

- **mode** — a *positional* dimension of a tensor (the TiledArray / cuTENSOR
  convention; "axis" is not used). A mode is the thing we batch.
- **index / label** — SeQuant's `Index`: the *labeled* quantity occupying a mode.

The **role** of a mode in a binary contraction `D = A * B` (TAPP Sec. 2.2):

| role | TAPP symbol | where it appears |
|---|---|---|
| contracted | P | both operands, summed; not on the result |
| free | F_A / F_B | one operand + the result |
| Hadamard | H | both operands + the result (element-wise) |
| isolated / reduction | I_A / I_B | one operand only, summed within it |
| isolated / broadcast | I_C | result only |

**external** (our umbrella term) = a mode present on the **result** = free (F)
union Hadamard (H) [union broadcast I_C].

> **Divergence from the TAPP appendix glossary (deliberate).** The TAPP appendix
> synonyms are traps for us: it equates *External* with *Free* (too narrow — it
> drops Hadamard) and equates *Batch* with *Hadamard* (the batched-GEMM meaning,
> not ours). In SeQuant, **external = result-present (free OR Hadamard)**, and
> **batch = memory-blocking** (below), never Hadamard.

## Batch vs slice

Two views of partitioning a mode's range into blocks. Keep them distinct:

- **batch** — *work* partitioning: the loop over blocks; the *outside-the-loop*
  view. A mode carrying such a loop is a **batched mode**.
- **slice** — *data* partitioning: one block's restricted data; the
  *inside-the-loop* view. The runtime op that restricts a tensor to one block is
  a **slice** (`slice_mode`).

A batched mode is partitioned into **`nbatches`** blocks; the batch loop iterates
over them, taking one **slice** of the data per iteration. (So the count is a
batch-count — `nbatches` — an outside/work view; "slice" names the per-block data,
never the count.)

## `BatchModeType`: combine semantics

How a batched mode's per-block partials are combined follows from its role:

- **`BatchModeType::Contracted`** — the batched mode is *summed* (contracted `P`,
  or isolated-reduction `I`). Per-block partials are summands → **accumulate**
  (`add_inplace`).
- **`BatchModeType::External`** — the batched mode is *result-present* (free `F`
  or Hadamard `H`). Per-block partials are disjoint slices of the result →
  **scatter** (`write_into_slice` into a pre-sized destination).

So `BatchModeType` is exactly the accumulate-vs-scatter combine, and it maps
one-to-one onto the summed-vs-result-present split of the TAPP roles.

## Two batching regimes: node-local vs forest-level

*Where* a mode can be batched to actually bound memory follows from its role —
this is load-bearing, not a naming detail:

- **Node-local (contracted).** A contracted mode is summed at one node; batching
  its sum *there* bounds that intermediate (partials accumulate). The node-local
  batch candidates at a node are exactly the modes contracted at it →
  **`contracted_here`**.
- **Forest-level (external).** An external mode is *never summed* — it is carried
  onto the result of a node and of **every ancestor up to the root**.
  Scatter-batching it at an inner node re-materializes that node whole and does
  nothing for the (larger) parents that also carry it. To bound memory you must
  wrap the **whole subtree / forest** in an outer block-loop and never
  materialize it whole. External batching is therefore **forest-level**, not
  node-local, and gets its own forest-scope candidate concept (not
  `contracted_here`).

`contracted_here` naming a contracted-only set is thus *correct*, not a
limitation: it is the node-local candidate set, and node-local batching is
contracted-only by construction. External modes are batchable, at forest scope.

## New <-> legacy code map (translation key)

The cost model and eval engine still carry `aux`-era identifiers — fossils from
when the DF/RI auxiliary index was the only batchable mode. The mechanism is now
general (any mode role is batchable), so the names mislead. New code uses the
left column; this table lets you read pre-rename code, commits, and the older
`doc/dev/specs/` designs.

| Concept (new) | Legacy identifier |
|---|---|
| batch mode | batch axis |
| `BatchModeType {Contracted, External}` | `AxisKind {Contracted, External}` |
| `batchable_modes` (the set of modes eligible for batching) | `ctx.aux` |
| `open_modes[n]` (modes live on intermediate `n`) | `open_aux[n]` |
| `contracted_here` (node-local batch candidates) | `Acand` |
| `batched_here` (modes batched at this node) | `aprime` (cost model), `batch_axes()` (eval) |
| `batched_enclosing` (modes batched by enclosing loops) | `B` |
| `nbatches[k]` (batches mode `k` is cut into) | `nbatch[k]` |
| `batched_here()` / `set_batched_here()` | `batch_axes()` / `set_batch_axes()` |
| `external` (mode role: result-present) | `spectator` |
| `is_external_mode` | `is_spectator_axis` |
| `batch_external_modes` (policy flag) | `batch_spectator_indices` |
| `reconstruct_batched_modes` (emits per-node `batched_here`) | `reconstruct_axes` |
| `batchable_mode_list` (builds `batchable_modes`) | `batchable_index_list` |

The code's **`spectator`** (in `is_spectator_axis`, `batch_spectator_indices`, and
the `2026-07-11` design) is the same concept as our **`external`** ("open on the
result, contracted nowhere") - rename to `external`. (`external` is preferred over
`spectator`, which wrongly implies "uninvolved in any product"; see the External
section above.)

Legacy shorthands to retire in prose: `aux` (means "any batchable mode"),
`open_aux` ("modes live on an intermediate"), `Acand` ("contracted-at-node
candidates"), `axis` (use **mode**).
