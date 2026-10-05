# Dry run

A standalone utility that plans an evaluation and meters the plan, with no
numerical data: it reads the equations and the sizes they are to be planned
with from a JSON _driver_, runs the same steps a host program runs before it
evaluates anything -- optimize (factorize and batch), binarize, schedule -- and
replays the schedule through the dry-run backend (`eval/backends/dryrun`),
which prices every operation from the sizes alone. A run takes seconds, so the
sizes, the memory budget and the batch sizes can be varied freely.

Unlike `cost_analysis`, which reports symbolic `AsyCost` polynomials for
unbatched trees, this utility understands CSV (proto-indexed) tensors and
batching: its numbers are those of the dry-run meter of the host program.

## Build & run

Built as part of SeQuant (target `dry_run`, with `SEQUANT_UTILITIES=ON`):

```bash
dry_run --driver plan.json
```

| Option | Meaning |
|---|---|
| `--driver <file>` | the driver (required) |
| `--sizes <file>` | take the sizes (`spaces`, `slices`, `csv`) from this file -- another driver, or a JSON file with just those members -- instead of the driver's own |
| `--peak-threshold <x>` | override `batch.peak_threshold`, the per-term modeled peak in bytes the optimizer batches under (`inf` = no budget) |
| `--flavours` | keep the Kramers flavours apart in the op-class table |
| `--trace <file>` | write the per-op trace of the metered replay |
| `--top <n>` | rows of the op-class table (default 12) |
| `--equation-fingerprints` | print the tree fingerprint of every equation |
| `--dump-tree <name>` | print the planned tree of the named equation, one node per line with the indices it is sliced on |

The exit status is 0 on success, 1 on an error, and 2 when the driver records
what its producer's own planning made and this replay made something else.

## Output

```
equations: 8  roots: 8  metered trees: 7
eval trees: fingerprint 0xc55895f902cf6e95
dry-run meter [ordered, combined 7 trees]: peak=… GB (dense model) | flops persistent=… volatile=… | time persistent=… volatile=… | builds=…

products by class (flops summed over every execution of the replay; * = depends on t):
       flops   share     runs   max result  class
   …
```

- **eval trees: fingerprint** identifies the planned trees (structure, layout
  and phases of every node). Two plans with one fingerprint are the same plan.
- **dry-run meter** is the report of `eval::dryrun::meter`: the peak of the
  dense storage model, the operation counts and model execution costs of the
  work that is independent of the volatile leaf (_persistent_) and of the work
  that depends on it (_volatile_), and the number of values built.
- **products by class** sums the replay's contractions by class -- the
  operands' and result's labels with index ordinals and Kramers flavours
  dropped -- so the instances of one contraction over blocks and batches add
  up. `runs` counts executions, so a contraction inside a batch loop counts
  once per batch.

## Driver

A driver is written by the program whose planning is to be replayed (MPQC:
`plan_sizes:save`); every member can then be edited by hand.

| Member | Meaning |
|---|---|
| `format`, `version` | `"sequant-dry-run"`, `1` |
| `context.convention` | the index-space registry: `mpqc-sr` (MPQC's single-reference LCAO convention) |
| `context.vbs` | whether the registry has a separate virtual basis |
| `context.field` | `real` or `complex`, for every index space |
| `context.strict_braket_symmetry`, `context.first_dummy_index_ordinal` | the corresponding `sequant::Context` properties |
| `spaces` | `[{key, size, type?, qns?}]`: the extent of each index space; `type`/`qns`, when given, are checked against the registry the convention builds |
| `slices` | `{key: [extents]}`: the batch partition of each batched space (uniform `target_size` blocks when absent) |
| `csv.moments` | `{rank: [M_0, …, M_4]}`: the power means of the per-cluster domain sizes of the CSV composites of each cluster rank |
| `csv.union_composites` | whether a spin-free composite spans both Kramers flavours' members (twice the domain) |
| `optimize` | the `OptimizeOptions` as the optimizer receives them: `objective_function` (`dense_flops`, `dense_size`, `dense_space_time[_batched]`, `dense_time_space[_batched]`), `volatile_weight`, `cse_subnet`, `peak_flops_tolerance`, `footprint_weight`, `prune_outer_products`, `term_batch_axes`, `roofline` |
| `batch` | the `BatchPolicy`: `enabled`, `scheduler` (`ordered` or `forest_descent`), `peak_threshold` (`null` = none), `accumulation_factor`, `persistent_only`, `spectator_indices`, and `spaces: {key: {target_size, contracted, external}}` |
| `cache` | `enabled`, `min_repeats`, `max_footprint` (elements; 0 = no cap) |
| `binarize` | `volatile_label` (the leaf being solved for), `fold_kramers_eval_leaves`, `kramers_blind_leaf_label`, `kramers_blind_loop_fusion`, `kramers_fold_intermediates` |
| `equations` | `[{name, head, rhs, tree_fingerprint?, dry_tree_fingerprint?}]`: each block as handed to the planner, in SeQuant serialization (V1), and the fingerprints of the trees the producer built for it (the ones it evaluates; the ones its own dry run metered). A scalar-headed block (an observable) is planned but not metered |
| `expected` | what the producer's planning made: `eval_tree_fingerprint` (all blocks), and `meter` when it ran its own dry run |

A composite index of cluster rank `r` in a group of `k` composites of one
cluster is sized by `max(2, M_k)` of its rank's moments, a pair (`r >= 2`)
never by less than its RMS (`k >= 2`); `csv.union_composites` doubles it for a
spin-free composite.

> [!IMPORTANT]
> The equations are fixed in the driver. Steps their producer took before the
> planner saw them -- for MPQC, the cost-aware density-fitting split, which
> itself reads the sizes and the memory budget -- are not redone when `--sizes`
> or `--peak-threshold` change those inputs. To plan for other sizes
> faithfully, have the producer write the driver with those sizes in force
> (MPQC: `plan_sizes:load` together with `plan_sizes:save`).

## Self-check

When the driver carries `expected` and nothing is overridden, the utility
compares its trees (block by block, by fingerprint) and its meter report with
the producer's and prints `MATCH` or `MISMATCH` for each; a mismatch means this utility and the producer
no longer plan alike (or the driver was edited).
