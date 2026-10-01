# Extended Wick theorem (generalized normal order) — design

Date: 2026-10-01
Status: approved in discussion, awaiting implementation plan

## Goal

Support the extended (Mukherjee/Kutzelnigg) Wick theorem relative to a general
multiconfigurational reference. Products of generalized-normal-ordered (GNO)
operator strings are reduced to sums of GNO strings multiplied by one-body
densities (γ), one-body hole densities (η) and density cumulants (κₖ, k ≥ 2),
all confined to the active (partially occupied) space.

Two uses, in delivery order:

1. **Reference expectation values** (full contractions): scalar expressions in
   γ, η, κₖ. Replaces the `ref_av` core-vacuum workaround in `mbpt`.
2. **Operator algebra in GNO** (partial contractions): products and
   commutators of GNO operators returning GNO operators, with cumulant-rank
   truncation. Needed for DSRG, canonical transformation, MR-UCC and their
   BCH expansions.

The design serves 2; milestone 1 is a filter on it.

## Decisions taken

- Cumulants of any rank, bounded by a user-set `max_cumulant_rank`
  (default: unbounded).
- Spin-orbital only. Two explicit hooks are left where spin-free would plug
  in; spin-free itself is a separate, later spec. Spin-free results meanwhile
  come from `spintrace` on spin-orbital output.
- The reference is `Vacuum::MultiProduct`, coexisting with the existing
  `SingleProduct` path, which is left untouched.
- **Implemented as post-processing of the standard Wick theorem**, not as a
  new engine. The engine changes are confined to the quasiparticle
  classifiers, the middle factor of `contract`, and skipping the pair-based
  connectivity filters. The recursion, topology and `reduce` logic are not
  modified.
- Nonorthogonal active orbitals are out of scope: spaces carrying protoindices
  must be disjoint from the active space.

## Background: what exists

- `WickTheorem<S>` (`SeQuant/core/wick.hpp`, `wick.impl.hpp`) enumerates
  pairwise contractions; each evaluates to an overlap `s` or Kronecker `δ`
  via `contract` (`wick.hpp:1626`). Non-pure spaces are projected with
  `δ(L,l)·s(l,r)·δ(r,R)` using temporary indices. `reduce` folds `s`/`δ`
  into index renames.
- `Vacuum::MultiProduct` exists (`SeQuant/core/attr.hpp:323`) but every
  classifier in `SeQuant/core/op.hpp:200-332` throws on it.
- `IndexSpaceRegistry` already separates vacuum occupancy from reference
  occupancy and defines hole/particle spaces; `make_mr_spaces`
  (`mbpt/convention.cpp:152`) gives core `O`, active `u`, hole `I = O ∪ u`,
  particle `A = u ∪ virt`, with `u = I ∩ A`.
- `ref_av` (`mbpt/op.cpp:1264-1480`) runs Wick relative to the core vacuum
  with all partial contractions and replaces each surviving normal operator
  by a full RDM γ projected onto `u`. No η, no cumulants.
- `mbpt/rdm.cpp` (`mbpt::decompositions`) converts cumulants (label `κ`) to
  densities (label `γ`) as a post-hoc rewrite.

## Theory used by the design

For GNO strings `{…}₁ {…}₂ … {…}ₙ`, the extended Wick theorem gives a sum
over all sets of disjoint contractions, where a contraction is either

- a **pair** joining a creator and an annihilator from different strings:
  `cre·ann → γ`, `ann·cre → η` (both over the active space; on core, γ = δ
  and on virtual, η = δ), or
- a **block** of k creators and k annihilators, k ≥ 2, drawn from at least
  two different strings, evaluating to the cumulant κₖ.

A contraction entirely within one string vanishes by definition of GNO. The
uncontracted operators remain as a GNO string; the sign is the parity of the
permutation that brings each contraction's legs together.

### No-double-counting rule

Each extended-Wick term corresponds to exactly **one** standard-Wick partial
contraction term: the one carrying exactly its pair contractions. Cumulant
blocks are built **only from operators that survive standard Wick**; a γ/η
pair emitted by Wick is never absorbed into a block, and a rank-1 "block" is
never formed in the pass. This makes the mapping one-to-one, so no 1/k!
weights are needed and `use_topology` remains valid.

Checks: ⟨{a†ₚa†_q a_r}{a_s}⟩ → the uncontracted term gives κ^{rs}_{pq}; the
p-s term leaves a rank-1 remainder and is dropped. ⟨{a†ₚa_q}{a†_r a_s}⟩ →
γη (both pairs) + κ₂ (no pairs); one-pair terms are dropped.

## A. Engine side

### Classifiers (`SeQuant/core/op.hpp`)

Fill the six `MultiProduct` cases. With R = `reference_occupied_space`
(core ∪ active; the γ-type contraction space) and U = `vacuum_unoccupied_space`
(active ∪ virtual; the η-type contraction space), so that active = R ∩ U and
the SR limit R = vacuum-occupied reproduces `SingleProduct` exactly:

| action | qp annihilator | qp creator |
|---|---|---|
| `Annihilate` | space ∩ U ≠ ∅ (pure: ⊆ U) | space ∩ R ≠ ∅ (pure: ⊆ R) |
| `Create` | space ∩ R ≠ ∅ (pure: ⊆ R) | space ∩ U ≠ ∅ (pure: ⊆ U) |

`qpannihilator_space`/`qpcreator_space` return the corresponding intersection.
`can_contract` is unchanged (left qp annihilator, right qp creator, spaces
meet), which yields exactly ann·cre over U and cre·ann over R. (The
registry's hole/particle spaces are *not* used: they describe where
excitation operators act, and in `make_mr_spaces` exclude the frozen
spaces that must still contract to δ.) An active operator is *both* a qp
creator and a qp annihilator; implementation must audit engine code that
assumes exclusivity (leftmost-op early return in `recursive_nontensor_wick`,
nop sorting).

### `contract` (`wick.hpp:1626`)

Structure unchanged. A `MultiProduct` branch picks the middle factor: `γ`
over the common space for cre·ann, `η` over the common space for ann·cre
(`s` otherwise). Projection Kroneckers are produced as today. γ/η are plain
tensors, so `reduce` only renames their indices; no changes to `reduce`.

A γ over R therefore stands for δ on core plus γ on active, and an η over U
for δ on virtual plus η on active; the pass splits them (B.1).

`SEQUANT_ASSERT` that when either index carries protoindices the common
space does not meet the active space.

### Connectivity

When the vacuum is `MultiProduct` the engine ignores `set_nop_connections` /
`set_nop_avoided_connections`; the pass enforces them (B.6), since a κ block
can connect nops that no pair connects.

### Deferred optimization

A "full contractions, but emit terms whose survivors are all active" engine
mode would prune earlier than running partial + filter. Not in milestone 1.

## B. The pass: `cumulant_expand`

Input: output of `WickTheorem::compute` with `full_contractions(false)`,
a Sum of Products `scalar × tensors × ≤1 NormalOperator` (survivors merged
and phased by `normalize`), plus provenance (B.2). Per term:

1. **Split mixed-space γ/η.** Pure core index → δ; pure active → keep;
   mixed (e.g. general `p`) → two-term sum via the `δ(p,u)·γ(u,…)` projection
   idiom, expanded and passed through `WickTheorem::reduce`. Likewise η with
   virtual ↔ core. Afterwards every γ/η is pure-active.
2. **Provenance.** Survivors keep their original `Index` (`normalize` only
   reorders). The wrapper canonicalizes the input itself, records
   `index → input nop ordinal`, and calls
   `compute(false, /*skip_input_canonicalization=*/true)`. Sum inputs are
   handled summand by summand. *To verify in planning:* nothing between input
   canonicalization and output renames a surviving index.
3. **Project survivors to active.** Pure core/virtual survivors cannot enter a
   block (full mode: term dies; partial mode: stays in the string). Mixed
   survivors split as in step 1.
4. **Enumerate block assignments.** Disjoint blocks among active survivors;
   each has k creators and k annihilators, 2 ≤ k ≤ `max_cumulant_rank`, and
   legs from ≥ 2 input nops. The rest stays in the string. Full mode keeps only
   assignments with an empty remainder.
5. **Emit.** One κₖ per block (creator legs and annihilator legs in the slot
   convention `ref_av` uses for γ, so that η = δ − γ holds slot-wise); the
   remainder rewrapped as `NormalOperator<S>` with vacuum `MultiProduct`;
   sign = parity of the permutation pulling each block's legs, in order, to
   the front of the survivor string.
6. **Filters.** Connectivity graph: a κ is a hyperedge among the nops its
   legs came from; resolved pair contractions appear as shared indices. Apply
   required/avoided connections; optionally drop disconnected terms.
7. Canonicalize and accumulate.

Cost is a sum over set partitions of the active survivors — factorial, as
accepted. `max_cumulant_rank` is the only pruning knob.

### Spin-free hooks

Both default to spin-orbital:

1. `contract`'s middle factor comes from a policy
   `(vacuum, left, right, common space) → ExprPtr`.
2. Block emission goes through `make_cumulant(rank, cre legs, ann legs) →
   ExprPtr`, plus a per-term weight hook (where a spin-free cycle factor
   would go).

## C. API and output tensors

```cpp
// SeQuant/core/wick_extended.hpp (symb target)
struct ExtendedWickOptions {
  bool full_contractions = true;
  std::optional<std::size_t> max_cumulant_rank;  // nullopt = unbounded
  bool eta_as_delta_minus_gamma = false;
  bool use_topology = true;
  // nop connections / avoided connections, as in WickTheorem
};
template <Statistics S>
ExprPtr extended_wick(ExprPtr input, const ExtendedWickOptions& = {});
template <Statistics S>
ExprPtr cumulant_expand(ExprPtr wick_output, /*provenance*/,
                        const ExtendedWickOptions&);  // exposed for tests
```

Requires `ctx.vacuum() == Vacuum::MultiProduct`, else throws
`sequant::Exception`. Bosons unsupported (as today).

**Tensors.** Labels shared with `mbpt/rdm.cpp` so its decompositions work on
the output unchanged: `γ` (1-RDM, Hermitian), `κ` (cumulants, rank ≥ 2,
Antisymm + Hermitian + column-symmetric — the same `TensorSymmetries` as
`rdm.cpp`). New: `η` (1-hole RDM, same symmetries as γ). `λ` is not used
because `mbpt` means CC deexcitation amplitudes by it. The label helpers move
out of `rdm.cpp`'s anonymous namespace into a shared header; `η` is
registered in the legacy op registry (`mbpt/context.cpp:148`).

**Rename.** `cumu_to_density` → `cumulant_to_density`; `cumu2_to_density`
and `cumu3_to_density` follow the same spelling (`cumulant2_to_density`,
`cumulant3_to_density`). Users: `rdm.{hpp,cpp}`, `tests/unit/test_mbpt.cpp`.
Own commit.

**mbpt.** `ref_av` dispatches on the context vacuum: `MultiProduct` →
`extended_wick` (+ `cumulant_to_density` behind an option for callers that
want densities only); `SingleProduct` → today's core-vacuum path, untouched.
`QuantumNumberChange`/`combine` (`mbpt/op.hpp:230`, `op.cpp:164`) gain a
`MultiProduct` branch only as far as `ref_av` screening needs; otherwise they
keep throwing.

## D. Tests and docs

Unit tests (`tests/unit/test_wick.cpp` new section "multiproduct vacuum";
`tests/unit/test_mbpt.cpp`):

1. Classifier truth table for `MultiProduct` (mirrors "Op contractions").
2. Identities via `EquivalentTo`: ⟨{a†a†a}{a}⟩ = κ₂; ⟨{a†a}{a†a}⟩ = γη + κ₂;
   a 3-body product with `max_cumulant_rank = 2` dropping exactly the κ₃
   terms; a partial-contraction product whose remainder is a GNO string.
3. Cross-check: for the MRSO cases (`test_mbpt.cpp:1298`),
   `extended_wick` → `cumulant_to_density` is equivalent to the current
   `ref_av` result.
4. Single-reference limit: with an empty active space, `MultiProduct` output
   equals `SingleProduct` output.
5. `use_topology` on/off agreement.

Docs: new section in `doc/developer/wick.rst` (classification table,
contraction values, the pass, no-double-counting rule); `MultiProduct`
paragraph in `doc/user/guide/context.rst`; compiled example under
`doc/examples/user/`. Fix the stray `Vacuum::SingleReference` in
`context.hpp:115` while that comment is touched.

## Commit structure (per AGENTS.md)

Separate commits, in order: rename `cumu*_to_density`; move RDM/cumulant
label helpers to a shared header and add `η`; `MultiProduct` classifiers;
`contract` middle-factor policy and connectivity skip; `cumulant_expand`
pass; `extended_wick` wrapper; `ref_av` dispatch; docs.
