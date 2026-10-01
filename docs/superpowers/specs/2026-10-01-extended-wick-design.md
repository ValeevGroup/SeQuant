# Extended Wick theorem (generalized normal order) — design

Date: 2026-10-01
Status: implemented (see "Deviations found during implementation")

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

## Background: the starting point

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
weights are needed. (The rule would also keep topological folding valid, but
`extended_wick` does not fold: see B.1.)

Checks: ⟨{a†ₚa†_q a_r}{a_s}⟩ → the uncontracted term gives κ^{rs}_{pq}; the
p-s term leaves a rank-1 remainder and is dropped. ⟨{a†ₚa_q}{a†_r a_s}⟩ →
γη (both pairs) + κ₂ (no pairs); one-pair terms are dropped.

## A. Engine side

### Classifiers (`SeQuant/core/op.hpp`)

The six `MultiProduct` cases use R = `reference_occupied_space`
(core ∪ active; the γ-type contraction space) and U = `vacuum_unoccupied_space`
(active ∪ virtual; the η-type contraction space), so that active = R ∩ U and
the SR limit R = vacuum-occupied reproduces `SingleProduct` exactly:

| action | qp annihilator | qp creator |
|---|---|---|
| `Annihilate` | space ∩ U ≠ ∅ (pure: ⊆ U) | space ∩ R ≠ ∅ (pure: ⊆ R) |
| `Create` | space ∩ R ≠ ∅ (pure: ⊆ R) | space ∩ U ≠ ∅ (pure: ⊆ U) |

`qpannihilator_space`/`qpcreator_space` return the corresponding intersection.
(The registry's hole/particle spaces are *not* used: they describe where
excitation operators act, and in `make_mr_spaces` exclude the frozen
spaces that must still contract to δ.) An active operator is *both* a qp
creator and a qp annihilator, so left-qp-annihilator/right-qp-creator alone
would also admit cre·cre and ann·ann pairs of active operators.
`can_contract` therefore additionally requires, under `MultiProduct`, that
the two operators have opposite actions; with the space condition (left qp
annihilator, right qp creator, spaces meet) this yields exactly ann·cre over
U and cre·ann over R. The other engine sites that branch on qp character
(leftmost-op early return in `recursive_nontensor_wick`, nop sorting) are
unaffected by an operator being both.

### `contract` (`wick.hpp:1626`)

Structure unchanged. A `MultiProduct` branch of the
`WickTheorem::contraction_value` hook picks the middle factor: `γ`
over the common space for cre·ann, `η` over the common space for ann·cre
(`s` otherwise). Projection Kroneckers are produced as today. γ/η are plain
tensors, so `reduce` only renames their indices; no changes to `reduce`.

A γ over R therefore stands for δ on core plus γ on active, and an η over U
for δ on virtual plus η on active; the pass splits them (B.2).

`SEQUANT_ASSERT` that when either index carries protoindices the common
space does not meet the active space.

### Connectivity

When the vacuum is `MultiProduct` the engine ignores `set_nop_connections` /
`set_nop_avoided_connections` (`WickTheorem::pairwise_connectivity()` is
false); the pass enforces them (B.6), since a κ block can connect nops that
no pair connects.

### Deferred optimization

A "full contractions, but emit terms whose survivors are all active" engine
mode would prune earlier than running partial + filter. Not in milestone 1.

## B. The pass: `extended_wick` and `cumulant_expand`

`extended_wick` runs the standard theorem with `full_contractions(false)` and
post-processes its output, a Sum of Products `scalar × tensors × ≤1
NormalOperator` (survivors merged and phased by `normalize`). Sum inputs are
expanded and handled summand by summand. Per input term:

1. **Provenance.** The wrapper canonicalizes the term and records
   `index → input nop ordinal` for every operator index; an index shared by
   two nops (a summed index) is first renamed in the later one, and the δ
   binding the two names multiplies the term after the filters (6), so it
   connects nothing. It then runs
   `WickTheorem` on the **bare operator sequence** (a
   `NormalOperatorSequence` of the term's nops) and multiplies the c-number
   factors back in afterwards. In that run every operator index is external,
   so no partially contracted term is canonicalized and no surviving index
   is renamed; running on the full Product would canonicalize the terms that
   carry c-number factors and rename the operator indices that are dummies
   of those factors, and `set_external_indices` does not prevent that.
   Consequence: **`use_topology` has no effect.** Topological pruning needs
   equivalent operators, i.e. dummy indices; provenance needs named ones.
2. **Split mixed-space γ/η.** Pure core index → δ; pure active → keep;
   mixed (e.g. general `p`) → sum via the `δ(p,u)·γ(u,…)` projection idiom.
   Likewise η with virtual ↔ core. Afterwards every γ/η is pure-active.
3. **Project survivors.** A mixed survivor is projected onto its active part
   and, with partial contractions, also onto its core and virtual parts, each
   projection δ-bound to the original index. A part that is not a registered
   space (e.g. {a, g} in `make_mr_spaces`) is split over its base spaces. The
   pieces of steps 2-3 pass through `WickTheorem::reduce` with every input
   index (operator and c-number) kept fixed, so each projected index stays
   δ-bound to an input op index and inherits its provenance.
4. **Enumerate block assignments** (`cumulant_expand`). Disjoint blocks among
   the *active* survivors; each has k creators and k annihilators,
   2 ≤ k ≤ `max_cumulant_rank`, and legs from ≥ 2 input nops. Core and
   virtual survivors are never legs (full mode: term dies; partial mode: they
   stay in the string). The rest stays in the string. Full mode keeps only
   assignments with an empty remainder.
5. **Emit.** One κₖ per block (bra = annihilator legs, ket = creator legs,
   the slot convention of γ, so that η = δ - γ holds slot-wise); the
   remainder rewrapped as `NormalOperator<S>` with vacuum `MultiProduct`;
   sign = parity of the permutation pulling each block's legs, in order, to
   the front of the survivor string.
6. **Filters.** Connectivity is hyperedge connectivity: a factor the theorem
   produces (γ, η, κ, δ or overlap) connects all the input nops its indices
   came from, so a κ connects every nop it has a leg from and a resolved pair
   contraction connects its two; any other tensor, e.g. a coefficient that
   spans several nops, connects none. Terms that miss a required connection
   or realize an avoided one are dropped.
7. Simplify and accumulate.

After all terms: every η is optionally rewritten as δ - γ
(`eta_as_delta_minus_gamma`), then `apply_dummy_deltas` reduces each term by
its own index counts, applying every δ over a summed index (e.g.
⟨h^p_q {a†_p}{a_q}⟩ gives h^O_O, not h^p_q δ^O_p δ^q_O), and the result is
simplified. A δ between two external indices is kept.

Cost is a sum over set partitions of the active survivors — factorial, as
accepted. `max_cumulant_rank` is the only pruning knob.

### Spin-free hooks

Both default to spin-orbital:

1. `contract`'s middle factor comes from
   `WickTheorem::contraction_value(bra, ket, left_is_annihilator, vacuum)`.
2. Block emission goes through `detail::block_value(block nop) → ExprPtr`
   (`density::make_cumulant` by default), plus a per-term weight hook
   `detail::term_weight` (where a spin-free cycle factor would go).

Spin-free itself is unsupported: `WickTheorem::compute` throws for a
`MultiProduct` vacuum with `SPBasis::Spinfree`.

## C. API and output tensors

```cpp
// SeQuant/core/wick_extended.hpp (symb target)
struct ExtendedWickOptions {
  bool full_contractions = true;
  std::optional<std::size_t> max_cumulant_rank;  // nullopt = unbounded
  bool eta_as_delta_minus_gamma = false;
  bool use_topology = true;  // documented as having no effect (B.1)
  container::svector<std::pair<std::size_t, std::size_t>> nop_connections;
  container::svector<std::pair<std::size_t, std::size_t>>
      nop_avoided_connections;
};
using OpProvenance = container::map<Index, std::size_t>;
template <Statistics S>
ExprPtr extended_wick(ExprPtr input, const ExtendedWickOptions& = {});
template <Statistics S>
ExprPtr cumulant_expand(const ExprPtr& wick_output,
                        const OpProvenance& provenance,
                        const ExtendedWickOptions&);  // exposed for tests
```

Requires `ctx.vacuum() == Vacuum::MultiProduct`, else throws
`sequant::Exception`. Bosons unsupported (as today). `use_topology` is kept
in the struct so that callers (e.g. `ref_av`) can forward their setting, but
it is not used.

**Tensors.** Labels shared with `mbpt/rdm.cpp` so its decompositions work on
the output unchanged: `γ{ann;cre}` (1-RDM; Hermitian, column-symmetric),
`κ{ann…;cre…}` (cumulants, rank ≥ 2; Antisymm + Hermitian +
column-symmetric). New: `η{ann;cre}` (1-hole RDM, same symmetries as γ). In
the output every γ, η and κ index is active. `λ` is not used because `mbpt`
means CC deexcitation amplitudes by it. The labels, symmetries and factories
(`density::make_rdm`, `make_hole_rdm`, `make_cumulant`, `rdm_from_nop`) live
in `SeQuant/core/density.hpp`, used by `contract`, the pass and `ref_av`
alike; `η` is registered in the legacy op registry (`mbpt/context.cpp`).

**Deserialization.** A tilde operator (`ã`/`b̃`) is deserialized with the
default context's vacuum when that is not `Physical`, and with
`SingleProduct` otherwise (`tilde_vacuum` in
`io/serialization/v1/ast_conversions.hpp`), so that under `MultiProduct`
parsed GNO operators carry the right vacuum.

**Rename.** `cumu_to_density` → `cumulant_to_density`; `cumu2_to_density`
and `cumu3_to_density` follow the same spelling (`cumulant2_to_density`,
`cumulant3_to_density`). Users: `rdm.{hpp,cpp}`, `tests/unit/test_mbpt.cpp`.
Own commit.

**mbpt.** `ref_av` (and `vac_av`, at the tensor and operator level) dispatch
on the context vacuum: `MultiProduct` → `extended_wick` with
`full_contractions = true`, forwarding the required and avoided
connections; `SingleProduct` → the core-vacuum path, untouched. The result
keeps its cumulants. `mbpt::decompositions::cumulants_to_densities(expr)`
rewrites every κ (recognized by label) of rank ≤ 3 into densities via
`cumulant_to_density`/`cumulant2_to_density`/`cumulant3_to_density`, then
expands and simplifies; a κ of rank > 3 throws `sequant::Exception`.
Multi-body γ from these helpers are antisymmetric (`density::cumulant_symmetries`),
matching the multi-body γ that the core-vacuum `ref_av` builds.

Under `MultiProduct` the mbpt operators (`ã`: h, f, g, t, …) are
generalized-normal-ordered relative to the reference, so `ref_av` evaluates
⟨{H}_ref {T}_ref⟩_ref, a different quantity from the `SingleProduct` path's
⟨{H}_core {T}_core⟩_ref. The two agree for products of elementary operators
(bare a†/a strings), whose normal ordering is immaterial.

`QuantumNumberChange::size()` and `combine()` (`mbpt/op.hpp`, `op.cpp`) gain
`MultiProduct` branches for `ref_av` screening: core and virtual base spaces
as for `SingleProduct`; an active base space bounds the number of removed
creator/annihilator pairs by min(total creators, total annihilators) of the
two sides, since a cumulant may take any equal number of creators and
annihilators from either side.

## D. Tests and docs

Unit tests (`tests/unit/test_wick.cpp` section "multiproduct vacuum";
`tests/unit/test_wick_extended.cpp`; `tests/unit/test_mbpt.cpp` section
"MRSO-MultiProduct"):

1. Classifier truth table for `MultiProduct` (mirrors "Op contractions"),
   including the opposite-action requirement of `can_contract`.
2. Identities via `EquivalentTo`: ⟨{a†a†a}{a}⟩ = κ₂; ⟨{a†a}{a†a}⟩ = γη + κ₂;
   a 3-body product with `max_cumulant_rank = 2` dropping exactly the κ₃
   terms; a partial-contraction product whose remainder is a GNO string;
   mixed-space inputs (general indices split into core δ + active γ, and
   projected survivors under partial contractions); δs over summed indices
   applied; connectivity, including connection through a κ; η = δ - γ.
3. Cross-check against the core-vacuum path for products of **elementary**
   operators (c-number coefficient × bare a†/a string), whose normal order is
   immaterial: `cumulants_to_densities(ref_av)` under `MultiProduct` is
   equivalent to `ref_av` under `SingleProduct` (after both are rewritten in
   base spaces), for one-body, two-body (up to κ₂) and two-body × one-body
   (up to κ₃) strings, without connectivity (which means different things on
   the two paths). Products of mbpt operators are not compared, since they
   are different operators on the two paths (C, **mbpt**).
4. `ref_av(h(2)·t(2))` runs under `MultiProduct` and returns a nonzero
   result, which reaches κ₄, so `cumulants_to_densities` throws on it.
5. Single-reference limit: with an empty active space, `MultiProduct` output
   equals `SingleProduct` output. (Not yet covered by a test.)

There is no `use_topology` on/off test for `extended_wick`, since the option
has no effect (B.1); the mbpt "topology on/off agree" section only guards
`ref_av`'s dispatch.

Docs: section "The extended Wick theorem" in `doc/developer/wick.rst`
(engine changes, the pass, no-double-counting rule); `MultiProduct` in the
vacuum list of `doc/user/guide/context.rst`; `ref_av`'s dispatch in
`doc/user/guide/operator.rst`; a getting-started subsection in
`doc/user/getting_started/wick.rst` with the compiled example
`doc/examples/user/getting_started/extended_wick.cpp`. Fix the stray
`Vacuum::SingleReference` in `context.hpp` while that comment is touched.

## Commit structure (per AGENTS.md)

Separate commits, in order: rename `cumu*_to_density`; move RDM/cumulant
label helpers to a shared header and add `η`; `MultiProduct` classifiers;
`contract` middle-factor policy and connectivity skip; `cumulant_expand`
pass; `extended_wick` wrapper; `ref_av` dispatch; docs.

## Deviations found during implementation

Relative to the design as first approved:

- `can_contract` is not unchanged: under `MultiProduct` it also requires
  opposite actions, since an active operator is both a qp creator and a qp
  annihilator (A).
- `extended_wick` runs `WickTheorem` on the bare operator sequence, not on
  the canonicalized Product with `skip_input_canonicalization`, because the
  latter renames surviving indices and breaks provenance; as a result
  `use_topology` has no effect (B.1).
- Cumulant legs are active operators only; mixed-space γ/η and survivors are
  split into pure pieces, over base spaces where a part is not registered
  (B.2-B.4).
- Connectivity is enforced in the pass as hyperedge connectivity (B.6).
- The δs over summed indices left by the projections are applied at the end
  (`apply_dummy_deltas`), after the optional η → δ - γ rewrite (B).
- Densities are not produced behind a `ref_av` option; the separate
  `cumulants_to_densities` (rank ≤ 3, throws above) does it (C).
- A deserialized tilde operator takes the context's vacuum (C).
- Under `MultiProduct` `ref_av` evaluates generalized-normal-ordered mbpt
  operators, so the cross-check with the core-vacuum path holds for
  elementary operators only (C, D.3); the `use_topology` test is dropped and
  the h(2)·t(2) test added (D).
- Found along the way and fixed: `antisymmetrize` preserves tensor
  symmetries and recognizes repeated products by index labels;
  `cumulant3_to_density` returns γ₃ minus its factorizable part; multi-body γ
  from the cumulant decompositions are antisymmetric.
