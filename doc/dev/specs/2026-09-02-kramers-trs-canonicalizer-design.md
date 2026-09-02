# Kramers time-reversal (TRS) canonicalizer — design

Branch: Kramers round 2 (builds on the conjugation-eval PR-2 work).
Date: 2026-09-02.

## Problem

On the Kramers-restricted CC/CSV path every tensor with Kramers-flavored
slot indices (`a↑`, `i↓`, CSV virtuals `a⇑<..>`) satisfies a time-reversal
(TRS) data identity: the spelling with EVERY flavored slot flipped equals
`phase · conj(T)`, with `phase = (−1)^{#slots flipped from down}`. This is the
`kramers_transform` convention (spinor.hpp:164-186) and it reproduces the
PNS coefficient signs (`C↓↓ = +conj C↑↑`, `C↓↑ = −conj C↑↓`) and the DF /
Fock identities (`g↓↓ = +conj g↑↑`, `g↓↑ = −conj g↑↓`) alike.

The canonicalizer does not know this identity. Consequences measured on
dch (PNS-MP1, 2026-09-02):

- the flavor-resolved residual carries all four C families as distinct
  leaves (R2: 400 of 800 C leaves down-row), so the conj-aware eval CSE
  never pairs cross-flavor subnets;
- the energy TRS fold pairs 24/36 terms; the 12 survivors are partners
  whose spellings differ by a component orientation the rebase could not
  choose;
- the within-block twin fold on the self-conjugate residual block finds
  0/48 literal pairs: the block's closure `R = S(R)` is realized through
  the C/f/g leaf identities, not a symbolic bijection.

`kramers_internal_rebase` is a partial, external post-pass: it flips whole
internal components by a lexicographic / canonical-config criterion and
freezes components touching externals. It cannot be idempotent with the
canonicalizer and cannot serve the eval leaf boundary.

User specification (2026-09-02): non-Kramers-canonical configs of ANY
tensor (C, g, f, t) must not appear in the equations — only their
conj+phase variants; the DF auxiliary index never carries a Kramers label.

## Design

### 1. A declared tensor symmetry

`enum class KramersSymmetry { Nonsymm, TimeReversal }` in `core/attr.hpp`
next to `ColumnSymmetry`; orthogonal to `BraKetSymmetry` (they compose:
braket-Conjugate fold = swap+conj, Kramers fold = flavor flip+phase·conj).

Stored per `Tensor` (`kramers_symmetry_`), threaded through both main
constructors, the `Hermiticity` constructors, `with_slots`, `_hash_value`,
`static_equal` / `static_less_than`; exposed on `AbstractTensor`
(`_kramers_symmetry()`) with a free `kramers_symmetry(const AbstractTensor&)`
wrapper (the canonicalizer only uses the free wrappers). Serialization gains
an optional fourth attribute letter in the `:A-C-S` spec (`T` = TimeReversal,
`N` = Nonsymm, absent = Nonsymm) so tests can be deserialize-driven.

Assignment: the mbpt `OpRegistry` defaults every registered op to
`TimeReversal` when the operator's slots are Kramers-flavored (assigned at
`OpMaker::operator()` next to `op_hermiticity`); `csv_transform` and
`density_fit` propagate the attribute from the source tensor to the C and
DF-factor leaves they mint. Reserved / operator tensors (`Â`, `Ŝ`,
transposition) are never foldable (`braket_orientation_pinned` guard).

### 2. Flavor flip as a core operation

`SeQuant/core` cannot depend on `domain/mbpt`, so the flip is a registry
query, populated by the convention that creates the spin clones
(`convention.cpp:56`): `IndexSpaceRegistry::kramers_partner(const
IndexSpace&) -> std::optional<IndexSpace>`. Spaces without a registered
partner (the DF aux `Κ`, any spin-bit-free space) are unflavored and are
never flipped. `Index` gets `kramers_flipped(const IndexSpaceRegistry&)`
returning the partner-space index with the same ordinal, proto indices
mapped recursively.

### 3. Single-tensor fold (eval leaf boundary)

`canonicalize_kramers(AbstractTensor&, bool fold) -> int phase` in
`tensor_canonicalizer.cpp`, modelled on `canonicalize_braket`: (i) unfold —
nothing to do, the marker is orthogonal; (ii) decide: if the tensor's FIRST
flavored slot (bra, ket, aux order after the braket fold) is down, flip every
flavored slot index to its partner, toggle the conj marker, and return
`(−1)^{#down slots flipped}`; otherwise return +1. `TensorBlockCanonicalizer::apply`
calls it after `canonicalize_braket`. The eval leaf constructor
(`EvalExpr::EvalExpr(Tensor const&)`) folds the result into the existing
`CanonTransform` as `{.phase = phase, .conj = true}` — no new field. Result:
`C↓↑` and `C↑↓` occupy ONE cache slot (transform `{conj, −1}`), likewise
`C↓↓`/`C↑↑` (`{conj, +1}`), and the MPQC server is only ever asked for the
up-row arrays.

**Eval-leaf fold is an explicit opt-in** (`CanonicalizeOptions::
fold_kramers_eval_leaves`, default No). Measured on dch (2026-09-02): with
the leaf fold on, PNS-MP1 drifts to −1.0115 and diverges. Cause: the
evaluator serves a leaf from `as_tensor()` and applies the leaf's
transform to the served value under the leaf's own annotation; a flavor
flip changes the LABELS (a↓_1 → a↑_1), which no transform can relabel, so
the down block is served and then conj/phase-transformed again. The leaf
fold is only correct once the provider aliases the Kramers partner block
(serving-level aliasing, T19), at which point the option can be enabled.

### 4. Network fold (symbolic equations)

Under `CanonicalizeOptions::fold_kramers` (mirrored as
`CanonicalizeSlotsOptions::fold_kramers`), `TensorNetworkV3::kramers_orient`
runs as a PRE-PASS before the ordinary (flavor-aware) canonicalization:

- the TRS identity is a whole-tensor identity, so flipping a tensor flips
  every flavored slot it has; every tensor sharing a flavored dummy (as a
  slot or as a proto index) must flip with it. The orientation unit is
  therefore a connected component of tensors joined by shared flavored
  indices (union-find);
- a component is pinned if it contains a named (external) index — their
  flavor is fixed — or a tensor without a Kramers identity
  (`kramers_foldable` false);
- a free component is flipped iff the flipped spelling has fewer
  down-first tensors, ties broken by the sorted `(label, per-slot flavor,
  marker)` fingerprint, which is index-label-independent so both spellings
  of a component land on the same orientation;
- flipped tensors get partner-space indices (proto indices recursively)
  and the conjugation marker; the phase `(-1)^(#down slots)` per flipped
  tensor multiplies into the byproduct (a closed component's phase is +1:
  each dummy is down in two slots);
- `canonicalize_slots` reports the flipped input ordinals as
  `kramers_flipped_tensors` next to `conjugated_tensors` and folds the
  phase into `phase`.

No flavor-blind graph colouring is needed (an earlier draft of this
section proposed it): once the pre-pass has oriented every free component,
the flavor-aware graph identifies the two spellings by construction, and
the conj colouring of the marker (#57) keeps `C·C*` and `C*·C` distinct.
Consequently the fold does not run inside `TensorBlockCanonicalizer` by
default (`fold_kramers(false)`): applied per tensor inside a network it
would flip one tensor's dummies but not its partner's. The eval leaf opts
in and, after the fold has fixed the hash and the `{conj, phase}`
transform, restores the as-written flavors (parents contract by label).

**Braket orientation is Kramers-aware.** Measured on dch: the network
fold flipped an all-down energy component to all-up, but the Conjugate
bra/ket fold that follows re-oriented C and g so a down index sat in the
bra again — "first flavored slot" is not invariant under the swap a
Hermitian tensor may make. `canonicalize_braket` therefore prefers, for
`TimeReversal` tensors, the orientation whose bra carries fewer
down-flavored indices (label/permutation invariant, value-exact via the
marker; ties fall through), and `kramers_down_first()` judges a spelling
on its braket-canonical orientation. Consequence: a Hermitian C↓↑ reaches
the up row by the braket move (conj, no phase) rather than the Kramers
flip; only Nonsymm-braket tensors (t) need the flip itself.

### 5. Scope and expectations

- Energy: the 12 unpaired terms are exactly component-orientation twins;
  they pair once the fold is idempotent with canonicalization.
- Residual blocks: all internal components become canonical; components
  anchored to externals keep the externals' flavors (correct — their
  images are the reconstructed partner blocks). Whether the self-conjugate
  block's within-block twins pair after the combined braket+Kramers folds
  is an open question to be MEASURED, and independently verified by the
  numerical term study (`t' == sign·conj(perm(t))` on evaluated terms).
- `kramers_internal_rebase` becomes redundant; it stays until measured
  equal, then is retired.

## Non-goals

No change to the tracer's kept-block set, to the eval closure gadgets
(KrRecon / KrStab / the twin-fold closure), or to the DF/CSV emission
structure. No Θ-style notation anywhere: the operator is 𝒯 / "time
reversal".
