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

### 4. Network fold (symbolic equations)

Under a new `CanonicalizeOptions::fold_kramers` (mirrored in
`CanonicalizeSlotsOptions` and `CreateGraphOptions`, forwarded into
`canonicalize_graph` which today receives only `method` and
`ignore_named_index_labels`):

- index-vertex colours become flavor-blind for flavored spaces (mask the
  spin qns bits in `VertexPainterImpl::operator()(const Index&)`), and the
  tensor core colour is normalized over flavor the way bra/ket ranks are;
  named (external) indices keep their identity — their flavor is fixed;
- two spellings that differ by flipping the flavors of a connected component
  of dummies now produce the same canonical graph;
- after the canonical labeling, each connected component of flavored dummies
  gets ONE orientation: components touching a named index are fixed by it;
  free components take the orientation with fewer leaves whose first
  flavored slot is down, tie-broken by the `(label, flavor bits)`
  fingerprint — bit-identical to `kramers_rebase_term` so the two agree;
- the byproduct is reported per tensor like `conjugated_tensors`:
  `kramers_flipped_tensors` (ordinals) plus the accumulated `phase`;
  `canonicalize` respells flipped tensors (partner-space indices, conj
  marker toggled) and multiplies the phase into the scalar.

The conj colouring of the marker (#57) is kept: `C·C*` and `C*·C` stay
distinct; a Kramers-flipped-and-marked tensor is coloured as marked.

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
