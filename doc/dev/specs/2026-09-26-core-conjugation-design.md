# Core conjugation: redesign of the symbolic conjugation layer

(Working title was "adjoint-only"; the model has two core states, adjointed and conjugated, see section 2.)

> **Status.** This document is the design record: the type argument, the
> prior art, the experiment, and the decisions. The maintained description of
> the model as it exists in the code is `doc/dev/conjugation-model.md`. Where
> the implementation settled something differently from the text below, an
> "As built" note follows the affected paragraph; the text itself is kept as
> written. The two documents this one supersedes or extends
> (`2026-09-21-tensor-value-modifiers.md`,
> `2026-09-22-conjugation-parity-and-signed-braket-symmetry.md`) are review
> documents and are not kept in the tree.

Design spec for reworking the symbolic half of PR #602 (`kshitij/feature/
conjugation-symbolic`, head `3ea0a6e93`) on top of #632 (`fix602/review-batch-1`,
head `7c51815be`). Distilled from the review of #602 on 2026-09-24/25 and the
experiment of 2026-09-26 (`review/experiment-adjoint-only-conj.patch` in the
review checkout). Supersedes the state model of
`2026-09-21-tensor-value-modifiers.md`; keeps the traits of
`2026-09-22-conjugation-parity-and-signed-braket-symmetry.md`.

## 1. Problem

#602 gives `Tensor` two array-level state bits, `conjugated_` (elementwise
conjugation) and `transposed_` (slots read in exchanged roles), the Klein
four-group `ValueModifier {None, Conjugate, Transpose, Adjoint}`, a
normalization table that reduces the bits against the tensor's symmetries, and
a canonicalization fold that rewrites one orientation of a Hermitian tensor as
the other orientation's elementwise conjugate with the *same* slots
(`t{a;i}` becomes `t^*{i;a}`). Three things are wrong with that model.

**The transpose is not a tensor over a complex basis.** A bra slot transforms
with `U†` under a basis change and a ket slot with `U`. `t^T{i;a}` puts the
old ket index in the bra position, where it transforms with `U^T`; the object
is not a (1,1) tensor over that basis and is the matrix of no operator. Only
over a real basis, where `U† = U^T`, is the transpose a respelling of the same
tensor, and there it coincides with the adjoint. That is why operators have an
adjoint and no transpose. In the IR a transpose is a layout permutation.

**The elementwise conjugate with the same slots is the same type error.**
`conj <i|t|a> = <a|t†|i>`: the conjugate of a matrix element is the adjoint's
matrix element with the slots exchanged. Writing it with the original slots,
`t^*{i;a}`, is non-covariant. This shows up in networks: the fold turns the
covariant `t{a;i} u{i;a}` into `t^*{i;a} u{i;a}`, where `i` sits in both bras.
The spelling is only "valid" because the marker secretly means the other
orientation, and it is why `TensorNetworkV3` grew `fold_conjugate_braket`, the
marker colouring, the unfold loop, the post-relabel refold, the second pass and
the relaxed dummy-edge assertion, why `binarize` counts slots on an "unfolded"
copy, and why the eval boundary disables the fold and re-lowers marked leaves.

**Two conjugations were conflated.** The operator-level conjugation
`K O K⁻¹` is basis-independent; the tensor's answer to it is the parity trait
(`±O`, or an unrelated operator). The array-level one, elementwise conjugation
of `<p|O|q>`, is `(O†)` in transposed spelling. They coincide only over a real
basis, where `Kp = p`. Over a complex basis the array conjugate is the matrix
of `K O K⁻¹` between the conjugated basis functions, not an operation on the
array in the same basis.

**What is symbolic, and what carries it.** `t{i;a} + t{a;i} = 2 Re t{i;a}`
for Hermitian `t` is a symbolic fact, and the `Re`/`Im` fold and
Hermitian-network recognition need it. It is carried by `Hermiticity`,
`ConjugationParity` and the adjoint state: the conjugate of a term is the
adjoint of every factor with the head's bra and ket exchanged, which stays a
covariant network, and the fold compares canonical forms of a term and of its
conjugate. The experiment confirmed it: with the marker fold switched off and
`conjugate(ExprPtr)` routed through the adjoint, 271 of 277 test cases pass and
the six failures assert the marker mechanism itself (`fold_conjugate_pairs`,
`hermitian_network_recognition`, the signed BTAS networks and the ITF fixtures
all pass).

**Where a conjugation state is still needed.** Anything with no bra/ket pair
to absorb conjugation into: `Constant`, `Variable`, `Power`, and a `Tensor`
with no bra and no ket slot (slotless or aux-only, `w{;;x}`). `conj λ` and
`conj w{;;x}` must be writable in arbitrary and canonical form. That state is
the adjoint in the case where the adjoint has nothing to exchange.

**Prior art.** xAct's xTensor treats complex conjugation as a property of the
tensor object (`Dagger`): a tensor is real (self-conjugate) or complex, in
which case its conjugate is a distinct object, and index types map to their
conjugate types (undotted to dotted spinor indices). SymPy's `physics.quantum`
makes `Dagger` the primitive on operators and keeps `conjugate`/`transpose`
for matrices, i.e. arrays. Relativistic quantum chemistry writes the
time-reversal partner of an orbital as a barred index. In all three,
conjugation acts on the object and on index types, never on an array in
place; the elementwise conjugate is derived and covariant only where the
conjugate index type is the original one. This spec follows that model.

## 2. The model

A `Tensor` is `<bra|O|ket>`: a core `O` with the traits `Hermiticity` and
`ConjugationParity`, slots, and two commuting Z2 states on the core, *adjointed*
(`†`, mark `⁺`) and *conjugated* (`*`, the complex conjugation `K O K⁻¹` in the
coordinate representation, the operation `ConjugationParity` is defined
against). Both are covariant in every basis because they act on the operator,
not on the array. Together they are the Klein four-group `{1, †, *, †*}`;
`†*` is `K O† K⁻¹`. `K` is not time reversal: for spin-½ particles `Θ = U_T K` with
`U_T = −iσ_y`, so `Θ` differs from `K` by a unitary spin rotation, `Θ² = −1`,
and the `Θ`-parity of an operator differs from its `K`-parity (`σ_x`, `σ_z`
are `K`-even and `Θ`-odd). Time reversal is `K` composed with that unitary,
and belongs to the Kramers work (#566).

### 2.1 The adjoint state

`t⁺{q;p}` denotes `conj t{p;q}`, the matrix element of `O†`. It is covariant
for every slot structure: rectangular tensors, half tensors
(`X⁺{;a;x} = conj X{a;;x}`), tensors with aux slots
(`B⁺{a;i;x} = conj B{i;a;x}`; aux slots are array-like and stay in place), and
bra/ket-less tensors, where the exchange is empty and `w⁺{;;x} = conj w{;;x}`.

The state is normalized against the hermiticity, as `adopt_adjoint_mark` and
`normalize_value_modifier` do today for the `⁺` case, and nothing else:

| `Hermiticity` | `t⁺` | state kept? | sign returned by `adjoint()` |
|---|---|---|---|
| `Hermitian` | `t` | no | +1 |
| `AntiHermitian` | `−t` | no | −1 |
| `NonHermitian` | a distinct array | yes | +1 |

For a bra/ket-less tensor the same table reads "real", "imaginary",
"unknown": a real aux-only array is declared `Hermitian`. The parity and the
basis field do not enter the state at all; they enter only the derived
observables (`BraKetSymmetry`, `ConjugationSymmetry`), exactly as in the
2026-09-22 spec, and those observables tell canonicalization which free or
signed slot exchanges are symmetries.

### 2.2 The conjugated state

`t꙳{p;q}` denotes `<p|K O K⁻¹|q>`, the matrix of the conjugated operator in
the same basis with the same slots. It is normalized against the parity:

| `ConjugationParity` | `t꙳` | state kept? | sign returned by `kconjugate()` |
|---|---|---|---|
| `Even` | `t` | no | +1 |
| `Odd` | `−t` | no | −1 |
| `None` | a distinct array | yes | +1 |

> **As built.** `Tensor::kconjugation_sign()` reads the parity directly, in
> every basis (`K O K⁻¹ = ±O` is a statement about the operator, not about
> the matrix), and not the basis-derived `ConjugationSymmetry`: `t꙳` on a
> default (`Even`) tensor is `t` over a complex basis as well, and a kept star
> requires parity `None`. For a tensor with no bra and no ket slot, whose
> adjoint is its conjugate, the hermiticity ("real", "imaginary", "unknown")
> normalizes `꙳` first and the parity applies where the hermiticity is
> indefinite.

Why the state is needed and not only the trait: in the real-basis,
complex-amplitude regime (magnetic response, finite-field or complex CC) the
perturbation is `K`-odd and the first-order response imaginary, both covered
by the parity, but higher-order and mixed-perturbation amplitudes are complex
with no definite parity; `t^*{a;i}` is then a distinct atom, covariant over a
real basis (`Kp = p`) with the slots in place, which is the spelling the
physics uses and which the adjoint spelling `t⁺{i;a}` would replace with a
slot exchange on every amplitude.

**Value versus operator.** The conjugate of the *value* `t{p;q}` is `t⁺{q;p}`
in every basis, so `conjugate(ExprPtr)` is `adjoint(ExprPtr)` on c-number
content and needs nothing about the basis. `t꙳{p;q}` is the *operator*
conjugated, `kconjugate`: over a real basis it is the conjugated value with
the slots in place (and the coset rule of 2.3 then prefers it), over a complex
basis it is a different operator's matrix and not the conjugated value. A star
written on a complex-basis tensor to mean "conjugate the array" denotes the
wrong object, and the deserializer cannot tell; this is inherent in the model
and must be documented on the mark.

`K` is antiunitary, so on a c-number coefficient it is complex conjugation.
On operator-valued content it acts through the coefficients:
`K (Σ t_{pq} a†_p a_q) K⁻¹ = Σ conj(t_{pq}) a†_{Kp} a_{Kq}`, and when the
index space `S` is closed under `K` the sum over `{Kp}` is a sum over `S`, so
re-indexing gives `Σ <p|K t K⁻¹|q> a†_p a_q = Σ t꙳{p;q} a†_p a_q`: the string
is unchanged and the star sits on the tensor. `NormalOperator::kconjugate()`
is therefore the identity (returns +1). The precondition is `K`-closure of
every index space: real-field spaces are closed; a Kramers-paired space is
(`K S = U_T⁻¹ Θ S = S`); a generic complex subspace is not, and there
`K E K⁻¹` has no expression in the same labels (the "no conjugate space"
case of 2.3). `†` and `K` commute.

### 2.3 Conjugation and transposition are derived

The elementwise conjugate of a matrix element has two spellings:

- always, `conj t{p;q} = t⁺{q;p}`;
- through the antiunitary's action on the basis: `conj t{p;q} =
  <Kp|K O K⁻¹|Kq>`. Over a real basis `Kp = p`, so `conj t{p;q} = t꙳{p;q}`.
  Over a Kramers-paired basis `Kp` is not a basis function but `Θp = U_T Kp`
  is, the barred partner `p̄`, and `conj t{p;q} = <p̄|Θ O Θ⁻¹|q̄>`: the
  relation goes through the `Θ`-conjugated core, which is the `*` state
  composed with the unitary `U_T`, and through the `Θ`-parity (#566's
  index-level relation, `conj t{p;q} = ±t{p̄;q̄}`, is that statement for a
  `Θ`-definite operator). Over a generic complex basis `Kp` lies in the
  conjugate space; until such a space exists no elementwise identity is
  available, and `ConjugationSymmetry::NonSymm` is correct.

Hence `ConjugationSymmetry` is the `K`-parity read through the basis (an
elementwise sign exactly where `Kp = p`), and stays a derived observable. Over a real basis with both traits indefinite,
`t⁺{q;p}` and `t꙳{p;q}` are two spellings of one value; the coset rule is to
eliminate `†` in favour of `*` there (`adjoint()` over a real basis exchanges
the slots and toggles `*`), so that slots stay in place and real-basis
expressions never show `⁺` where a star is meant.

> **As built.** `normalize_states()` runs the identification in both
> directions: where the hermiticity is definite it consumes the star at the
> hermiticity's sign and exchanges the bundles instead (ahead of the parity),
> so a real-basis `kconjugate` and `conjugate` of a Hermitian or
> anti-Hermitian tensor produce one spelling.

- `conjugate(ExprPtr)` is the complex conjugate of the value: on c-number
  content it is `adjoint(ExprPtr)` (one implementation; the head exchanged for
  an open expression), and it throws on operator-valued content, which has no
  value to conjugate. `kconjugate(ExprPtr)` is conjugation of the operator,
  `K E K⁻¹`: it toggles `*` on every tensor factor and conjugates every
  coefficient, order preserving. Over a real basis the two coincide (the
  coset rule turns the adjoint spelling into stars); on scalars all of
  `adjoint`, `conjugate`, `kconjugate` coincide.
- Transposition does not exist symbolically. Over a real basis it is `†K`
  with the slots exchanged; over a complex basis it is not a tensor; on
  scalars and bra/ket-less tensors it is the identity; in the IR it is a
  permutation.
- `Constant` conjugates its value. `Variable` and `Power` keep their
  `conjugated_` bit, which is their `*` state (a scalar has no `†` distinct
  from `*`); `adjoint()` and `kconjugate()` on them toggle it and return +1.

### 2.4 What is removed from `Tensor`

`ValueModifier`, `ConjugateModifier`, `TransposeModifier` and their `*`
operators; the array-level `conjugated_`/`transposed_` bits and
`conjugated()`, `transposed()`, `value_modifier()`; `conjugate()`,
`transpose()`, `conjugate_transpose()`, `set_value_modifier()`,
`normalize_value_modifier()`; `value_oriented()` and `ValueOriented`;
`transpose_label` (`conjugate_label` becomes the star mark of section 4);
`even_parity_conjugation_symmetry()` stays only if the painter still needs it
(section 5). `AbstractTensor` is unchanged (it already carries no conjugation
members; `_conjugation_parity()` and `_conjugation_symmetry()` stay).

### 2.5 What `Tensor` keeps or gains

- `bool adjointed() const`, `bool kconjugated() const` (the two states);
  `std::int8_t adjoint() override` (`[[nodiscard]]`, toggles `†`, exchanges
  bra and ket, normalizes, returns the sign; over a real basis applies the
  coset rule of 2.3); `std::int8_t kconjugate() override` (toggles `*`,
  normalizes against the parity, returns the sign); `decorated_label()`
  (`label()` plus the marks as set, `⁺` first). There is no
  `Tensor::conjugate()`: the conjugate of a tensor's value is its adjoint.
- `std::int8_t set_states(bool adjointed, bool kconjugated)`: sets the states
  on a tensor whose slots are already the intended ones, without exchanging
  them; normalizes and returns the sign. Needed by `OpMaker` (which builds the
  adjoint operator's tensor from a `⁺`-marked label with the slots as given),
  by `with_slots`, by the `spin.cpp`/`csv.cpp` rebuild sites and by the
  deserializer.
- Constructors adopt trailing marks in the label, `⁺` and the star in either
  order and at most one of each, into the states (`adopt_marks`, generalizing
  `adopt_adjoint_mark`), normalize, and throw for marks whose normalization
  carries −1 (`sequant::adjoint`/`conjugate(const ExprPtr&)` are the way to
  build minus a tensor). A repeated mark (`t⁺⁺`, two stars) is malformed.
- `with_slots` carries both states and re-normalizes them against the
  rebuilt tensor's traits; a rebuild whose normalization would cost a sign
  throws (the #628 rule).
- `symmetries()` unchanged (traits only).

## 3. Operations on expressions

- The `Expr` interface: two virtuals, `[[nodiscard]] std::int8_t adjoint()`
  and `[[nodiscard]] std::int8_t kconjugate()`, each returning the sign the
  caller must apply. `Tensor`: `adjoint()` toggles `†`, exchanges bra and ket,
  normalizes against the hermiticity; `kconjugate()` toggles `*`, normalizes
  against the parity; queries `adjointed()`, `kconjugated()`. `Variable`,
  `Power`: one flag `conjugated_`, since a scalar's adjoint is its conjugate;
  both operations toggle it and return +1; the query is `conjugated()`.
  `Constant`: no flag; both operations conjugate the value. `Product`, `Sum`:
  `adjoint()` maps to the factors with reversal, `kconjugate()` without; a
  factor's −1 folds into the scalar. `NormalOperator`: `adjoint()` as today,
  `kconjugate()` the identity (2.2).
- `sequant::adjoint(const ExprPtr&)`: unchanged. Every `Expr::adjoint()`
  returns a sign; a −1 is wrapped as `Product{−1, ·}`; `Product::adjoint`
  reverses factors, `Sum::adjoint` maps summands.
- `sequant::conjugate(const ExprPtr&)`: the complex conjugate of the value
  of a c-number expression, implemented as `adjoint(ExprPtr)` (one
  implementation): `Constant`, `Variable`, `Power` conjugate in place, a
  `Tensor` is `adjoint()`ed (sign wrapped), `RealPart`/`ImagPart` are
  invariant, `Sum` maps summands, `Product` conjugates the scalar and adjoints
  every factor (the factor order of a c-number product is immaterial, and the
  result is rebuilt with the constructors' default flattening, so the
  involution holds up to that flattening); operator-valued content throws.
  For a c-number network this is the adjoint of each factor, and the result
  is covariant; for an expression with open indices the head's bra and ket
  are exchanged (`ResultExpr` callers exchange the head).
- `fold_conjugate_pairs`, `is_hermitian_network`: conjugate values through
  `sequant::adjoint`, as they already do; keep.
  Fix A-I1 in the same change: emit `real_part(canon[i])`, and give
  `RealPart`/`ImagPart` a `canonicalize` override and subexpression iteration
  so relabeling reaches the inner expression. Serialization of the two nodes is
  a separate decision (today `serialize` throws on them).
- `sequant::conjugate(const ExprPtr&)`: the value's conjugate, `adjoint` on
  c-number content, a throw on operator-valued content.
- `sequant::kconjugate(const ExprPtr&)`: `kconjugate()` on every `Tensor`
  (sign wrapped), the flag on `Variable`/`Power`, the value on `Constant`,
  the identity on operators, order preserving. On c-number content it is
  well-defined in every basis, since it acts on the operator: over a complex
  basis the surviving `t꙳` atoms are new arrays unless the parity is
  definite. On operator-valued content it requires every index of the
  operators to be over a `K`-closed space (2.2): a real-field space, or a
  space with a registered Kramers partner once #566 provides that registry;
  otherwise it throws `sequant::Exception`.

> **As built.** `sequant::kconjugate(const ExprPtr&)` dispatches to
> `Expr::kconjugate()` on a clone, as `sequant::adjoint(const ExprPtr&)`
> dispatches to `Expr::adjoint()`, rather than carrying a body of its own;
> `Sum` and `Product` recurse through the free function. A summand of a `Sum`
> whose K-conjugate carries a sign stays a `Product{−1, summand}` (a `Sum` has
> no scalar to fold into). The K-closure check walks the atoms of the
> operator-valued expression, sees proto indices, and covers
> `NormalOperator` and `NormalOperatorSequence`; a non-normal-ordered
> `Operator<S>` and an `mbpt::Operator` (whose K-conjugate is itself) pass
> without it.
- `is_hermitian_network` and the `Re`/`Im` algebra in `complex.*`: unchanged.
- The generic `sequant::adjoint(T&&)` template: as in #632 (throws for a sign
  it cannot hold).

## 4. Spelling, parsing, printing

Both states are spelled as trailing marks *in the label*, as `⁺` is today,
on tensors and on variables alike: the deserializer needs no grammar for
them, the constructor adopts them, and the serializer writes
`decorated_label()`. The star mark is `꙳`, U+A673 SLAVONIC ASTERISK
(`sequant::conjugate_label`): Unicode has no spacing superscript asterisk
(the candidates are the combining `⃰` U+20F0, which renders on top of the
preceding glyph, the baseline `∗` U+2217, the small form `﹡` U+FE61 and the
raised cross `˟` U+02DF), and `꙳` is a spacing character that renders as a
star at x-height. It is a single code unit in every `wchar_t` width, so the
trailing-mark adoption is the same `label.back()` test as for `⁺`. Its font
coverage is narrower than `⁺`'s, which does not matter for a serialization
spelling since `to_latex` maps it (below); `∗` is the fallback if it ever
does.

| object | state | serialized | LaTeX |
|---|---|---|---|
| `Tensor` | `†` | `t⁺{…}` | `{t^{\dagger}}` |
| `Tensor` | `*` | `t꙳{…}` | `{t^{*}}` |
| `Tensor` | `†*` | `t⁺꙳{…}` | `{t^{\dagger *}}` |
| `Variable`, `Power` | conjugated (`*`) | `x꙳` | `{x^{*}}` |

Parser rules:

- A tensor or variable name ending in `⁺`, `꙳`, or both (either order) is the
  bare label plus the states; a mark whose normalization carries −1 becomes
  `Product{−1, ·}` at conversion (the #628 rule); a repeated mark is malformed.
- The `^*` and `^T` modifier grammar is removed with `ast::Tensor::modifier`
  and the variable's `^*`; `t^*{…}`, `t^T{…}` and `x^*` are
  `SerializationError`s (the first and last with a message naming the mark). A
  modifier is never seen on an operator name, so the asserts of D-4/E-I3 go
  with the grammar.
- `to_latex` strips the marks from the label and renders them as
  superscripts, `^{\dagger}` for the adjoint and `^{*}` for the star, in that
  order (today `⁺` is passed through as the raw character, which pdfLaTeX
  rejects; the change moves the `{\hat{f⁺}}`-style expectations in
  `test_mbpt.cpp`, listed under the test migration). Notation: this is the
  physics convention, `†` adjoint and `*` conjugation; operator theory's
  `A^*` for the adjoint is not used. For scalars the two coincide, so `z^{*}`
  reads the same either way.
- The star means conjugation of the operator, which over a real basis is the
  elementwise conjugate (the reading #602 gave `^*` there); see "value versus
  operator" in 2.2 for the complex-basis caveat that the mark's documentation
  carries.
- The fourth annotation letter (parity) and the trait letters are unchanged.

> **As built.** The marks are split off the name before the operator and
> reserved-label lookups, so an operator name carrying a mark is a
> `SerializationError` and reserved labels keep their symmetries. A `⁺` on a
> variable name is a `SerializationError` (a variable has one conjugation
> flag; `⁺` has no separate meaning on a scalar). `ast::Variable` holds the
> name as a plain string. A removed spelling (`t^*{…}`, `t^T{…}`, `x^*`) fails
> with the parser's generic positioned error, not with a message naming the
> mark.

`ast::Tensor::modifier` is removed with its grammar branch.

## 5. Hashing, equality, ordering, colouring

- `Tensor::memoizing_hash`: bare label, slots, the two states when set (one
  numeric term, as the modifier contributes today), the symmetry attributes
  including the conjugation symmetry (as today).
- `static_equal`: compares both states and everything it compares today
  minus the removed bits. `static_less_than`: label, then states
  (`t < t⁺ < t꙳ < t⁺꙳`), then the rest as today.
- `hash_terminal_tensor` (eval leaves): bare label, slots, the two states,
  and the conjugation symmetry when it is not `NonSymm` (E-I4). Two leaves
  that `static_equal` distinguishes never share a cache slot.

> **As built.** `hash_terminal_tensor` keys a leaf by the bare label, the
> slot layout, the two states and the conjugation symmetry when it is
> `AntiSymm` (E-I4), and by none of the other traits. The `AntiSymm` term is
> the one that carries value information: over a real basis an odd-parity
> array is imaginary where an even-parity one is real, so two same-label
> leaves of different declared parity must not share a cache slot. `Symm` and
> `NonSymm` add no term, so every other leaf keeps its hash. The term arrived
> in a follow-up batch after the final review, and the ITF fixtures are
> unchanged by it.
- Vertex painter: the tensor core's colour includes the two states and the
  conjugation symmetry when it differs from the even-parity default (as
  today); the array-level `Conjugate`/`Transpose` colour terms go. Bra and ket bundles of a `Conjugate`/`AntiConjugate` tensor are always
  coloured distinctly; only `Symm`/`Antisymm` bundles may be coloured
  interchangeably, and only where the caller consumes the sign
  (`fold_signed_braket`, unchanged).

## 6. Canonicalization

- `DefaultTensorCanonicalizer::canonicalize_braket(t, fold_signed)`: folds
  only the plain exchanges. `Symm` swaps freely, `Antisymm` swaps with −1 into
  the returned sign; `Conjugate` and `AntiConjugate` tensors are left in their
  orientation. The `fold_conjugate` parameter, the unfold-first step, the
  marker-placing swap and the label tie-break for `(Anti)Conjugate` ties go
  (B-I3 becomes moot). Space-tie handling for `Symm`/`Antisymm` is unchanged.
- `braket_foldable(t)` is `braket_swap_sign(t._braket_symmetry()).has_value()`
  for c-number tensors and for operators (whose `_swap_bra_ket` #632 provides);
  `braket_conjugate_foldable` and `as_cnumber_tensor` gates around it go.
- `TensorNetworkV3`: `CreateGraphOptions::fold_conjugate_braket` and
  `CanonicalizeSlotsOptions::fold_conjugate_braket` are removed (B-I7 resolved:
  `canonicalize_slots` is orientation-sensitive for `(Anti)Conjugate` tensors,
  as on master). In `canonicalize_graph` the unfold loop, the
  `(Anti)Conjugate` branch of the graph-dictated reorientation, the
  post-relabel refold of `(Anti)Conjugate` tensors and the second pass for them
  are removed; the `Symm`/`Antisymm` branch, the sign accumulation into the
  parity byproduct, and the whole-bundle swap for column-nonsymmetric tensors
  (#628) stay. The column-symmetry gate and the dummy-edge assertion return to
  gating on `Symm`/`Antisymm` only.

> **As built.** With that gate the strict-braket dummy-edge assertion fires
> for a ket-ket contraction of a `Conjugate` tensor (`C{a_1;p_1} C꙳{a_2;p_1}`):
> a `Conjugate` tensor's bundles are not interchangeable, and the Gram overlap
> is spelled with the adjoint, `C{a_1;p_1} C⁺{p_1;a_2}`. A ket-ket spelling
> inherited from the marker fold fails loudly rather than evaluating to a
> silent value.
- The `TensorBlockCanonicalizer` constructor loses `fold_conjugate_braket`
  and keeps `fold_signed_braket` (off at the eval boundary, per #628).
- Canonical forms: the orientation of a `(Anti)Conjugate` tensor is part of
  its spelling, and two orientations are two values. `t⁺` on a `NonHermitian`
  tensor and `t꙳` on a parity-`None` tensor are distinct atoms (their own
  colour and hash). Over a real basis the coset rule of 2.3 makes `t⁺{q;p}`
  and `t꙳{p;q}` one spelling. The relabeling
  invariance that B-I1 probed is the ordinary graph canonicalization with
  distinct bundle colours; keep its probe as a test.
- `Sum` canonicalization and `simplify` are unchanged; conjugate pairs are
  folded by `fold_conjugate_pairs` (section 3), not by spelling.

## 7. The eval boundary

- `binarize(Tensor)`: an unmarked leaf is a leaf; a `†` (`NonHermitian`)
  leaf is `EvalOp::Adjoint` over the bare leaf, as today (`make_adjoint_node`).
  A `*` leaf (parity `None`) over a real basis is the elementwise conjugate of
  the bare leaf and lowers to `EvalOp::Conjugate` (elementwise, no
  permutation; the op #603 introduces, needed in every backend); over a
  complex basis it names a distinct array the yielder serves. The array-level
  `Conjugate`/`Transpose` cases, the marked-leaf lowering to a value
  orientation with a `Constant(−1)`, and the sign hoist in `binarize(Product)`
  are removed. The `unfolded` slot-counting copy in `binarize(Product)` goes.

> **As built.** No new `EvalOp`: a `꙳` leaf over a real basis lowers to the
> existing `Adjoint` kernel with an identity layout (`canon_ix` equal to the
> leaf's own canonical indices), which is a pure elementwise conjugation. The
> `⁺` case obtains its bare leaf through `Tensor::adjoint()` on a copy
> (bundles exchanged back, state cleared, sign asserted +1) and is served the
> same way in every basis, so a `⁺` that an index-mutating API leaves on a
> real-basis tensor is value-correct; a `꙳` under the `⁺` stays on the operand
> leaf as its own array. The sign hoist and the `unfolded` copy named above
> did not exist at the base; only comments referred to them.
- `EvalExpr(Tensor)`: block canonicalization with `fold_signed_braket = false`
  (the #628 rule; a leaf's phase is a cache-orientation round trip). Leaves keep
  their as-written orientation; the `Symm` swap is the only fold.
- The `Adjoint` op in the backends is unchanged and stays the only
  conjugation the IR performs.
- Storage orientation is out of scope: `t{i;a}` and `t{a;i}` of a Hermitian
  `t` are two leaves, and a yielder serves both. Lowering one as an `Adjoint`
  read of the other's stored block ("lazy conj") is the follow-up the earlier
  specs deferred, and it belongs in `binarize`, not in canonicalization.
- The pre-existing engine hazard (a flat leaf reordered with a sign by the
  block canonicalizer evaluates with the sign dropped) is untouched here.

## 8. Export, mbpt, rewrite rules

- Exporters name a tensor by `decorated_label()`, and neither `⁺` nor `꙳` is
  a legal identifier character in ITF, Julia or Python, so a shared helper
  maps the marks to the suffixes `_adj` and `_conj`, in that order, in every
  generator; today `⁺` is passed through and Python's sanitizer turns it into
  `_`. ITF's import-name map keys by the `Tensor` value, which includes the
  states, and needs no change. D-2 disappears with the array-level states. `ReorderingContext::rewrite` copies the adjointed state via
  `set_adjointed` and drops its "unreachable" comment.

> **As built.** The shared helper is `SeQuant/core/export/marked_name.hpp`
> (`export_label`, `export_name`, `fold_marks_into_label`). The export
> preprocessing folds the marks into the label unconditionally for every
> tensor it sees (label := `export_label(tensor)`, states cleared), before any
> context rewrite and regardless of `enable_rewriting`, so a marked tensor is
> a distinct array to every label-keyed map (declarations, reference counts,
> load strategies, import names) and to the backends' own label matching: the
> ITF two-electron remap never sees a marked integral, and `g⁺` exports as
> `g_adj`, not as a J/K integral. `ReorderingContext::rewrite` bakes the
> suffixes into the rebuilt label instead of copying the states: its rebuild
> has no bra/ket slots, and on such a tensor `set_states` takes the coset
> branch and the default parity then clears the star, which would silently
> drop the mark.
- `mbpt::OpMaker` strips a trailing `⁺` from the operator label, builds the
  tensor with the slots it computed, and applies `set_adjointed(true)`,
  wrapping a −1 in a `Product` (as today through `set_value_modifier`).
  `Operator::adjoint()` unchanged.
- `spin.cpp` and `rules/csv.cpp` rebuild sites copy the adjointed state
  (`with_slots`, or `set_adjointed` where they construct by hand). The
  `value_oriented` calls in `rules/{csv,df,thc}.cpp` go: a `⁺` tensor's slots
  are as written. Whether a rule should match a `⁺` tensor at all (D-3) is
  decided separately; this spec only removes the machinery.

> **As built.** `csv_transform_impl` prepends the sign `set_states()` returns
> to the rebuilt tensor, so it lands on the returned `Product`'s scalar: where
> the csv slots are over a different field, the normalization can trade `⁺`
> for `꙳` and clear it against an odd parity, consuming a sign that a `Tensor`
> cannot hold. Carrying a state across a field change presumes real
> coefficients, which the site's comment states.
- `biorthogonalization.cpp`'s `WK_biorthogonalization_filter_impl` groups
  terms by an orientation-sensitive `canonicalize_slots` hash again, as on
  master.

## 9. Serialization

Unchanged except section 4: both states are trailing label marks, the
`^*`/`^T` tensor grammar is removed, and the fourth annotation letter stays
the `K`-parity.

## 10. Behaviour changes (for the PR description)

1. `Tensor` no longer has `conjugate()`, `transpose()`, the modifier enums or
   `value_oriented()`; `kconjugate()` toggles the operator-conjugation state
   (normalized by the parity); `adjointed()`, `kconjugated()` and
   `set_states()` replace `value_modifier()`/`set_value_modifier()`.
   `sequant::conjugate(ExprPtr)` conjugates the value (the adjoint on c-number
   content) instead of setting an array-level bit.
2. The conjugated state is spelled as a trailing `꙳` in the label, on tensors
   and variables, not `^*`; the `^*`/`^T` grammar is gone. The star means
   conjugation of the operator, which over a real basis is elementwise conjugation
   with the slots in place and over a complex basis is not. `to_latex` renders
   both marks as superscripts. `sequant::conjugate(ExprPtr)` on a slotted tensor
   returns the adjoint spelling (`conj t{i;a}` is `t⁺{a;i}`, which for a
   Hermitian `t` is `t{a;i}`), so expressions that printed `^*` from the fold
   print the exchanged orientation instead (the `θ` and `h` cases in
   `test_mbpt.cpp`).
3. Canonicalization keeps the orientation of a Hermitian (`Conjugate`)
   tensor: the two orientations are two values and canonicalize separately;
   `Symm`/`Antisymm` folds are unchanged. Canonical spellings of expressions
   with complex-basis Hermitian tensors change where the fold used to
   reorient them.
4. `canonicalize_slots` is orientation-sensitive for `Conjugate` tensors
   (`C{a;b}` and `C{b;a}` get different hashes), as on master.
5. `^T` and `^*` are input errors; the marks are the states.
6. The eval leaf hash includes the conjugation symmetry.
   *As built:* the leaf hash carries the two states and the `AntiSymm`
   conjugation symmetry (section 5); `Symm` and `NonSymm` add no term, and
   the ITF fixtures are unchanged by it.
7. `TensorBlockCanonicalizer(bool fold_conjugate_braket, bool
   fold_signed_braket)` becomes `TensorBlockCanonicalizer(bool
   fold_signed_braket)`; `CreateGraphOptions`/`CanonicalizeSlotsOptions` lose
   `fold_conjugate_braket`.
8. The ITF fixtures are expected unchanged (the fold never applied over the
   real-field spaces the drivers use); if item 6 moves a CSE tie-break, that
   is the maintainer's regeneration, not the implementer's.

## 11. Test migration

Cases the experiment shows will change, and what they become:

| test | now asserts | becomes |
|---|---|---|
| `conj_tensor_marker_roundtrip` | `conjugate()` keeps slots, sets a bit | `conjugate(ExprPtr)` exchanges slots and sets the adjointed state; round trip through `⁺` |
| `conjugate_braket_fold_per_tensor` | both orientations fold onto one spelling with a marker | neither orientation moves; `Symm` case keeps its free swap |
| `canonicalize_signed_braket`, three sections | marker placement and one canonical form for two orientations | anti-Hermitian: `d{i;a} u` and `d{a;i} u` canonicalize to distinct forms whose difference is the conjugate pair; Hermitian with default column symmetry: two forms, and `fold_conjugate_pairs` folds their sum |
| `signed_normalization` "odd parity over a real field" | `conjugate()` spells `−t{p;q}` | it spells `t{q;p}`; assert value equality through canonicalization of both |
| `tensor_network_shared` "TN isomorphism", "conjugate braket fold" | `canonicalize_slots` identifies the orientations of a Hermitian `C` | they are distinct; the `Symm` isomorphism stays |
| `test_mbpt.cpp:807`, `:1322` | `^*` spellings | the exchanged-orientation spellings |
| `value_modifier_normalization`, `signed_eval_boundary`, the `^T` round trips, the `Transpose` sections in `test_conjugation.cpp`, the marked-leaf cases in `test_eval_expr.cpp`/`test_eval_btas.cpp`, the `Conjugate`-lowering cases in `test_eval_ta.cpp`/`test_eval_tapp.cpp` | the removed states and lowerings | removed, or rewritten against `⁺` where they pinned a sign |

Roughly 157 references in `test_conjugation.cpp` and 30 to 35 each in
`test_tensor.cpp`, `test_eval_expr.cpp`, `test_mbpt.cpp` and `test_parse.cpp`
touch the removed API; most are mechanical renames. New tests: `conjugate` of
a network is the covariant adjoint spelling; `⁺` on a bra/ket-less tensor is
its conjugate and a `Hermitian` one normalizes it away; `^*` input re-spells;
`^T` input rejected; relabeling invariance of canonical forms with Hermitian
cycles (the B-I1 probe); the leaf-hash separation (E-I4 probe).

## 12. Review items this resolves or leaves

Resolved by construction: B-I3, B-I6, B-I7, D-2, the `^*`/`^T` half of D-4,
E-I4 (section 5), C-I2 (no hoist to test). Fixed alongside: A-I1 (section 3).
Left as they are: D-3 (rule matching on `⁺`), E-I2 (documentation), the
description rewrite, the pre-existing leaf-phase hazard, the hidden `[.]`
failures, the `Antisymm` exploitation in `canonicalize_slots`. #603 and #566
build on the conjugated state: the eval-side `EvalOp::Conjugate`, and time reversal
`Θ = U_T K` with its parity and the Kramers partner map on operators and
indices.

## 13. Relation to #603 and #566

- #566 (`kshitij/feature/kramers-tracing-round2`) declares
  `KramersSymmetry {Nonsymm, TimeReversal}` on `Tensor`, a `Θ`-parity trait
  with only the even value, plus a partner map on index spaces
  (`kramers_partner`, `Index::kramers_flipped`). Both are orthogonal to the
  two core states and thread through the constructors, `with_slots`, hashing,
  equality and the painter the way `column_symmetry_` does. Unchanged by this
  spec.
- #566's single-tensor fold (`canonicalize_kramers`, run after
  `canonicalize_braket` in `TensorBlockCanonicalizer::apply`) flips every
  flavored slot and "toggles the conj marker", i.e. it is written on #602's
  array-level bit. Its identity `T{p̄;q̄} = phase · conj T{p;q}` reads
  `T{p̄;q̄} = phase · T⁺{q;p}` here: flip the flavored slots, exchange bra and
  ket, toggle `†` (normalized by the hermiticity, so `g`, `f` become
  `phase · g{q;p}` and an amplitude becomes a `⁺` atom). Covariant, no marker,
  and served at the eval boundary by the existing `Adjoint` op. The fold's
  shape (a per-tensor pass returning a sign into the phase byproduct) is kept.
- Annotation letters: #628 uses the fourth position for the `K`-parity
  (`E`/`O`/`N`); #566 uses the fourth position for `KramersSymmetry`
  (`T`/`N`). The `Θ`-parity moves to a fifth position, emitted only when not
  defaulted.
- #603 (conjugation in second quantization and eval): `kconjugate` of an
  operator-valued expression conjugates the coefficients and leaves the
  operator string alone (2.2), over any `K`-closed basis; `Θ` over a
  Kramers-paired basis is that composed with the flavor flip and phase. The eval side supplies `EvalOp::Conjugate` for `K` leaves over
  a real basis (section 7). `Θ` on an expression is `kconjugate` composed
  with the flavor flip and phase.

## 14. Documentation

The model has one rule that everything else serves: *`adjoint` is the
Hermitian adjoint, `conjugate` is the complex conjugate of a value (on
c-number content the same operation as `adjoint`, with the slots exchanged),
and `kconjugate` is complex conjugation of the operator (`K O K⁻¹`), which
over a real basis is the elementwise conjugate of the array and over a
complex basis is not.* Three audiences, three documents:

- **Users** (`mbpt`): `adjoint()` on operators, `conjugate` on variables,
  constants and c-number expressions. They meet the marks only in output,
  explained by one table: `t⁺` the adjoint, `t꙳` the conjugated operator
  (over real orbitals, the conjugated array), `x꙳` the conjugate of a scalar.
  `kconjugate` is not part of the user vocabulary.
- **Developers of the symbolic layer**: the class documentation of `Tensor`,
  `Expr::adjoint`, `Expr::conjugate`, and one page under `doc/` holding the
  two states, the two traits with their tables (2.1, 2.2), the rule above,
  the real-basis coset rule and the value-versus-operator paragraph. About a page,
  condensed from section 2.
- **Rationale**: this spec, kept under `doc/dev/specs` as the record of why
  the array-level bits were removed (the type argument, the prior art, the
  experiment). Not required reading for using or maintaining the code.

## 15. Out of scope

Storage orientation at the eval boundary (lazy conj); Kramers-structured
bases and index-level conjugation; complex-symmetric bilinear tensors (two
covariant slots with a permutation symmetry, representable today by slot
structure, not by an exchange symmetry); serialization of `RealPart`/`ImagPart`;
`TensorNetworkV1`/`V2` deprecation.

## 16. Open decisions

1. (decided) Bra/ket-less tensors: the exchange is empty, so `w⁺{;;x}` and
   `w꙳{;;x}` are one value; the coset rule keeps `꙳` there, as over a real
   basis, and `⁺` normalizes into it.
2. (decided) `NormalOperator::kconjugate()` is the identity, in every basis
   (2.2); `kconjugate(ExprPtr)` checks `K`-closure of the operators' index
   spaces (real field now; registered Kramers partner once #566 lands).

Decided: the conjugated state is required (magnetic response, section 2.2);
the real-basis coset rule prefers the star (section 2.3), and bra/ket-less
tensors keep the star too; the state is named "conjugated"; the mark is `꙳`
(U+A673), a trailing label mark like `⁺`, on tensors and variables alike;
scalars have one flag, `conjugated()`, and no `adjointed()`; the operator
conjugation is `kconjugate`, the value conjugation `conjugate`;
`NormalOperator::kconjugate()` is the identity; exporters map
the marks to `_adj`/`_conj` suffixes; `to_latex` renders `⁺` as `^{\dagger}`
and the star as `^{*}`.

As built, in addition: the star normalizes by the parity in every basis
(2.2); `⁺` over a real basis is traded for the star with the slots exchanged
back, and a bra/ket-less `⁺` for the star with nothing to exchange (2.3); the
eval leaf hash carries the states only (5); no eval op beyond `Adjoint` (7);
the exporters fold the marks into the label ahead of every label-keyed step
(8). The page `doc/dev/conjugation-model.md` states the resulting model.
