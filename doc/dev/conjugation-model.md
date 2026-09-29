# The conjugation model: adjointed and K-conjugated tensors

A `Tensor` is `<bra|O|ket>`: the matrix, over the slots' basis, of a core
`O` with two traits, `Hermiticity` and `ConjugationParity`. Two commuting Z2
states act on the core, not on the array, so both are covariant in every
basis and for every slot structure (rectangular and half tensors, aux slots,
bra/ket-less tensors alike):

| state        | mark | denotes                                                                  | via            |
|--------------|------|--------------------------------------------------------------------------|----------------|
| adjointed    | `⁺`  | `t⁺{q;p} = conj t{p;q}`, the matrix of `O†`                              | `adjoint()`    |
| K-conjugated | `꙳`  | `t꙳{p;q} = ⟨p\|K O K⁻¹\|q⟩`, the conjugated operator's matrix, slots in place | `kconjugate()` |

`K` is complex conjugation in the coordinate representation, the operation
`ConjugationParity` is defined against; it is not time reversal
(`Θ = U_T K`, the Kramers work). The states form the Klein four-group
`{1, ⁺, ꙳, ⁺꙳}` and live on `Tensor` only: operator-valued tensors
(`NormalOperator`) have an adjoint and are K-invariant.

## Normalization by the traits

Normalization against the traits runs at construction, in `set_states()`,
`adjoint()`, `kconjugate()`, `set_label()` / `adopt_marks()`, `with_slots()`,
and the slot-mutating `transform_indices()`, `set_bra()`, `set_ket()` and
`set_aux()`, so a state set through any of those always denotes a genuinely
distinct array. A slot mutation can put the tensor on a basis of another field,
where the coset rule identifies `⁺` with `꙳` and a definite trait consumes the
mark outright, so a slot mutation derives `BraKetSymmetry` and the elementwise
`ConjugationSymmetry` again from the traits and the new slots' field and then
normalizes the states against them: a mutated tensor is the tensor a
construction over those slots would have given. The `AbstractTensor`
primitives `_swap_bra_ket()` and `_bra_mutable()` / `_ket_mutable()` reconcile
nothing; normalization and the tensor-network machinery drive them themselves.
The sign consumed is returned to the caller; `sequant::adjoint` /
`kconjugate(const ExprPtr&)` turn it into a scalar factor, while a
constructor, `with_slots()` and a slot mutation refuse a mark whose
normalization carries −1 — a refused slot mutation leaves the tensor as the
call found it, down to the slots' field.

| `Hermiticity`   | `t⁺`             | kept? | sign returned by `adjoint()` |
|-----------------|------------------|-------|------------------------------|
| `Hermitian`     | `t`              | no    | +1                           |
| `AntiHermitian` | `−t`             | no    | −1                           |
| `NonHermitian`  | a distinct array | yes   | +1                           |

| `ConjugationParity` | `t꙳`             | kept? | sign returned by `kconjugate()` |
|---------------------|------------------|-------|---------------------------------|
| `Even`              | `t`              | no    | +1                              |
| `Odd`               | `−t`             | no    | −1                              |
| `None`              | a distinct array | yes   | +1                              |

The parity is a property of the operator (`K O K⁻¹ = ±O`), so `꙳` normalizes
against it in every basis. The elementwise `ConjugationSymmetry` is a
different, derived observable (the parity read through the basis field; a
bra/ket-less tensor is `NonSymm` in every field), consumed only by the
tensor's hash and equality and by the vertex painter's colour. Which slot
exchanges are symmetries is `BraKetSymmetry`'s business, derived from the
hermiticity, the parity and the field. A bra/ket-less tensor (`w{;;x}`) has
an empty exchange, so its adjoint is its conjugate: there the hermiticity
reads "real", "imaginary", "unknown" and normalizes `꙳` ahead of the parity.

The defaults are `NonHermitian` and `Even` (a real operator): `t⁺` on a
default tensor is kept, `t꙳` on a default tensor is `t` in every basis, and a
kept star requires `ConjugationParity::None` (serialized as the fourth
annotation letter, `t꙳{a_1;i_1}:N-N-N-N`).

## Value versus operator

- The complex conjugate of the value `t{p;q}` is `t⁺{q;p}`, in every basis.
- `t꙳{p;q}` is the operator conjugated. Over a real basis (`Kp = p`) it is the
  conjugated value with the slots in place; over a complex basis it is the
  matrix of a different operator and not the conjugated array. A star written
  on a complex-basis tensor to mean "conjugate the array" denotes the wrong
  object, and the deserializer cannot tell.
- Transposition is not a state: over a real basis it is `⁺꙳` with the slots
  exchanged, over a complex basis it is not a tensor, in the IR it is a
  layout permutation.

## The coset rule

Over a real basis, or with no bra/ket slot, `t⁺{q;p}` and `t꙳{p;q}` are two
spellings of one value, and `normalize_states()` identifies them in both
directions. Where the hermiticity is indefinite (a `⁺` it does not clear) the
`⁺` is traded for a `꙳` with the bundles exchanged back; the parity then
normalizes that star. So `adjoint()` on such a tensor leaves the slots in
place: over real orbitals `adjoint(t{a;i})` is `t꙳{a;i}` for parity `None` and
`t{a;i}` for the default parity (a K-even operator has a real matrix there);
over complex orbitals it is `t⁺{i;a}`.

Where the hermiticity is definite it wins, ahead of the parity: a `꙳` is
`t⁺{q;p}`, which the hermiticity reduces to `±t{q;p}`, so the mark is consumed
at that sign and the bundles are exchanged. Over real orbitals
`kconjugate(h{p_1;p_2})` of a Hermitian, parity-`None` `h` is `h{p_2;p_1}` and
`kconjugate(d{p_1;p_2})` of an anti-Hermitian `d` is `−d{p_2;p_1}` — the same
atoms `conjugate` produces, so one value has one spelling.

## Vocabulary

| operation      | on a `Tensor`                                                                  | on a scalar (`Constant`, `Variable`, `Power`)                          | on operator-valued content                                                                                                                     |
|----------------|--------------------------------------------------------------------------------|------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| `adjoint(e)`   | exchanges bra and ket, toggles `⁺`, normalizes (coset rule over a real basis)  | the complex conjugate; `Variable`/`Power` keep one `conjugated()` flag | the adjoint; a `Product` reverses its factors                                                                                                  |
| `conjugate(e)` | the conjugate of the value: `adjoint`, one implementation                      | the same                                                               | throws `sequant::Exception` (an operator has no value)                                                                                         |
| `kconjugate(e)`| toggles `꙳`, slots in place, normalizes by the parity; factor order kept       | the same                                                               | the identity on the string; throws if a `NormalOperator`/`NormalOperatorSequence` index (proto indices included) is over a space `K` does not close, i.e. not a real-field space; `Operator<S>` and `mbpt::Operator` are not checked |

Each `Expr::adjoint()` / `Expr::kconjugate()` returns the sign it consumed;
the free functions wrap a −1 as a scalar factor and rebuild `Sum`s and
`Product`s with the constructors' default flattening. `conjugate` of an
expression with open indices exchanges the head's bra and ket, so a
`ResultExpr` caller must exchange the head itself. On scalars the three
coincide. On c-number content over a real basis `conjugate` and `kconjugate`
denote one value and spell it alike: a definite hermiticity consumes the star
and exchanges the bundles (for a Hermitian `h{i;a}` both give `h{a;i}`), an
indefinite one keeps the star on the slots as written (both give `t꙳{a;i}`).
`fold_conjugate_pairs` and `is_hermitian_network` conjugate values through
`adjoint`.

## Spelling

Both states are trailing marks of the label, in either order and at most one
of each; `label()` is the bare array name, `decorated_label()` the printed
one, and hashing, equality, ordering and graph colouring combine the bare
label with the states. The star is `꙳`, U+A673 (`sequant::conjugate_label`).

| object     | state | serialized     | LaTeX             | exported name |
|------------|-------|----------------|-------------------|---------------|
| `Tensor`   | `⁺`   | `t⁺{…}`        | `{t^{\dagger}}`   | `t_adj`       |
| `Tensor`   | `꙳`   | `t꙳{…}`        | `{t^{*}}`         | `t_conj`      |
| `Tensor`   | `⁺꙳`  | `t⁺꙳{…}`       | `{t^{\dagger *}}` | `t_adj_conj`  |
| `Variable` | `꙳`   | `x꙳`           | `{{x}^{*}}`       |               |
| `Power`    | `conjugated()` flag, no mark | `(x^(2))^*` | `{{{x}^{2}}^{*}}` |    |

The deserializer splits the marks off a tensor or variable name and applies
them through `set_states()`, so a mark whose normalization carries −1 becomes
`Product{−1, ·}`; a repeated mark, a `⁺` on a variable name and any mark on an
operator name are `SerializationError`s, and `^*`/`^T` are not grammar on
tensor and variable names (`^*` is the `Power` spelling). An `mbpt::Operator`
carries the `⁺` in its own label, whose LaTeX prints the bare label with a
dagger, braced so that a rank superscript attaches to the group
(`{{\hat{t}^{\dagger}}^{1}}`, or `{\hat{f}^{\dagger}}` where no rank is
printed); the tensor form's amplitude prints as `{t^{\dagger}}`.

## The eval boundary

`EvalExpr(Tensor)` respells a leaf as the array a provider serves and records
what maps that array to the value the leaf denotes: a `CanonTransform`
`{phase, conj, braket_swap}`, applied on retrieval and excluded from the leaf's
own slot hash. `decode_leaf_states` is exhaustive over the two states in three
cases:

- `⁺`: `t⁺{q;p}` is `conj t{p;q}`, so the array is the bare `t` -- the
  bundles are exchanged back and the state cleared, which costs no sign -- and
  the transform is `{conj, braket_swap}`.
- `꙳` over a real basis: the elementwise conjugate of the same array with
  the slots in place, so the array is the bare `t` and the transform is
  `{conj}` over the identity layout.
- `꙳` over a complex basis: the matrix of another operator, an array of its
  own, so the state stays on the stored spelling and the transform is trivial.

`t⁺꙳` over a complex basis composes the first case with the third: the
array is `t꙳`, the transform `{conj, braket_swap}`. The decoded spelling is
then block-canonicalized with `fold_signed_braket` off -- indices are
reordered within each bundle, the `Symm` bra/ket swap is the only exchange, and
a leaf keeps its as-written orientation -- and that sign composes into the same
transform, which is applied once on the way out of the leaf fetch. A
tensor-of-tensors leaf, one with proto indices, runs that same normalization
and then, in addition, takes its slot hash and a further phase from
`TensorNetwork::canonicalize_slots`. `denoted_expr()` inverts the
decoding, respelling the stored array as written. `binarize(Tensor)` is
therefore a plain leaf in every case.

A transform hoists out of a product or a sum exactly when it is a pure
conjugation (`hoistable`: `conj` without `braket_swap`), because elementwise
conjugation distributes over contraction and addition while a bra/ket exchange
respells the node's own result. A transform carrying an exchange salts the
parent's hash instead, so the adjoint of a contraction keeps a slot of its own.

The leaf hash keys an array by bare label, slot layout, the states left on the
stored spelling and an `AntiSymm` conjugation symmetry, not by the other
traits. It is label-blind, so a flat Hermitian tensor's two orientations share
one slot and each node's annotations carry the wiring; a tensor-of-tensors
leaf's hash is its canonical labeling, which colours a `Conjugate` tensor's
bundles apart, so the same Hermitian pair hashes apart there -- a missed cache
hit, never a wrong value.

## Export and mbpt

The exporters fold the marks into the label before any context rewrite
(`fold_marks_into_label`, suffixes `_adj` then `_conj`), so a marked tensor is
an array of its own to every label-keyed map and backend match: `g⁺` is
exported as `g_adj` and is not remapped to a J/K integral; the index
reordering bakes the suffixes into the rebuilt label.

A conjugated `Variable` is spelled `conj(x)` by the generators that have a
conjugation spelling for a scalar (text, Julia, Python-einsum). ITF has none,
so there the same folding names it: a conjugated `x` is the object `x_conj[]`
wherever it appears -- declaration, load, value and drop alike.

`adjoint` of an `mbpt::Operator` marks its label `t⁺`, inverts its
quantum-number action and regenerates its tensor form as
`sequant::adjoint(tensor form)`, which is where an anti-Hermitian operator's
−1 comes from. Separately, `OpMaker` given an already-marked name (an
operator registered as `v⁺`) strips the mark, builds the tensor with the
slots it computed and applies `set_states(true, false)`, wrapping a −1 in a
`Product`. The `spin` and `csv` rebuild sites carry both states through
`with_slots()` / `set_states()`; a sign consumed by a rebuild across a field
change lands on the term's scalar.

## What a user meets

- `t꙳` on a default-parity tensor is `t`; declare `ConjugationParity::None`
  to keep it. A marked tensor is an array of its own to the exporters, named
  `t_adj` / `t_conj`; at evaluation its states decode into the retrieval
  transform as described above.
- Over real orbitals `adjoint` and `conjugate` spell a star with the slots in
  place; over complex orbitals, `t⁺` with bra and ket exchanged.
- The canonicalizer keeps the orientation of a Hermitian (`Conjugate`)
  tensor: `t{i;a}` and `t{a;i}` are two values whose sum folds through
  `fold_conjugate_pairs`, not by spelling.
