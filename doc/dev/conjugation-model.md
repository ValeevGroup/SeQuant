# The conjugation model: adjointed and K-conjugated tensors

A `Tensor` is `<bra|O|ket>`: the matrix, over the slots' basis, of a core
`O` with two traits, `Hermiticity` and `ConjugationParity`. Two commuting Z2
states act on the core, not on the array, so both are covariant in every
basis and for every slot structure (rectangular and half tensors, aux slots,
bra/ket-less tensors alike):

| state        | mark | denotes                                                                  | via            |
|--------------|------|--------------------------------------------------------------------------|----------------|
| adjointed    | `⁺`  | `t⁺{q;p} = conj t{p;q}`, the matrix of `O†`                              | `adjoint()`    |
| K-conjugated | `꙳`  | `t꙳{p;q} = <p|K O K⁻¹|q>`, the conjugated operator's matrix, slots in place | `kconjugate()` |

`K` is complex conjugation in the coordinate representation, the operation
`ConjugationParity` is defined against; it is not time reversal
(`Θ = U_T K`, the Kramers work). The states form the Klein four-group
`{1, ⁺, ꙳, ⁺꙳}` and live on `Tensor` only: operator-valued tensors
(`NormalOperator`) have an adjoint and are K-invariant.

## Normalization by the traits

After every mutation (construction from a marked label, `adjoint()`,
`kconjugate()`, `set_states()`, `with_slots()`) each state is normalized
against its trait, so a set state always denotes a genuinely distinct array.
The sign consumed is returned to the caller; `sequant::adjoint` /
`kconjugate(const ExprPtr&)` turn it into a scalar factor, and a constructor
or `with_slots()` refuses a mark whose normalization carries −1.

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
spellings of one value. `normalize_states()` enters this branch on the
hermiticity alone (a `⁺` the hermiticity does not clear) and trades it for a
`꙳` with the bundles exchanged back; the parity then normalizes that star.
So `adjoint()` on such a tensor leaves the slots in place: over real orbitals
`adjoint(t{a;i})` is `t꙳{a;i}` for parity `None` and `t{a;i}` for the default
parity (a K-even operator has a real matrix there); over complex orbitals it
is `t⁺{i;a}`.

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
denote one value, and spell it alike only where the coset branch runs
(hermiticity indefinite): for a Hermitian `h{i;a}` `conjugate` gives
`h{a;i}` and `kconjugate` gives `h{i;a}`, equal through the declared
symmetry. `fold_conjugate_pairs` and `is_hermitian_network` conjugate values
through `adjoint`.

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
| `Variable` | `꙳`   | `x꙳`           | `{x^{*}}`         |               |
| `Power`    | `꙳`   | `(x^(2))^*`    | `{{{x}^2}^{*}}`   |               |

The deserializer splits the marks off a tensor or variable name and applies
them through `set_states()`, so a mark whose normalization carries −1 becomes
`Product{−1, ·}`; a repeated mark, a `⁺` on a variable name and any mark on an
operator name are `SerializationError`s, and `^*`/`^T` are not grammar on
tensor and variable names (`^*` is the `Power` spelling). An `mbpt::Operator`
carries the `⁺` in its own label, which its LaTeX keeps as written
(`{\hat{t⁺}}`); the tensor form's amplitude prints as `{t^{\dagger}}`.

## The eval boundary

`EvalExpr(Tensor)` block-canonicalizes a leaf with `fold_signed_braket` off:
indices are reordered within each bundle, the `Symm` bra/ket swap is the only
exchange, and a leaf keeps its as-written orientation. `binarize(Tensor)`
serves the states as IR ops over the bare leaf, which carries no sign:

- `⁺`: `EvalOp::Adjoint` (permute bra/ket, then conjugate) over the bare
  array, in every basis; a `꙳` under it stays on the operand leaf.
- `꙳` over a real basis: the Adjoint kernel with an identity layout, the
  elementwise conjugate of the bare array.
- `꙳` over a complex basis: a plain leaf the yielder serves as its own array.

The leaf hash keys an array by label, slot layout and states, not by the
traits; an Adjoint node hashes the bare leaf salted by the op. The Adjoint op
is the only conjugation the IR performs.

## Export and mbpt

The exporters fold the marks into the label before any context rewrite
(`fold_marks_into_label`, suffixes `_adj` then `_conj`), so a marked tensor is
an array of its own to every label-keyed map and backend match: `g⁺` is
exported as `g_adj` and is not remapped to a J/K integral; the index
reordering bakes the suffixes into the rebuilt label.

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
  to keep it. A marked tensor is a separate array, named `t_adj` / `t_conj`
  by the exporters and served by the yielder as described above.
- Over real orbitals `adjoint` and `conjugate` spell a star with the slots in
  place; over complex orbitals, `t⁺` with bra and ket exchanged.
- The canonicalizer keeps the orientation of a Hermitian (`Conjugate`)
  tensor: `t{i;a}` and `t{a;i}` are two values whose sum folds through
  `fold_conjugate_pairs`, not by spelling.
