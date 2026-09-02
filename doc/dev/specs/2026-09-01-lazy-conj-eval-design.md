#Lazy - conj eval redesign(conjugation PR 2) — design

Date: 2026-09-01. Builds on the symbolic conjugation redesign
(`kshitij/feature/conjugation-symbolic`, PR #602, base commit
`13dd902c7`). Draft PR for evaluation;
branch
`kshitij / feature / conjugation - eval`,
    PR based on the PR -
        1 branch and retargeted to master when #602 merges.

        ##Problem

            PR 1 made the elementwise -
        conjugation marker first -
        class in the symbolic layer but kept the eval layer on the structural
            model
    : conjugation is served by an `EvalOp::Adjoint` IR
      node(own cache slot,
           own dispatch stage `NeedLeftAdj`, `Constant{1}` sentinel,
           eager whole - array
`Result::adjoint`),
    while the canonicalization * phase *
        is already served lazily on
        retrieval(`apply_phase`, cache holds the pre - phase value)
            .Consequences,
    confirmed by review of the PR -
        1 branch :

    -the eval boundary disables the conjugate
    -
    braket fold,
    so leaf serving is orientation - sensitive
    : a starred Conjugate spelling and the plain swapped spelling of the same
      value are served differently(review finding F1, CONFIRMED),
    and off - pipeline orientation aliasing is possible(F3);
- `fold_conjugate_braket` is threaded through five gating points that must stay
    mutually consistent(F4);
- marked -
    Nonsymm leaves cannot be served at all — PR 1 added explicit
  `std::logic_error` gates in `binarize(Tensor)` (
        including before the '⁺' label channel)pending this redesign;
- every adjoint -
    served leaf pays a separate cache slot plus an eager
  `.conj()` materialization per node.

    ##Architecture:one retrieval -
                       time canonical transform

                       Generalize the phase
                       - on - retrieval model to carry conjugation and bra <
                       ->ket
transposition:

    - `EvalExpr::canon_phase_` generalizes to a canonical transform
  `canon_transform_ =
        {phase : int8_t, conj : bool, braket_swap : bool}` — the map from the *
            cached canonical * result to *
            this node's denoted* value.Excluded from the
            node's own hash, exactly like phase today. - `apply_phase` in `evaluate` (
                                                             eval.hpp)becomes
  `apply_canon_transform`.Backend
hook:one new virtual
  `Result::apply_transform(phase, conj, swap - annotation)`;
the TA implementation fuses relabel + `.conj()` + scale in one pass;
the default composes `mult_by_phase` with the existing
        pieces.- `EvalOp::Adjoint`,
    the `NeedLeftAdj` stage, the `Constant{1}` sentinel convention,
    and `make_adjoint_node` are retired.

            Leaf boundary(`EvalExpr(Tensor)` / `binarize(Tensor)`)
    :

      -**Conjugate
            * *leaves : the braket fold is re
        - enabled at the eval boundary(`fold_conjugate_braket` returns to one
                  setting everywhere — F4 dissolves)
              .Both orientations land on one cache slot spelled canonically;
the folded spelling carries `{ conj, braket_swap }`.
  Yielders are only ever asked for the canonical spelling — F1 and F3
  dissolve.
- **Nonsymm starred** leaves (`t^*`): same slot as `t`, transform
  `{
  conj
}
`.The PR - 1 throws in `binarize` are replaced by serving.- **'⁺' -
    labeled Nonsymm adjoints** : label stripped at the boundary,
    same slot as bare `t`, transform `{
  conj, braket_swap
}
`.The '⁺' + marker combination composes to `{ braket_swap }` alone (pure
  transpose) — that refusal also becomes serving.
- **Symm** markers: dropped (identity in value), as today.

`sequant::value_oriented` is unchanged: its Nonsymm throw remains
correct for its surviving (symbolic, slot-rebuilding) consumers — the
mbpt rules;
the eval - layer call sites are replaced by transform -
    based logic.

        ##Cache /
        CSE identity and the hash contract

            Identity splits in two;
this is the load - bearing invariant.

    - **Slot identity**(what the cache stores under)
    : a node's own hash covers only its canonical spelling;
its own transform is excluded.So
  `t`/`t ^ *` and both orientations of a Conjugate tensor share one slot
                    holding the canonical value
                        .PR 1's Nonsymm marker salt in
  `hash_terminal_tensor` moves out of the slot hash into the transform.-
                **Structural identity**(what a parent is)
    : a child's *phase* is multiplicatively hoistable(`A·(-B) = -(A·B)`) and
    stays out of hash combination,
folding into the parent's own transform — unchanged. Conjugation obeys a finer
    rule — *
    *uniform conj hoists,
mixed conj salts * * — because elementwise conj distributes over contraction and
    addition(`(A·B) ^ * = A ^ *·B ^ *`, `(sum_i T_i) ^ * = sum_i T_i ^ *`,
  `(c·X) ^ * = conj(c)·X ^ *`)
    : -if EVERY tensor child of a product / sum carries the conj
      bit(and the scalar is conjugated correspondingly),
the conjugation is a whole - node transform : strip the child marks,
    conjugate the scalar, set `{
  conj
}` on the node's own transform, and hash the UNMARKED
    spelling. Recursively, whole-subtree conjugation bubbles to the
    root, so a conjugated TN hashes onto its unconjugated
    counterpart's slot: `B^*·A^*` (reordering already canonicalized
    away) is a cache HIT on the `A·B` slot, served as one retrieval
    conj — no recomputation. This is the mechanism the TRS/Kramers
    tracing fold (PR 3) relies on to serve \mathcal{T}-partner
    (time-reversal) intermediates
    from one evaluation;
it is a named requirement of this PR,
    not an optimization.- with MIXED marks the conjugation is NOT a whole -
        node transform(`C·C ^ *` is no transform of `C·C`),
    so each marked child contributes `slot_hash(+)
            transform_salt(conj, braket_swap)` to the parent
        's hash — what the retired Adjoint node' s
    `hash::combine(h, EvalOp::Adjoint)` encoded structurally.The C·C
        ^ *vs C ^
        *·C regression tests must pass unchanged.
  `braket_swap` does not hoist in
         this PR(adjoint of a product also reverses the contraction frame);
hermitian - network recognition stays with the symbolic layer's `is_hermitian_network()`.

    Consequences :

    -One stored canonical array per slot;
each consumer applies its own transform on retrieval.TA has only
    eager `.conj()`,
    so a retrieval materializes a transformed copy — the same per
        - use cost as the Adjoint node today,
    but without duplicate cache entries.A TA lazy conj
        - view is a named future optimization,
    out of scope.-
        PR 1's `conjugated_tensors` report(test - only until now)
            gains its intended consumer : deriving a ToT leaf's transform from
  `canonicalize_slots`.

                                              **Conj
        -
        aware CSE(in scope).*
            *The hoisting rule applies recursively at EVERY node,
    so common - subexpression identity is conj - aware at all granularities,
    not just whole terms :

    -inside a mixed product,
    a uniformly conjugated subgroup hoists at its own node
    : binarizing `A
      ^
      *·B ^
      *·C` forms the `(A ^ *·B ^ *)` intermediate with hash
          == the canonical `A·B` slot and transform `{
  conj
}
` — a cache HIT when `A·B` (or its reordering) was computed anywhere,
    in this term or another;
- `CacheManager` keys and the intermediate hashes (`imed_hashes`)
  inherit this from the hash contract, so cross-term CSE spans conj
  variants with no separate machinery — PR 1's deferred symbolic
  rewrite `A^*B^*·C -> (BA)^*·C` is realized as slot identity at
  binarize time rather than as an expression transform (the symbolic
  named hook stays for pretty-printing/export);
- the single -
    term optimizer's repeated-subnet detection is already fold -
    aware(audit item 5),
    so its groupings and the eval - side slots agree; steering the contraction-ORDER search toward conj-reusable
  groupings (a cost-model credit for conj-related cached intermediates)
  is a stretch goal within this PR, exercised by the mixed-product
  test below.

## Symbolic-surface audit -> eval obligations

Every conj construct PR 1 introduced, checked against this design
(gap audit, 2026-09-01):

1. **Re/Im nodes (MAJOR).** `simplify()` in a complex field folds
   `A + A^* -> 2 Re(A)` and `A - A^* -> 2 i Im(A)` behind `SimplifyOptions::FoldConjugatePairs`,
   whose default is parked at `No` explicitly "until the evaluation
   layer understands them" (options.hpp); `binarize(ExprPtr)` throws
   `Exception("unsupported expression")` on a RealPart/ImagPart node.
   PR 2 adds first-class unary eval nodes for Re/Im — they are
   PROJECTIONS (not invertible), so unlike Adjoint they remain IR
   nodes: `EvalOp::RealPart`/`ImagPart`, backend
   `Result::real_part()/imag_part()` (TA: elementwise), inner operand
   served from its own (shared)
cache slot.PR 2 then flips the FoldConjugatePairs default to `Yes`,
    completing the owner's A6 decision(auto - fold in complex -
                                       field simplify())
            .2. *
        *Conjugated scalar leaves(`Variable`, `Power`).**Export already prints
                                                        them(`wrap_conj`),
    but eval cannot serve them
    : a marked Variable
      /
      Power hashes to its own slot and no conj is ever applied.Under
      this design they follow the same transform model as tensor leaves
    : slot = unmarked spelling,
      transform = `{
  conj
}`
   (`conj(b^n) = conj(b)^n` for integer Power), and
   `Result::apply_transform` must therefore work for scalar-backed
   results too.
3. **Anti-conjugate forward-compat (PR 3).** The plan records the
   owner directive to add a sign-carrying fold
   (`T{
  p;
  q} = -conj(T{
  q;
  p})`, Hermiticity::AntiHermitian) when TRS
   lands. The transform accommodates it as `{
phase:
  -1, conj, braket_swap
}
` with no structural change — noted so nothing in PR 2 precludes it.4. *
    *Export.**Today `EvalOp::Adjoint` exports only the transpose,
    with a documented real - field -
        only limitation(conj is not representable in the exported IR)
            .With Adjoint retired,
    export consumes the canonical transform instead : `braket_swap` as the same index reordering, `conj` via
    the generators' existing `wrap_conj` where supported,
    else the same documented real
        - field limitation — no regression against today.5. *
              *Optimizer consistency.*
                  * `single_term_detail` subnet identity already
                  runs `canonicalize_slots` with the fold ON,
    so conj - related subnets dedupe there;
once the eval boundary re - folds(F4 dissolution),
    optimizer subnet identity and eval slot identity share one fold semantics by
        construction.6. *
        *Reality recognition **(`N = conj(N)` up to canonicalization = > real,
                                the Kramers energy fold basis) stays
        symbolic(`is_hermitian_network()`);
eval needs nothing beyond conj hoisting; covered by a worked example.

## Testing (owner requirements, recorded in the PR-1 plan)

- Unit tests for evaluation WITH conj.
- Verification + worked examples of each PR-1 symbolic relation showing
  how it is exploited at eval time — per relation: the symbolic form,
  the eval strategy, a numeric test.
- **Exhaustive over PR 1's symbolic catalogue**: a new
  `test_eval_conjugation.cpp` mirrors EVERY test case of PR 1's
  `test_conjugation.cpp` (and the conj sections of `test_eval_expr`)
  with an eval-level numeric validation, maintained as a one-to-one
  mapping table in the test file; a purely-symbolic invariant with no
  eval observable (e.g. serialization roundtrip, with_slots attribute
  carriage, TN slot determinism) gets an explicit `n/a` entry with
  one-line justification instead of silent omission. Mapping classes:
*marker / involution / roundtrip cases->slot -
    identity invariants(`(
        A ^ *)^*` shares `A`'s slot; marked-leaf transforms compose to identity);
*per - symmetry fold identities(`T {
         q;
         p
       } = conj(T {
                                  p;
                                  q
                                })`)
           ->numeric equality of served values on random Hermitian data,
    per BraKetSymmetry class;
* `adjoint_conjugate_transpose_relations` (Klein four - group)->transform
    - composition tests : '⁺' alone,
    marker alone, '⁺' + marker(pure transpose),
    each served value checked against the directly - computed reference;
*2Re / 2iIm fold cases(incl.mixed sums, bystanders, i - rotation,
                       scalar hoisting) -> fold
        - on
    == fold - off numeric equality;
* `hermitian_network_recognition`/ reality->the energy - reality worked example;
* `conjugate_free_function_total` ->every node kind conj -
    evaluated(Constant, Variable, Power, Tensor, Sum, Product, Re / Im).-
    Existing conj identity
    regressions(cache ON == cache OFF, C·C ^ *vs C ^ *·C) pass unchanged.-
    Conj - hoisting cache reuse : evaluate `A·B`,
    then request
  `conjugate(A·B)` (spelled `B ^ *·A ^ *`) — assert a cache hit on the
  `A·B` slot and `result == conj(A·B)` numerically; same for a sum of
  uniformly conjugated products (the Kramers-tracing shape).
- Conj-aware CSE at intermediate granularity: with `A·B` cached,
  evaluate `A^*·B^*·C` — assert the `(A^*·B^*)` intermediate is a
  cache hit on the `A·B` slot and the result equals the reference
  computed without reuse.
- Re/Im evaluation: `Re(A) + i·Im(A) == A` on random complex data;
`2·Re(A)` from the auto - fold evaluates equal to `A + A ^
    *` computed with the fold off;
inner `A` slot shared between `Re(A)` and other consumers of `A`.-
            Conjugated scalar leaves
    : `x ^
        *` (Variable) and `conj(b ^ n)` (Power)
                                  served as conj of the unmarked leaf's value. -
                              Energy -
                              reality worked example
    : a self
      -
      conjugate scalar network evaluates with zero imaginary part(numeric),
    matching
  `is_hermitian_network()` recognition(symbolic).- Export / serialized eval
        - tree goldens churn(Adjoint nodes vanish from tree shapes);
regenerated with justification in the PR description.-
    MPQC smoke before the draft goes up
    : on the current MPQC branch(`kshitij / refactor / kramer - pairs -
                                 separation`),
    repin SeQuant to the PR - 2 branch and re -
        run the certified decks(h2o - 0.23967684425, dch - 1.04169026053)
            .

        ##Out of scope

        - TA lazy conj -
        view(retrieval stays eager per use).- TRS / korbit folding(PR 3).

## Implementation deltas (recorded 2026-09-01, Plan A complete)

Refinements the implementation forced on the letter (not the spirit) of
the contract above; the plan doc's "Execution deviations" section has
the play-by-play.

- **Marker composition is syntactic and slot-free**: the elementwise
  marker contributes a PURE {conj} bit; orientation deltas come from the
  canonicalizer fold alone. (The draft's marker-unfold would have
  re-aliased starred-canonical spellings with their plain form.)
- **Hoisting is per-PREFIX of the left fold**: a uniformly conjugated
  factor prefix strips its factors' conj salts (per-prefix, keeping
  C.C^* vs C^*.C order-insensitive) and records {conj} on its node --
  this is what makes the buried-intermediate CSE work; whole-product
  uniformity is the special case.
- **Evaluation invariant**: hand-ups are DENOTED values; the cache holds
  CANONICAL data. The leaf evaluator serves the canonical spelling and
  the leaf hand-up applies the transform; finish_phase_b's store
  re-applies it (an involution), leaving canonical data in the cache.
- **ToT leaves carry array-faithful Nested (outer;inner) indices** via
  tot_indices -- the md named list (proto-only constituents, named
  order) stays internal to slot canonicalization.
- **Export**: scalar leaves re-materialize the conj marker
  (denoted_expr); tensor leaves keep the documented transpose-only
  real-field limitation.
- **symmetrize/antisymmetrize wrappers** use the DENOTED bra rank
  (braket_swap-aware).
- **Test-fixture contract**: yield keys and stored arrays are
  canonical-spelling shaped; literal spellings are served through the
  leaf's transform (the fixture caches the transformed variant).
