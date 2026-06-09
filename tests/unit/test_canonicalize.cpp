#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/v1.hpp>
#include <SeQuant/core/tensor_network/v2.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/v3.hpp>
#include <SeQuant/core/utility/tensor.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include "data/sf_r2_direct_real_inc.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

TEST_CASE("canonicalization", "[algorithms]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  auto isr = sequant::mbpt::make_legacy_spaces();
  mbpt::add_pao_spaces(isr);
  auto ctx = get_default_context();
  ctx.set(isr);
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("Tensors") {
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                           Symmetry::Nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_3", L"p_4"},
                           Symmetry::Nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::Nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                           Symmetry::Nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::Symm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}:S"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::Antisymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("-g{p1,p2;p3,p4}:A"));
    }

    // aux indices
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_2"}, aux{L"p_3"},
                           Symmetry::Nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1;p2;p3}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::Nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1,p2;p4,p3;p5}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::Symm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1,p2;p3,p4;p5}:S"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::Antisymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("-B{p1,p2;p3,p4;p5}:A"));
    }
  }

  SECTION("Products") {
    // P.S. ref outputs produced with complete canonicalization
    auto ctx = get_default_context();
    ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
    auto _ = set_scoped_default_context(ctx);

    {
      auto input =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_5"}, ket{L"a_1"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_5", L"a_2"},
                     Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo("Ŝ{a1,a2;i1,i2} f{a3;i3} t{i3;a2} t{i1,i2;a1,a3}"));
    }
    {
      auto input =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo(
              "Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} t{i_2;a_3} t{i_1,i_3;a_1,a_2}"));
    }
    {  // Azam's example:
      // two intermediates that are equivalent modulo permutation of columns of
      // named indices. as of https://github.com/ValeevGroup/SeQuant/pull/349
      // canonicalization does not recognize them as identical because by
      // default named index labels are ignored (this is done to make canonical
      // TNs to be independent of external index renamings). However in
      // the context of a sum external index labels are meaningful and should be
      // accounted.
      for (auto ignore_named_index_labels : {true, false}) {
        auto input1 =
            deserialize(L"1/2 t{a3,a1,a2;i4,i5,i2}:N-C-S g{i4,i5;i3,i1}:N-C-S");
        //      auto input1 = deserialize(L"1/2
        //      t{a1,a2,a3;i5,i2,i4}:N-C-S g{i4,i5;i3,i1}:N-C-S");
        auto input2 =
            deserialize(L"1/2 t{a1,a3,a2;i5,i4,i2}:N-C-S g{i5,i4;i1,i3}:N-C-S");
        canonicalize(
            input1,
            {.method = CanonicalizationMethod::Topological,
             .ignore_named_index_labels =
                 static_cast<CanonicalizeOptions::IgnoreNamedIndexLabel>(
                     ignore_named_index_labels)});
        canonicalize(
            input2,
            {.method = CanonicalizationMethod::Topological,
             .ignore_named_index_labels =
                 static_cast<CanonicalizeOptions::IgnoreNamedIndexLabel>(
                     ignore_named_index_labels)});
        if (ignore_named_index_labels)
          REQUIRE(input1 != input2);
        else
          REQUIRE(input1 == input2);
      }
    }

    {  // Product containing Variables
      auto q2 = ex<Variable>(L"q2");
      q2->adjoint();
      auto input =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, Symmetry::Nonsymm) *
          q2 * ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Variable>(L"p") *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) *
          ex<Variable>(L"q1") *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(input,
                   SimplifiesTo("p q1 q2^* Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "t{i_2;a_3} t{i_1,i_3;a_1,a_2}"));
    }
    {  // Product containing adjoint of a Tensor
      auto f2 = ex<Tensor>(L"f", bra{L"a_1", L"a_2"}, ket{L"i_5", L"i_2"},
                           Symmetry::Nonsymm, BraKetSymmetry::Nonsymm);
      f2->adjoint();
      auto input1 =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) * f2;
      canonicalize(input1);
      REQUIRE_THAT(input1,
                   SimplifiesTo("Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "f⁺{i_1,i_3;a_1,a_2}:N-N-S t{i_2;a_3}"));
      auto input2 =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) * f2 *
          ex<Variable>(L"w") * ex<Constant>(rational{1, 2});
      canonicalize(input2);
      REQUIRE_THAT(input2,
                   SimplifiesTo("1/2 w Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "f⁺{i_1,i_3;a_1,a_2}:N-N-S t{i_2;a_3}"));
    }
    // with aux indices
    {
      auto input =
          ex<Constant>(rational{1, 2}) *
          ex<Tensor>(L"B", bra{L"p_2"}, ket{L"p_4"}, aux{L"p_5"},
                     Symmetry::Nonsymm) *
          ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_3"}, aux{L"p_5"},
                     Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm);
      canonicalize(input);
      // because bra and ket are in same space dummy renaming flips the bra and
      // ket even though the tensors are not bra-ket symmetric
      REQUIRE_THAT(
          input, EquivalentTo("1/2 t{p1;p3} t{p2;p4} B{p3;p1;p5} B{p4;p2;p5}"));
    }
    // with bra-ket symmetry
    {
      // Tensor's BraKetSymmetry is per-tensor (Symm passed explicitly below);
      // no Context manipulation needed.
      // TN is invariant wrt flipping one if the tensors
      // N.B. it's not possible purely to canonicalize each tensor since bra and
      // ket slots are equivalent, only the overall TN topology determines
      // whether bra/ket swap should occur for each tensor
      auto input = ex<Constant>(rational{1, 2}) *
                   ex<Tensor>(L"B", bra{L"p_2"}, ket{L"p_1"}, aux{L"p_5"},
                              Symmetry::Nonsymm, BraKetSymmetry::Symm) *
                   ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_2"}, aux{L"p_5"},
                              Symmetry::Nonsymm, BraKetSymmetry::Symm);
      REQUIRE_THAT(input, EquivalentTo("1/2 B{p1;p2;p5}:N-S B{p1;p2;p5}:N-S"));
    }
    // SF R2 ±pair extracted from the real-field CCSD doubles. Under
    // make_min_sr_spaces + Real-field + Spinfree + SingleProduct (the srcc
    // SF context), the two terms are swap∘column-equivalent for a Symm braket
    // and must collapse to a single 16· term. Before the bundle-vertex fix at
    // v3.cpp (canonical_bra_ket_bundle_order was being read from the last
    // bra/ket slot's canon_perm value instead of from the bra/ket bundle
    // vertex), bliss's input-order-dependent labeling among same-color
    // column-bundle vertices made these two orientations canonicalize to
    // different symbolic forms and not merge under min_sr.
    {
      auto sr_reg = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
      std::vector<std::wstring> keys;
      for (auto const& s : *sr_reg) keys.push_back(s.base_key());
      for (auto const& k : keys)
        if (auto* sp = sr_reg->retrieve_ptr(k)) sp->field(Field::Real);
      // Disable strict bra↔ket-symmetry policy: this expression has a_3 in
      // g.bra and t.bra under one term's orientation (a bra-bra contraction,
      // legitimate for Symm-braket g), which the default-context Conjugate
      // policy would reject.
      auto srcc_resetter = set_scoped_default_context(
          Context({.index_space_registry_shared_ptr = sr_reg,
                   .vacuum = Vacuum::SingleProduct,
                   .spbasis = SPBasis::Spinfree})
              .set(AssertStrictBraKetSymmetry::No));
      auto input = deserialize(
          L"8 * Ŝ{i_1,i_2;a_1,a_2}:N-C-S * g{i_3,a_1;a_3,i_1}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S "
          L"+ 8 * Ŝ{i_1,i_2;a_1,a_2}:N-C-S * g{i_1,a_3;a_1,i_3}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S");
      simplify(input);
      const std::size_t n =
          input ? (input->is<Sum>() ? input->size() : std::size_t{1}) : 0;
      INFO("pair-1 under srcc context → " << n << " terms (1=merged, 2=not)");
      REQUIRE(n == 1);
    }
    // Structural invariant for the same ±pair: build each term as its own
    // TensorNetworkV3 and verify that the canonicalized bliss::Graph objects
    // returned by canonicalize_slots() compare equal via Graph::cmp. This
    // checks the graph encoding (colors + topology) independently of the
    // downstream slot/column/braket reordering applied at v3.cpp:352-411 —
    // if cmp != 0 we'd know the encoding (not the consumer) is at fault.
    {
      auto sr_reg = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
      Context ctx_min = get_default_context();
      ctx_min.set(sr_reg);
      ctx_min.set(AssertStrictBraKetSymmetry::No);
      auto resetter = set_scoped_default_context(ctx_min);
      auto exA = deserialize(
          L"8 * Ŝ{i_1,i_2;a_1,a_2}:N-C-S * g{i_3,a_1;a_3,i_1}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S");
      auto exB = deserialize(
          L"8 * Ŝ{i_1,i_2;a_1,a_2}:N-C-S * g{i_1,a_3;a_1,i_3}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S");
      REQUIRE(exA);
      REQUIRE(exB);
      TensorNetworkV3 tnA(exA);
      TensorNetworkV3 tnB(exB);
      TensorNetworkV3::NamedIndexSet named{Index(L"i_1"), Index(L"i_2"),
                                           Index(L"a_1"), Index(L"a_2")};
      auto mdA = tnA.canonicalize_slots({}, &named);
      auto mdB = tnB.canonicalize_slots({}, &named);
      REQUIRE(mdA.graph);
      REQUIRE(mdB.graph);
      const int cmpAB = mdA.graph->cmp(*mdB.graph);
      INFO("canonical bliss graph cmp(A, B) = " << cmpAB << " (0 = equal)");
      REQUIRE(cmpAB == 0);
    }
    // PNO-CCSD duplicate-intermediate regression (real PAO/PNO/Κ spaces).
    // A real DF factor g(μ̃,μ̃,Κ) transformed PAO->PNO on its bra vs its ket leg
    // yields equivalent half-transformed intermediates that must dedup to one,
    // both as leaf coefficients (C{a;μ̃} ≡ C{μ̃;a}) and as g·C products. Before
    // the fix the eval cache compared stored tensors by bra/ket slot order and
    // saw the two orientations as distinct, recomputing the (large)
    // intermediate.
    {
      auto sr_reg = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
      mbpt::add_pao_spaces(sr_reg);  // μ̃ (PAO)
      mbpt::add_df_spaces(sr_reg);   // Κ  (DF aux)
      std::vector<std::wstring> keys;
      for (auto const& s : *sr_reg) keys.push_back(s.base_key());
      for (auto const& k : keys)
        if (auto* sp = sr_reg->retrieve_ptr(k)) sp->field(Field::Real);
      Context ctx = get_default_context();
      ctx.set(sr_reg);
      ctx.set(AssertStrictBraKetSymmetry::No);
      auto resetter = set_scoped_default_context(ctx);

      auto graph_cmp = [](std::wstring a, std::wstring b) {
        auto exA = deserialize(a), exB = deserialize(b);
        REQUIRE(exA);
        REQUIRE(exB);
        TensorNetworkV3 tnA(exA), tnB(exB);
        auto mdA = tnA.canonicalize_slots();
        auto mdB = tnB.canonicalize_slots();
        REQUIRE(mdA.graph);
        REQUIRE(mdB.graph);
        return mdA.graph->cmp(*mdB.graph);
      };
      auto nodes_equal = [](std::wstring a, std::wstring b) {
        auto na = binarize(deserialize(a)), nb = binarize(deserialize(b));
        TreeNodeEqualityComparator<std::remove_cvref_t<decltype(na)>> eq;
        // sanity: the eval-node hashes fold (the comparator must be consistent)
        CHECK(na->hash_value() == nb->hash_value());
        return eq(na, nb);
      };

      // the proto-indexed leaf coefficient folds under bra<->ket swap ...
      CHECK(graph_cmp(L"C{a_1<i_1,i_2>;μ̃_1}:N-S-S",
                      L"C{μ̃_1;a_1<i_1,i_2>}:N-S-S") == 0);
      // ... and so does the g·C product network ...
      CHECK(graph_cmp(L"g{μ̃_1;μ̃_2;Κ_1}:N-S-S * C{a_1<i_1,i_2>;μ̃_1}:N-S-S",
                      L"g{μ̃_1;μ̃_2;Κ_1}:N-S-S * C{μ̃_2;a_1<i_1,i_2>}:N-S-S") ==
            0);

      // ... and, crucially, the bra- vs ket-transform g·C eval nodes compare
      // EQUAL (so the cache deduplicates them): both when the surviving PNO
      // external is the same (a_1) and when it differs but shares the space and
      // pair domain (a_1 vs a_4 — residual target vs internal dummy in CCSD).
      CHECK(nodes_equal(L"g{μ̃_1;μ̃_2;Κ_1}:N-S-S * C{a_1<i_1,i_2>;μ̃_1}:N-S-S",
                        L"g{μ̃_1;μ̃_2;Κ_1}:N-S-S * C{μ̃_2;a_1<i_1,i_2>}:N-S-S"));
      CHECK(nodes_equal(L"g{μ̃_1;μ̃_2;Κ_1}:N-S-S * C{a_1<i_1,i_2>;μ̃_1}:N-S-S",
                        L"g{μ̃_1;μ̃_2;Κ_1}:N-S-S * C{μ̃_2;a_4<i_1,i_2>}:N-S-S"));
    }
    // Top-level regression: the 113-term SF R2 direct-path (real field)
    // residual extracted byte-for-byte from `srcc 2 t std sf real`. The
    // spin-traced reference collapses to 110 terms; before the fix the direct
    // path produced 113 and refused to merge the 3 ±pairs that exercised the
    // bra↔ket-swap canonicalization under bliss's ambiguous labeling of
    // same-color bundle vertices.
    {
      auto sr_reg = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
      std::vector<std::wstring> keys;
      for (auto const& s : *sr_reg) keys.push_back(s.base_key());
      for (auto const& k : keys)
        if (auto* sp = sr_reg->retrieve_ptr(k)) sp->field(Field::Real);
      auto srcc_resetter = set_scoped_default_context(
          Context({.index_space_registry_shared_ptr = sr_reg,
                   .vacuum = Vacuum::SingleProduct,
                   .spbasis = SPBasis::Spinfree})
              .set(AssertStrictBraKetSymmetry::No));
      auto input = deserialize(tests::data::sf_r2_direct_real());
      REQUIRE(input);
      const std::size_t n_before =
          input->is<Sum>() ? input->size() : std::size_t{1};
      REQUIRE(n_before == 113);
      simplify(input);
      const std::size_t n_after =
          input ? (input->is<Sum>() ? input->size() : std::size_t{1}) : 0;
      INFO("after simplify: " << n_after
                              << " terms (expected 110 if merge works)");
      REQUIRE(n_after == 110);
    }
  }

  SECTION("Sum of Variables") {
    {
      auto input =
          ex<Variable>(L"q1") + ex<Variable>(L"q1") + ex<Variable>(L"q2");
      simplify(input);
      canonicalize(input);
      REQUIRE_THAT(input, EquivalentTo("q2 + 2 q1"));
    }

    {
      auto input =
          ex<Variable>(L"q1") * ex<Variable>(L"q1") + ex<Variable>(L"q2");
      simplify(input);
      canonicalize(input);
      REQUIRE_THAT(input, EquivalentTo("q2 + q1 * q1"));
    }
  }

  SECTION("Powers in Products") {
    const auto f = deserialize(L"f{p_1;p_2}:A-C-S * ã{p_2;p_1}");
    const auto t1 = deserialize(L"t{a_1;i_1}:A-C-S * ã{i_1;a_1}");
    const auto pw = ex<Power>(ex<Variable>(L"x"), rational{1, 2});

    auto expr1 = f * t1 * pw;
    simplify(expr1);
    REQUIRE_THAT(expr1, EquivalentTo(L"x^(1/2) * ã{p_2;p_1} * ã{i_1;a_1} * "
                                     L"t{a_1;i_1}:A-C-S * f{p_1;p_2}:A-C-S"));

    auto expr2 = ex<Constant>(rational{1, 2}) * f * t1 * pw;
    simplify(expr2);
    REQUIRE_THAT(expr2, EquivalentTo(L"1/2 x^(1/2) * ã{p_2;p_1} * ã{i_1;a_1} * "
                                     L"t{a_1;i_1}:A-C-S * f{p_1;p_2}:A-C-S"));

    auto expr3 = f * t1 * ex<Power>(2, 3);
    simplify(expr3);
    REQUIRE_THAT(expr3, EquivalentTo(L"8 * ã{p_2;p_1} * ã{i_1;a_1} * "
                                     L"t{a_1;i_1}:A-C-S * f{p_1;p_2}:A-C-S"));
  }

  SECTION("Sum of Powers") {
    const auto vx = ex<Variable>(L"x");

    // x^{1/2} + x^{1/2} = 2 * x^{1/2}
    auto pw1 = ex<Power>(vx, rational{1, 2});
    auto pw2 = ex<Power>(vx, rational{1, 2});
    auto sum_expr = pw1 + pw2;
    simplify(sum_expr);
    REQUIRE(sum_expr->is<Product>());
    REQUIRE(sum_expr->as<Product>().scalar() == 2);

    // 2 x^{1/2} + 3 x^{1/2} = 5 x^{1/2}
    auto s1 = ex<Constant>(2) * ex<Power>(vx, rational{1, 2});
    auto s2 = ex<Constant>(3) * ex<Power>(vx, rational{1, 2});
    auto sum2 = s1 + s2;
    simplify(sum2);
    REQUIRE(sum2->is<Product>());
    REQUIRE(sum2->as<Product>().scalar() == 5);
  }

  SECTION("Sum of Products") {
    // P.S. ref outputs produced with complete canonicalization
    auto ctx = get_default_context();
    ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
    auto _ = set_scoped_default_context(ctx);

    {
      // CASE 1: Non-symmetric tensors
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm);
      simplify(input);
      canonicalize(input);
      REQUIRE_THAT(input,
                   EquivalentTo("g{p_1,p_2;p_3,p_4} t{p_3;p_1} t{p_4;p_2}"));
    }

    // CASE 2: Symmetric tensors
    {
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::Symm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::Symm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(input, EquivalentTo("g{p2,p3;p1,p4}:S t{p1;p2} t{p4;p3}"));
    }

    // Case 3: Anti-symmetric tensors
    {
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(input, EquivalentTo("g{p2,p3;p1,p4}:A t{p1;p2} t{p4;p3}"));
    }

    // Case 4: permuted indices
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"i_1"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::Antisymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         Symmetry::Antisymm);
      canonicalize(input);
      REQUIRE(input->size() == 1);
      REQUIRE_THAT(input,
                   EquivalentTo("g{i3,i4;i1,a3}:A t{a2;i3} t{a1,a3;i2,i4}:A"));
    }

    // Case 4: permuted indices from CCSD R2 biorthogonal configuration
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"i_1"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::Nonsymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         Symmetry::Nonsymm);

      canonicalize(input);
      REQUIRE(input->size() == 1);
      REQUIRE_THAT(input,
                   EquivalentTo("g{i3,i4;i1,a3} t{a2;i4} t{a1,a3;i3,i2}"));
    }

    {  // Case 5: CCSDT R3: S3 * F * T3

      {  // Terms 1 and 6 from spin-traced result
        auto input =
            ex<Constant>(-4) *
                ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_2", L"i_4"}, Symmetry::Nonsymm) +
            ex<Constant>(-4) *
                ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_4", L"i_3"}, Symmetry::Nonsymm);
        canonicalize(input);
        REQUIRE_THAT(
            input,
            EquivalentTo(
                "-8 Ŝ{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2}"));
      }

      {
        auto term1 =
            ex<Constant>(-4) *
            ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_3", L"i_2", L"i_4"}, Symmetry::Nonsymm);
        auto term2 =
            ex<Constant>(-4) *
            ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_2", L"i_4", L"i_3"}, Symmetry::Nonsymm);
        canonicalize(term1);
        canonicalize(term2);
        REQUIRE_THAT(term1,
                     EquivalentTo("-4 Ŝ{i_1,i_2,i_3;a_1,a_2,a_3} f{i_4;i_3} "
                                  "t{a_1,a_2,a_3;i_1,i_4,i_2}"));
        REQUIRE_THAT(term2,
                     EquivalentTo("-4 Ŝ{i_1,i_2,i_3;a_1,a_2,a_3} f{i_4;i_3} "
                                  "t{a_1,a_2,a_3;i_1,i_4,i_2}"));
        auto sum_of_terms = term1 + term2;
        simplify(sum_of_terms);
        REQUIRE_THAT(
            sum_of_terms,
            EquivalentTo(
                "-8 Ŝ{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2}"));
      }

      {  // Terms 2 and 4 from spin-traced result
        auto input =
            ex<Constant>(2) *
                ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_4", L"i_2"}, Symmetry::Nonsymm) +
            ex<Constant>(2) *
                ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Nonsymm);
        canonicalize(input);
        REQUIRE_THAT(
            input, EquivalentTo(
                       "4 Ŝ{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i4,i1,i2}"));
      }
    }

    // Case 6: Case 4 w/ aux indices
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"B", bra{L"i_3"}, ket{L"a_3"}, aux{L"p_5"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"B", bra{L"i_4"}, ket{L"i_1"}, aux{L"p_5"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::Nonsymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"B", bra{L"i_3"}, ket{L"i_1"}, aux{L"p_5"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"B", bra{L"i_4"}, ket{L"a_3"}, aux{L"p_5"},
                         Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, Symmetry::Nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         Symmetry::Nonsymm);

      canonicalize(input);
      simplify(input);
      REQUIRE_THAT(
          input,
          EquivalentTo("t{a2;i4} t{a1,a3;i3,i2} B{i3;i1;p5} B{i4;a3;p5}"));
    }
  }
}

TEST_CASE("braket_symmetric_half_tensor_canonicalization", "[algorithms]") {
  using namespace sequant;

  auto canon_hash = [](std::wstring spec) {
    auto e = deserialize(spec);
    ExprPtrList tl{e};
    TensorNetwork tn(tl);
    return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels())
        .hash_value();
  };

  // A braket-symmetric tensor is invariant under bra<->ket exchange, so a
  // half-tensor with the orbital in bra must canonicalize identically to the
  // form with the orbital in ket (regression: the vertex painter's tensor
  // "shade" hashed bra_rank/ket_rank in fixed order, distinguishing them).
  CHECK(canon_hash(L"X{a1;;i1}:N-S-N") == canon_hash(L"X{;a1;i1}:N-S-N"));
  // Without braket symmetry the two forms must remain distinct.
  CHECK(canon_hash(L"X{a1;;i1}:N-N-N") != canon_hash(L"X{;a1;i1}:N-N-N"));
}
