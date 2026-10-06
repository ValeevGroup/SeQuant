#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>  // mbpt::Spin

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
#include <string_view>
#include <type_traits>
#include <vector>

// the `particle_symmetric` symmetry pack (column = Symm) is defined in
// catch2_sequant.hpp and shared across the MBPT test TUs
using namespace sequant::tests;

TEST_CASE("canonicalization", "[algorithms]") {
  using namespace sequant;

  auto isr = sequant::mbpt::make_legacy_spaces();
  mbpt::add_pao_spaces(isr, mbpt::Spin::null);
  auto ctx = get_default_context();
  ctx.set(isr);
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("Tensors") {
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                           particle_symmetric);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_3", L"p_4"},
                           particle_symmetric);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           particle_symmetric);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                           particle_symmetric);
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
                           particle_symmetric);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1;p2;p3}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, particle_symmetric);
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
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_5"}, ket{L"a_1"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_5", L"a_2"},
                     particle_symmetric);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo("Ŝ{a1,a2;i1,i2} f{a3;i3} t{i3;a2} t{i1,i2;a1,a3}"));
    }
    {
      auto input =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     particle_symmetric);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo(
              "Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} t{i_2;a_3} t{i_1,i_3;a_1,a_2}"));
    }
    {  // Azam's example:
      // two spellings of one product, which differ by permutations of the
      // columns of its column-symmetric tensors; they canonicalize alike
      // whether or not named index labels are ignored
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
        REQUIRE(input1 == input2);
      }
    }

    {  // Product containing Variables
      auto q2 = ex<Variable>(L"q2");
      q2->adjoint();
      auto input =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          q2 * ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Variable>(L"p") *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, particle_symmetric) *
          ex<Variable>(L"q1") *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     particle_symmetric);
      canonicalize(input);
      REQUIRE_THAT(input,
                   SimplifiesTo("p q1 q2^* Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "t{i_2;a_3} t{i_1,i_3;a_1,a_2}"));
    }
    {  // Product containing adjoint of a Tensor
      auto f2 = ex<Tensor>(L"f", bra{L"a_1", L"a_2"}, ket{L"i_5", L"i_2"},
                           Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                           ColumnSymmetry::Symm);
      f2->adjoint();
      auto input1 =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, particle_symmetric) * f2;
      canonicalize(input1);
      REQUIRE_THAT(input1,
                   SimplifiesTo("Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "f⁺{i_1,i_3;a_1,a_2}:N-N-S t{i_2;a_3}"));
      auto input2 =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, particle_symmetric) * f2 *
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
                     particle_symmetric) *
          ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_3"}, aux{L"p_5"},
                     particle_symmetric) *
          ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric);
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
                              Symmetry::Nonsymm, BraKetSymmetry::Symm,
                              ColumnSymmetry::Symm) *
                   ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_2"}, aux{L"p_5"},
                              Symmetry::Nonsymm, BraKetSymmetry::Symm,
                              ColumnSymmetry::Symm);
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
          L"8 * Ŝ{i_1,i_2;a_1,a_2} * g{i_3,a_1;a_3,i_1}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S "
          L"+ 8 * Ŝ{i_1,i_2;a_1,a_2} * g{i_1,a_3;a_1,i_3}:N-S-S "
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
          L"8 * Ŝ{i_1,i_2;a_1,a_2} * g{i_3,a_1;a_3,i_1}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S");
      auto exB = deserialize(
          L"8 * Ŝ{i_1,i_2;a_1,a_2} * g{i_1,a_3;a_1,i_3}:N-S-S "
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
      mbpt::add_pao_spaces(sr_reg, mbpt::Spin::any);  // μ̃ (PAO)
      mbpt::add_df_spaces(sr_reg);                    // Κ  (DF aux)
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
        auto exA = deserialize(a), exB = deserialize(b);
        REQUIRE(exA);
        REQUIRE(exB);
        // binarize(ExprPtr) is deprecated (builds a positional head); here we
        // only need the eval node to compare, so suppress as other eval tests
        // do
        SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
        auto na = binarize(exA);
        auto nb = binarize(exB);
        SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
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
                         particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric);
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
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::Symm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric);
      canonicalize(input);
      REQUIRE_THAT(input, EquivalentTo("g{p2,p3;p1,p4}:S t{p1;p2} t{p4;p3}"));
    }

    // Case 3: Anti-symmetric tensors
    {
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, particle_symmetric);
      canonicalize(input);
      REQUIRE_THAT(input, EquivalentTo("g{p2,p3;p1,p4}:A t{p1;p2} t{p4;p3}"));
    }

    // Case 4: permuted indices
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"i_1"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::Antisymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                         Symmetry::Antisymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, particle_symmetric) *
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
                         particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         particle_symmetric) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                         particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         particle_symmetric);

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
                           ket{L"a_1", L"a_2", L"a_3"}, particle_symmetric) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}, particle_symmetric) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_2", L"i_4"}, particle_symmetric) +
            ex<Constant>(-4) *
                ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, particle_symmetric) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}, particle_symmetric) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_4", L"i_3"}, particle_symmetric);
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
                       ket{L"a_1", L"a_2", L"a_3"}, particle_symmetric) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}, particle_symmetric) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_3", L"i_2", L"i_4"}, particle_symmetric);
        auto term2 =
            ex<Constant>(-4) *
            ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, particle_symmetric) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}, particle_symmetric) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_2", L"i_4", L"i_3"}, particle_symmetric);
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
                           ket{L"a_1", L"a_2", L"a_3"}, particle_symmetric) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}, particle_symmetric) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_4", L"i_2"}, particle_symmetric) +
            ex<Constant>(2) *
                ex<Tensor>(reserved::symm_label(), bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, particle_symmetric) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}, particle_symmetric) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_3", L"i_4"}, particle_symmetric);
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
                         particle_symmetric) *
              ex<Tensor>(L"B", bra{L"i_4"}, ket{L"i_1"}, aux{L"p_5"},
                         particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         particle_symmetric) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"B", bra{L"i_3"}, ket{L"i_1"}, aux{L"p_5"},
                         particle_symmetric) *
              ex<Tensor>(L"B", bra{L"i_4"}, ket{L"a_3"}, aux{L"p_5"},
                         particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, particle_symmetric) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         particle_symmetric);

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
    return tn.canonicalize_slots(get_default_context().cardinal_tensor_labels())
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

TEST_CASE("context_tensor_canonicalizers", "[algorithms]") {
  using namespace sequant;

  auto Q = [](std::initializer_list<std::wstring_view> b,
              std::initializer_list<std::wstring_view> k) {
    return ex<Tensor>(L"Q", bra(b), ket(k), Symmetry::Antisymm);
  };

  SECTION("Tensor::canonicalize uses the context's default canonicalizer") {
    {
      auto scoped = set_scoped_default_context(
          Context(get_default_context())
              .set_tensor_canonicalizer(
                  L"", std::make_shared<NullTensorCanonicalizer>()));
      auto q = Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"});
      CHECK(q->as<Tensor>().canonicalize() == nullptr);
      CHECK(*q == *Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"}));
    }
    auto q = Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"});
    const auto bp = q->as<Tensor>().canonicalize();
    REQUIRE(bp);
    CHECK(bp->as<Constant>().value<int>() == -1);
    CHECK(*q == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
  }

  SECTION("TN canonicalization uses the context's label canonicalizer") {
    auto canonicalized = [&Q] {
      TensorNetwork tn(ExprPtrList{Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"})});
      const auto byproduct =
          tn.canonicalize(get_default_context().cardinal_tensor_labels(),
                          {.method = CanonicalizationMethod::Complete});
      const int phase = byproduct ? byproduct->as<Constant>().value<int>() : 1;
      return std::make_pair(std::dynamic_pointer_cast<Expr>(tn.tensors().at(0)),
                            phase);
    };
    {
      // a canonicalizer registered for label Q, which leaves Q as is and
      // counts its uses
      struct CountingCanonicalizer : NullTensorCanonicalizer {
        std::shared_ptr<int> calls = std::make_shared<int>(0);
        ExprPtr apply(AbstractTensor&) const override {
          ++*calls;
          return {};
        }
      };
      const auto counting = std::make_shared<CountingCanonicalizer>();
      auto scoped = set_scoped_default_context(
          Context(get_default_context())
              .set_tensor_canonicalizer(L"Q", counting));
      const auto [tensor, phase] = canonicalized();
      CHECK(*counting->calls > 0);
      // the graph-based pass alone places the named indices in label order
      CHECK(*tensor == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
      CHECK(phase == -1);
    }
    const auto [tensor, phase] = canonicalized();
    CHECK(*tensor == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
    CHECK(phase == -1);
  }

  SECTION("DefaultTensorCanonicalizer uses the context's index comparer") {
    const auto default_cmp = TensorCanonicalizer::default_index_comparer();
    {
      auto scoped = set_scoped_default_context(
          Context(get_default_context())
              .set_index_comparer(
                  [default_cmp](const Index& idx1, const Index& idx2) {
                    return default_cmp(idx2, idx1);
                  }));
      auto q = Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"});
      q->as<Tensor>().canonicalize();
      CHECK(*q == *Q({L"i_2", L"i_1"}, {L"a_2", L"a_1"}));
    }
    auto q = Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"});
    CHECK(q->as<Tensor>().canonicalize() == nullptr);
    CHECK(*q == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
  }

  SECTION("Product canonicalization follows the context's cardinal labels") {
    auto canonical_labels = [] {
      auto product = ex<Tensor>(L"Z", bra{L"i_1"}, ket{L"a_1"}) *
                     ex<Tensor>(L"Y", bra{L"a_1"}, ket{L"i_1"});
      canonicalize(product, {.method = CanonicalizationMethod::Complete});
      REQUIRE(product->is<Product>());
      std::vector<std::wstring> labels;
      for (const auto& factor : product->as<Product>().factors())
        labels.emplace_back(factor->as<Tensor>().label());
      return labels;
    };
    const std::vector<std::wstring> y_first{L"Y", L"Z"};
    const std::vector<std::wstring> z_first{L"Z", L"Y"};
    CHECK(canonical_labels() == y_first);
    {
      auto scoped = set_scoped_default_context(
          Context(get_default_context()).set_cardinal_tensor_labels({L"Z"}));
      CHECK(canonical_labels() == z_first);
    }
    CHECK(canonical_labels() == y_first);
  }
}

TEST_CASE("canonicalization_zero_by_symmetry", "[algorithms][canonicalize]") {
  using namespace sequant;

  auto isr = sequant::mbpt::make_legacy_spaces();
  auto ctx = get_default_context();
  ctx.set(isr);
  auto ctx_resetter = set_scoped_default_context(ctx);

  auto is_zero = [](const ExprPtr& expr) {
    return expr->is<Constant>() && expr->as<Constant>().is_zero();
  };

  // a1<->a2 is an automorphism; it permutes the (symmetric) bra of t and
  // the (antisymmetric) creators of the operator, so the term equals minus
  // itself
  SECTION("symmetric tensor contracted with fermionic operator") {
    auto make = [](Symmetry symm) {
      return ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}, symm) *
             ex<FNOperator>(cre({L"a_1", L"a_2"}), ann({L"i_1", L"i_2"}));
    };
    auto zero = make(Symmetry::Symm);
    canonicalize(zero);
    REQUIRE(zero->is_zero());
    // a zero product keeps no factors, so zeros of different inputs agree
    REQUIRE(zero->as<Product>().factors().empty());
    {
      ExprPtr other_zero =
          ex<Tensor>(L"t", bra{L"a_2", L"a_1"}, ket{L"i_2", L"i_1"},
                     Symmetry::Symm) *
          ex<FNOperator>(cre({L"a_2", L"a_1"}), ann({L"i_1", L"i_2"}));
      canonicalize(other_zero);
      REQUIRE(*other_zero == *zero);
      REQUIRE(other_zero->hash_value() == zero->hash_value());
    }
    simplify(zero);
    REQUIRE(is_zero(zero));

    auto nonzero = make(Symmetry::Antisymm);
    simplify(nonzero);
    REQUIRE(!nonzero->is_zero());
  }

  SECTION("symmetric tensor contracted with antisymmetric tensor") {
    // dummy indices only
    {
      auto expr = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"h", bra{L"p_3", L"p_4"}, ket{L"p_1", L"p_2"},
                             Symmetry::Symm);
      simplify(expr);
      REQUIRE(is_zero(expr));
    }
    // named indices are not permuted, the contracted pair is
    {
      auto make = [] {
        return ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                          Symmetry::Antisymm) *
               ex<Tensor>(L"s", bra{L"p_5", L"p_6"}, ket{L"p_1", L"p_2"},
                          Symmetry::Symm);
      };
      auto expr = make();
      simplify(expr);
      REQUIRE(is_zero(expr));

      // in a sum named index labels are meaningful
      auto sum = make() + ex<Tensor>(L"f", bra{L"p_3", L"p_4"},
                                     ket{L"p_5", L"p_6"}, Symmetry::Antisymm);
      simplify(sum);
      REQUIRE_THAT(sum, EquivalentTo("f{p3,p4;p5,p6}:A"));
    }
  }

  SECTION("nonzero terms are unchanged") {
    // swapping identical antisymmetric tensors (with their slots) is even
    {
      auto expr = ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                             Symmetry::Antisymm);
      simplify(expr);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo("t{a1,a2;i1,i2}:A t{i1,i2;a1,a2}:A"));
    }
    // CCD quadratic term: nontrivial automorphisms, all of them even
    {
      auto expr = ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"t", bra{L"a_3", L"a_4"}, ket{L"i_1", L"i_2"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_3", L"i_4"},
                             Symmetry::Antisymm);
      simplify(expr);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo("g{i3,i4;a3,a4}:A t{a3,a4;i1,i2}:A "
                                      "t{a1,a2;i3,i4}:A"));
    }
    // a lone antisymmetric tensor: swapping its named indices is not a
    // symmetry of the term
    {
      auto expr = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                             Symmetry::Antisymm);
      simplify(expr);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo("g{p1,p2;p3,p4}:A"));
    }
    // automorphisms permuting antisymmetric bundles of unequal or odd size
    // in pairs are even
    for (std::wstring spec :
         {L"t{a1,a2;i1}:A u{;a1,a2}:A", L"X{a1,a2,a3;i1}:A Y{;a1,a2,a3}:A"}) {
      auto expr = deserialize(spec);
      simplify(expr);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo(spec));
    }
    {
      auto sum = deserialize(L"t{a1,a2;i1}:A u{;a1,a2}:A + w{;i1}:A");
      simplify(sum);
      REQUIRE(sum->is<Sum>());
      REQUIRE(sum->size() == 2);
    }
  }

  // a bra<->ket symmetric tensor's bra bundle can map onto its ket bundle:
  // p1<->p3, p2<->p4 maps g onto itself (+1) and is odd on u
  SECTION("automorphism exchanging the bra and ket of a tensor") {
    auto make = [](Symmetry w_symm) {
      return ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, aux{},
                        Symmetry::Antisymm, BraKetSymmetry::Symm) *
             ex<Tensor>(L"u", bra{L"p_1", L"p_3"}, ket{}, Symmetry::Antisymm) *
             ex<Tensor>(L"w", bra{L"p_2", L"p_4"}, ket{}, w_symm);
    };
    auto zero = make(Symmetry::Symm);
    simplify(zero);
    REQUIRE(is_zero(zero));

    // w antisymmetric too: the exchange is even
    auto nonzero = make(Symmetry::Antisymm);
    simplify(nonzero);
    REQUIRE(!nonzero->is_zero());
  }

  // only topological canonicalization looks for automorphisms
  SECTION("rapid canonicalization does not detect a zero") {
    auto expr = ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                           Symmetry::Symm) *
                ex<FNOperator>(cre({L"a_1", L"a_2"}), ann({L"i_1", L"i_2"}));
    expr->canonicalize(CanonicalizeOptions::default_options().copy_and_set(
        CanonicalizationMethod::Rapid));
    REQUIRE(!expr->is_zero());
    canonicalize(expr);
    REQUIRE(expr->is_zero());
  }

  // odd-size bundles: antisymmetric in a1,a2,a3 against symmetric
  SECTION("odd-size antisymmetric bundle contracted with symmetric one") {
    auto expr = deserialize(L"X{a1,a2,a3;i1}:A Y{;a1,a2,a3}:S");
    simplify(expr);
    REQUIRE(is_zero(expr));
  }
}

TEST_CASE("canonicalize_named_index_automorphisms", "[algorithms]") {
  using namespace sequant;

  // Topological canonicalization with named index labels ignored: automorphisms
  // of the network that permute named indices leave their placement to the
  // labels, so that every spelling of an expression canonicalizes alike and
  // canonicalizing again changes nothing (#666)
  const CanonicalizeOptions opts{
      .method = CanonicalizationMethod::Topological,
      .ignore_named_index_labels =
          CanonicalizeOptions::IgnoreNamedIndexLabel::Yes};
  // pairs of spellings of one expression
  for (const auto& [x, y] :
       std::initializer_list<std::pair<std::wstring_view, std::wstring_view>>{
           // identical tensors
           {L"X{i_1;a_1} * X{i_2;a_2}", L"X{i_2;a_2} * X{i_1;a_1}"},
           {L"X{i_1;a_1} * X{i_2;a_2} * Y{i_3;a_3}",
            L"Y{i_3;a_3} * X{i_2;a_2} * X{i_1;a_1}"},
           // antisymmetric slots, dummies renamed
           {L"g{i_1,i_2;a_3,a_4}:A * t{a_3,a_4;i_3,i_4}:A",
            L"t{a_5,a_6;i_3,i_4}:A * g{i_1,i_2;a_5,a_6}:A"},
           {L"f{i_1;a_3} * t{a_3,a_1;i_2,i_3}:A",
            L"-1 f{i_1;a_5} * t{a_1,a_5;i_2,i_3}:A"},
           // columns of a column-symmetric tensor
           {L"X{i_3,i_5;a_1,a_2}:N-N-S", L"X{i_5,i_3;a_2,a_1}:N-N-S"},
           {L"X{i_3,i_5,i_6;a_1,a_2,a_4}:A-N-S",
            L"X{i_6,i_3,i_5;a_4,a_1,a_2}:A-N-S"}}) {
    auto canonicalized = [&opts](std::wstring_view input) {
      ExprPtr e = deserialize(input);
      // a lone tensor is canonicalized as a network only within a Product
      if (!e->is<Product>()) e = ex<Product>(ExprPtrList{e});
      canonicalize(e, opts);
      return e;
    };
    const auto cx = canonicalized(x);
    CAPTURE(x, y, serialize(cx));
    REQUIRE(cx == canonicalized(y));
    auto again = cx->clone();
    canonicalize(again, opts);
    REQUIRE(again == cx);
  }
}

TEST_CASE("current_contexts_version", "[algorithms]") {
  using namespace sequant;
  const auto version = current_contexts_version();

  // any change to the canonicalization configuration of the context in
  // effect for any Statistics changes the version, for as long as it is in
  // effect
  auto changed_by = [&](auto&& modify, Statistics s = Statistics::Arbitrary) {
    Context ctx(get_default_context(s));
    modify(ctx);
    auto resetter = set_scoped_default_context({{s, ctx}});
    return current_contexts_version() != version;
  };
  CHECK(
      changed_by([](Context& ctx) { ctx.set_cardinal_tensor_labels({L"Z"}); }));
  CHECK(changed_by([](Context& ctx) {
    ctx.set_tensor_canonicalizer(L"version_test",
                                 std::make_shared<NullTensorCanonicalizer>());
  }));
  CHECK(changed_by([](Context& ctx) {
    ctx.set_index_comparer(TensorCanonicalizer::default_index_comparer());
  }));
  for (auto s : {Statistics::FermiDirac, Statistics::BoseEinstein,
                 Statistics::Arbitrary}) {
    CHECK(changed_by(
        [](Context& ctx) {
          ctx.set_cardinal_tensor_labels(container::vector<std::wstring>{L"Z"});
        },
        s));
  }
  // settings that canonicalization does not read do not change it
  CHECK(!changed_by([](Context& ctx) {
    ctx.set(ctx.vacuum() == Vacuum::Physical ? Vacuum::SingleProduct
                                             : Vacuum::Physical);
  }));
  // the scoped contexts have ended
  CHECK(current_contexts_version() == version);

  // an unmodified copy of a context keeps the version
  CHECK(!changed_by([](Context&) {}));

  // a change of the default context changes it until the default is restored
  const auto bose_einstein = get_default_context(Statistics::BoseEinstein);
  // a copy given its own registry, which is a configuration of its own
  auto changed = bose_einstein;
  changed.set(IndexSpaceRegistry(*bose_einstein.index_space_registry()));
  set_default_context(changed, Statistics::BoseEinstein);
  CHECK(current_contexts_version() != version);
  set_default_context(bose_einstein, Statistics::BoseEinstein);
  CHECK(current_contexts_version() == version);
}

TEST_CASE("canonicalize_canonical", "[algorithms]") {
  using namespace sequant;

  // only Complete canonicalization marks its result
  auto complete_ctx = get_default_context();
  complete_ctx.set(
      CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  auto complete_resetter = set_scoped_default_context(complete_ctx);

  auto make_product = [] {
    return ex<Constant>(rational{1, 2}) *
           ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      Symmetry::Antisymm) *
           ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_2", L"i_1"},
                      Symmetry::Antisymm);
  };
  auto make_sum = [&] {
    return make_product() +
           ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_1"}) *
               ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}) +
           ex<Constant>(3);
  };
  const auto opts = CanonicalizeOptions::default_options();

  SECTION("canonicalize marks its result") {
    for (auto e : {make_product(), make_sum()}) {
      canonicalize(e);
      REQUIRE(e->is_canonical());
      REQUIRE(e->is_canonical(opts));
    }
    // also when the byproduct is absorbed into a new Product
    auto t = ex<Tensor>(L"g", bra{L"i_2", L"i_1"}, ket{L"a_1", L"a_2"},
                        Symmetry::Antisymm);
    canonicalize(t);
    REQUIRE(t.is<Product>());
    REQUIRE(t->is_canonical());
    REQUIRE(t->as<Product>().factor(0)->is_canonical());
  }

  SECTION("rapid canonicalization does not mark") {
    auto e = make_sum();
    canonicalize(e, opts.copy_and_set(CanonicalizationMethod::Rapid));
    REQUIRE(!e->is_canonical());
    e->rapid_canonicalize();
    REQUIRE(!e->is_canonical());
    // and invalidates the mark
    canonicalize(e);
    REQUIRE(e->is_canonical());
    e->rapid_canonicalize();
    REQUIRE(!e->is_canonical());
  }

  SECTION("mutation invalidates, canonicalize marks again") {
    auto e = make_product();
    canonicalize(e);
    REQUIRE(e->is_canonical());
    e->as<Product>().append(1, ex<Tensor>(L"f", bra{L"i_3"}, ket{L"i_4"}));
    REQUIRE(!e->is_canonical());
    canonicalize(e);
    REQUIRE(e->is_canonical());
    e->as<Product>().scale(2);
    REQUIRE(!e->is_canonical());
    canonicalize(e);
    REQUIRE(e->is_canonical());
    // in-place mutation of a factor leaves the memoized hash of the Product
    // stale, but must neither throw nor leave the mark valid
    e->hash_value();
    auto& t = e->as<Product>().factor(1)->as<Tensor>();
    t.transform_indices(
        container::map<Index, Index>{{Index{L"i_1"}, Index{L"i_5"}}});
    t.reset_tags();
    REQUIRE(!e->is_canonical());
    canonicalize(e);
    REQUIRE(e->is_canonical());
  }

  SECTION("in-place mutation of a summand's factor invalidates") {
    auto e = make_sum();
    canonicalize(e);
    e->hash_value();
    // the first summand that is a Product, wherever canonical order puts it
    auto product = ranges::find_if(e->as<Sum>().summands(), [](const auto& s) {
      return s->template is<Product>();
    });
    REQUIRE(product != ranges::end(e->as<Sum>().summands()));
    auto& t = (*product)->as<Product>().factor(0)->as<Tensor>();
    t.transform_indices(
        container::map<Index, Index>{{Index{L"i_1"}, Index{L"i_5"}}});
    t.reset_tags();
    REQUIRE(!e->is_canonical());
    REQUIRE_NOTHROW(canonicalize(e));
    REQUIRE(e->is_canonical());
  }

  SECTION("canonicalizing a canonical expression is a no-op") {
    for (auto e : {make_product(), make_sum()}) {
      REQUIRE(count_product_canonicalizations([&] { canonicalize(e); }) > 0);
      const auto* ptr = e.get();
      const auto latex = to_latex(e);
      const auto hash = e->hash_value();
      REQUIRE(count_product_canonicalizations([&] { canonicalize(e); }) == 0);
      REQUIRE(count_product_canonicalizations([&] { e->canonicalize(); }) == 0);
      REQUIRE(e.get() == ptr);
      REQUIRE(to_latex(e) == latex);
      REQUIRE(e->hash_value() == hash);
    }
    // simplify() canonicalizes, so a simplified expression is not
    // canonicalized again
    auto e = make_sum();
    simplify(e);
    REQUIRE(e->is_canonical());
    REQUIRE(count_product_canonicalizations([&] { simplify(e); }) == 0);
  }

  SECTION("Topological canonicalization alone does not mark") {
    auto topological_ctx = get_default_context();
    topological_ctx.set(
        CanonicalizeOptions{.method = CanonicalizationMethod::Topological});
    auto topological_resetter = set_scoped_default_context(topological_ctx);
    // a lone tensor and the same tensor scaled are spelled differently by
    // Topological canonicalization alone, so a summand it leaves behind must
    // still merge with a fresh copy
    ExprPtr z = deserialize(L"X{i_1;a_1} - X{i_6,i_3,i_5;a_1,a_2,a_4}:A-C-S");
    simplify(z);
    REQUIRE(!z->is_canonical());
    ExprPtr d = z - deserialize(serialize(z));
    REQUIRE(simplify(d) == ex<Constant>(0));
  }

  SECTION("canonicalizing a Sum leaves its canonical summands alone") {
    auto e = make_sum();
    canonicalize(e);
    const auto extra = ex<Tensor>(L"h", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
    e->as<Sum>().append(extra);
    REQUIRE(!e->is_canonical());
    // only the new summand is canonicalized, in the rapid and in the full pass
    REQUIRE(count_product_canonicalizations([&] { canonicalize(e); }) == 2);
    REQUIRE(e->is_canonical());
    // the result is that of canonicalizing from scratch
    auto reference = make_sum() + extra->clone();
    canonicalize(reference);
    REQUIRE(to_latex(e) == to_latex(reference));
  }

  SECTION("a change of global state invalidates") {
    auto e = make_sum();
    canonicalize(e);
    REQUIRE(e->is_canonical());
    {
      auto ctx = get_default_context();
      auto labels = ctx.cardinal_tensor_labels();
      labels.push_back(L"canonical_test");
      ctx.set_cardinal_tensor_labels(std::move(labels));
      auto _ = set_scoped_default_context(ctx);
      REQUIRE(!e->is_canonical());
      REQUIRE(count_product_canonicalizations([&] { canonicalize(e); }) > 0);
      REQUIRE(e->is_canonical());
    }
    REQUIRE(!e->is_canonical());
    canonicalize(e);
    REQUIRE(e->is_canonical());

    {
      auto _ = set_scoped_default_context(
          Context(get_default_context())
              .set_tensor_canonicalizer(
                  L"canonical_test",
                  std::make_shared<DefaultTensorCanonicalizer>()));
      REQUIRE(!e->is_canonical());
      canonicalize(e);
      REQUIRE(e->is_canonical());
    }
    REQUIRE(!e->is_canonical());

    // a setting canonicalization does not read does not invalidate
    canonicalize(e);
    {
      auto ctx = get_default_context();
      ctx.set(ctx.vacuum() == Vacuum::Physical ? Vacuum::SingleProduct
                                               : Vacuum::Physical);
      auto _ = set_scoped_default_context(ctx);
      REQUIRE(e->is_canonical());
    }
  }

  SECTION("different options invalidate") {
    auto e = make_sum();
    canonicalize(e);
    const auto other_opts =
        opts.copy_and_set(container::set<Index>{Index{L"i_1"}});
    REQUIRE(!e->is_canonical(other_opts));
    REQUIRE(count_product_canonicalizations(
                [&] { canonicalize(e, other_opts); }) > 0);
    REQUIRE(e->is_canonical(other_opts));
    REQUIRE(!e->is_canonical(opts));
  }

  SECTION("a clone of a canonical expression is canonical") {
    for (auto e : {make_product(), make_sum()}) {
      canonicalize(e);
      auto c = e->clone();
      REQUIRE(c->is_canonical());
      REQUIRE(count_product_canonicalizations([&] { canonicalize(c); }) == 0);
      REQUIRE(c == e);
    }
  }
}

TEST_CASE("canonicalize_tensor_with_repeated_index", "[algorithms]") {
  using namespace sequant;

  for (const auto method : {CanonicalizationMethod::Topological,
                            CanonicalizationMethod::Complete}) {
    auto ctx = get_default_context();
    ctx.set(CanonicalizeOptions{.method = method});
    auto resetter = set_scoped_default_context(ctx);

    // two spellings of one tensor network, differing in dummy names
    auto require_equal = [](const ExprPtr& x, const ExprPtr& y) {
      REQUIRE(simplify(x - y) == ex<Constant>(0));
    };
    // a slot index repeated in bra and ket
    require_equal(ex<Tensor>(L"h", bra{L"i_7"}, ket{L"i_7"}),
                  ex<Tensor>(L"h", bra{L"i_1"}, ket{L"i_1"}));
    // ... next to a named index
    require_equal(ex<Tensor>(L"h", bra{L"i_7", L"i_2"}, ket{L"i_7", L"a_1"}),
                  ex<Tensor>(L"h", bra{L"i_1", L"i_2"}, ket{L"i_1", L"a_1"}));
    // a repeated index carrying protoindices
    require_equal(ex<Tensor>(L"X", bra{Index(L"a_1", {L"i_1"})},
                             ket{Index(L"a_1", {L"i_1"})}),
                  ex<Tensor>(L"X", bra{Index(L"a_2", {L"i_1"})},
                             ket{Index(L"a_2", {L"i_1"})}));
    // an index repeated only as a protoindex: its lone and scaled spellings
    // merge
    {
      auto x = [] {
        return ex<Tensor>(L"X", bra{Index(L"a_1", {L"i_1"})},
                          ket{Index(L"a_2", {L"i_1"})});
      };
      REQUIRE(simplify(x() - ex<Constant>(2) * x()) ==
              simplify(ex<Constant>(-1) * x()));
    }
    // a lone and a scaled spelling merge
    REQUIRE(simplify(ex<Tensor>(L"h", bra{L"i_7"}, ket{L"i_7"}) -
                     ex<Constant>(2) *
                         ex<Tensor>(L"h", bra{L"i_1"}, ket{L"i_1"})) ==
            simplify(ex<Constant>(-1) *
                     ex<Tensor>(L"h", bra{L"i_1"}, ket{L"i_1"})));
  }
}
