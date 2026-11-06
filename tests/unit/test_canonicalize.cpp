#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
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

#include <algorithm>
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

  auto ibr = sequant::mbpt::make_legacy_spaces();
  mbpt::add_pao_basis(ibr, mbpt::Spin::any);
  auto ctx = get_default_context();
  ctx.set(ibr);
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("Tensors") {
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                           particle_symmetric);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_3", L"p_4"},
                           particle_symmetric);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           particle_symmetric);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                           particle_symmetric);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::Symm);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}:S"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::Antisymm);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("-g{p1,p2;p3,p4}:A"));
    }

    // aux indices
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_2"}, aux{L"p_3"},
                           particle_symmetric);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("B{p1;p2;p3}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, particle_symmetric);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("B{p1,p2;p4,p3;p5}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::Symm);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("B{p1,p2;p3,p4;p5}:S"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::Antisymm);
      canonicalize(op, deduce_named_indices);
      REQUIRE_THAT(op, SimplifiesTo("-B{p1,p2;p3,p4;p5}:A"));
    }
  }

  SECTION("Products") {
    // P.S. ref outputs produced with complete canonicalization
    {
      auto input =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_5"}, ket{L"a_1"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_5", L"a_2"},
                     particle_symmetric);
      canonicalize(input, deduce_named_indices);
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
      canonicalize(input, deduce_named_indices);
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
        const CanonicalizeOptions opts{
            .method = CanonicalizationMethod::Topological,
            .ignore_named_index_labels =
                static_cast<CanonicalizeOptions::IgnoreNamedIndexLabel>(
                    ignore_named_index_labels)};
        auto input1 =
            deserialize(L"1/2 t{a3,a1,a2;i4,i5,i2}:N-C-S g{i4,i5;i3,i1}:N-C-S");
        //      auto input1 = deserialize(L"1/2
        //      t{a1,a2,a3;i5,i2,i4}:N-C-S g{i4,i5;i3,i1}:N-C-S");
        auto input2 =
            deserialize(L"1/2 t{a1,a3,a2;i5,i4,i2}:N-C-S g{i5,i4;i1,i3}:N-C-S");
        canonicalize(input1, deduce_named_indices, opts);
        canonicalize(input2, deduce_named_indices, opts);
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
      canonicalize(input, deduce_named_indices);
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
      canonicalize(input1, deduce_named_indices);
      REQUIRE_THAT(input1,
                   SimplifiesTo("Ŝ{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "f⁺{i_1,i_3;a_1,a_2}:N-N-S t{i_2;a_3}"));
      auto input2 =
          ex<Tensor>(reserved::symm_label(), bra{L"a_1", L"a_2"},
                     ket{L"i_1", L"i_2"}, particle_symmetric) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, particle_symmetric) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, particle_symmetric) * f2 *
          ex<Variable>(L"w") * ex<Constant>(rational{1, 2});
      canonicalize(input2, deduce_named_indices);
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
      canonicalize(input, deduce_named_indices);
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
        if (!k.empty()) sr_reg->field(k, Field::Real);
      // Disable strict bra↔ket-symmetry policy: this expression has a_3 in
      // g.bra and t.bra under one term's orientation (a bra-bra contraction,
      // legitimate for Symm-braket g), which the default-context Conjugate
      // policy would reject.
      auto srcc_resetter = set_scoped_default_context(
          Context({.index_basis_registry_shared_ptr = sr_reg,
                   .vacuum = Vacuum::SingleProduct,
                   .spbasis = SPBasis::Spinfree})
              .set(AssertStrictBraKetSymmetry::No));
      auto input = deserialize(
          L"8 * Ŝ{i_1,i_2;a_1,a_2} * g{i_3,a_1;a_3,i_1}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S "
          L"+ 8 * Ŝ{i_1,i_2;a_1,a_2} * g{i_1,a_3;a_1,i_3}:N-S-S "
          L"* t{a_2,a_3;i_2,i_3}:N-N-S");
      simplify(input, deduce_named_indices);
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
      auto mdA = tnA.canonicalize_slots({}, named);
      auto mdB = tnB.canonicalize_slots({}, named);
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
      mbpt::add_pao_basis(sr_reg, mbpt::Spin::any);  // μ̃ (PAO)
      mbpt::add_df_spaces(sr_reg);                   // Κ  (DF aux)
      std::vector<std::wstring> keys;
      for (auto const& s : *sr_reg) keys.push_back(s.base_key());
      for (auto const& k : keys)
        if (!k.empty()) sr_reg->field(k, Field::Real);
      Context ctx = get_default_context();
      ctx.set(sr_reg);
      ctx.set(AssertStrictBraKetSymmetry::No);
      auto resetter = set_scoped_default_context(ctx);

      auto graph_cmp = [](std::wstring a, std::wstring b) {
        auto exA = deserialize(a), exB = deserialize(b);
        REQUIRE(exA);
        REQUIRE(exB);
        TensorNetworkV3 tnA(exA), tnB(exB);
        auto mdA = tnA.canonicalize_slots({}, deduce_named_indices);
        auto mdB = tnB.canonicalize_slots({}, deduce_named_indices);
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
        if (!k.empty()) sr_reg->field(k, Field::Real);
      auto srcc_resetter = set_scoped_default_context(
          Context({.index_basis_registry_shared_ptr = sr_reg,
                   .vacuum = Vacuum::SingleProduct,
                   .spbasis = SPBasis::Spinfree})
              .set(AssertStrictBraKetSymmetry::No));
      auto input = deserialize(tests::data::sf_r2_direct_real());
      REQUIRE(input);
      const std::size_t n_before =
          input->is<Sum>() ? input->size() : std::size_t{1};
      REQUIRE(n_before == 113);
      simplify(input, deduce_named_indices);
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
      simplify(input, deduce_named_indices);
      canonicalize(input, deduce_named_indices);
      REQUIRE_THAT(input, EquivalentTo("q2 + 2 q1"));
    }

    {
      auto input =
          ex<Variable>(L"q1") * ex<Variable>(L"q1") + ex<Variable>(L"q2");
      simplify(input, deduce_named_indices);
      canonicalize(input, deduce_named_indices);
      REQUIRE_THAT(input, EquivalentTo("q2 + q1 * q1"));
    }
  }

  SECTION("Powers in Products") {
    const auto f = deserialize(L"f{p_1;p_2}:A-C-S * ã{p_2;p_1}");
    const auto t1 = deserialize(L"t{a_1;i_1}:A-C-S * ã{i_1;a_1}");
    const auto pw = ex<Power>(ex<Variable>(L"x"), rational{1, 2});

    auto expr1 = f * t1 * pw;
    simplify(expr1, deduce_named_indices);
    REQUIRE_THAT(expr1, EquivalentTo(L"x^(1/2) * ã{p_2;p_1} * ã{i_1;a_1} * "
                                     L"t{a_1;i_1}:A-C-S * f{p_1;p_2}:A-C-S"));

    auto expr2 = ex<Constant>(rational{1, 2}) * f * t1 * pw;
    simplify(expr2, deduce_named_indices);
    REQUIRE_THAT(expr2, EquivalentTo(L"1/2 x^(1/2) * ã{p_2;p_1} * ã{i_1;a_1} * "
                                     L"t{a_1;i_1}:A-C-S * f{p_1;p_2}:A-C-S"));

    auto expr3 = f * t1 * ex<Power>(2, 3);
    simplify(expr3, deduce_named_indices);
    REQUIRE_THAT(expr3, EquivalentTo(L"8 * ã{p_2;p_1} * ã{i_1;a_1} * "
                                     L"t{a_1;i_1}:A-C-S * f{p_1;p_2}:A-C-S"));
  }

  SECTION("Sum of Powers") {
    const auto vx = ex<Variable>(L"x");

    // x^{1/2} + x^{1/2} = 2 * x^{1/2}
    auto pw1 = ex<Power>(vx, rational{1, 2});
    auto pw2 = ex<Power>(vx, rational{1, 2});
    auto sum_expr = pw1 + pw2;
    simplify(sum_expr, deduce_named_indices);
    REQUIRE(sum_expr->is<Product>());
    REQUIRE(sum_expr->as<Product>().scalar() == 2);

    // 2 x^{1/2} + 3 x^{1/2} = 5 x^{1/2}
    auto s1 = ex<Constant>(2) * ex<Power>(vx, rational{1, 2});
    auto s2 = ex<Constant>(3) * ex<Power>(vx, rational{1, 2});
    auto sum2 = s1 + s2;
    simplify(sum2, deduce_named_indices);
    REQUIRE(sum2->is<Product>());
    REQUIRE(sum2->as<Product>().scalar() == 5);
  }

  SECTION("Sum of Products") {
    // P.S. ref outputs produced with complete canonicalization
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
      simplify(input, deduce_named_indices);
      canonicalize(input, deduce_named_indices);
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
      canonicalize(input, deduce_named_indices);
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
      canonicalize(input, deduce_named_indices);
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
      canonicalize(input, deduce_named_indices);
      // the summands merge into one term, which canonicalization returns
      REQUIRE(input->is<Product>());
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

      canonicalize(input, deduce_named_indices);
      // the summands merge into one term, which canonicalization returns
      REQUIRE(input->is<Product>());
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
        canonicalize(input, deduce_named_indices);
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
        canonicalize(term1, deduce_named_indices);
        canonicalize(term2, deduce_named_indices);
        REQUIRE_THAT(term1,
                     EquivalentTo("-4 Ŝ{i_1,i_2,i_3;a_1,a_2,a_3} f{i_4;i_3} "
                                  "t{a_1,a_2,a_3;i_1,i_4,i_2}"));
        REQUIRE_THAT(term2,
                     EquivalentTo("-4 Ŝ{i_1,i_2,i_3;a_1,a_2,a_3} f{i_4;i_3} "
                                  "t{a_1,a_2,a_3;i_1,i_4,i_2}"));
        auto sum_of_terms = term1 + term2;
        simplify(sum_of_terms, deduce_named_indices);
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
        canonicalize(input, deduce_named_indices);
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

      canonicalize(input, deduce_named_indices);
      simplify(input, deduce_named_indices);
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
    return tn
        .canonicalize_slots(default_cardinal_tensor_labels(),
                            deduce_named_indices)
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

TEST_CASE("explicit_presentation", "[algorithms][canonicalize]") {
  using namespace sequant;

  auto Q = [](std::initializer_list<std::wstring_view> b,
              std::initializer_list<std::wstring_view> k) {
    return ex<Tensor>(L"Q", bra(b), ket(k), Symmetry::Antisymm);
  };
  auto with = [](CanonicalizationPresentation presentation) {
    return CanonicalizeOptions{
        .presentation = std::make_shared<const CanonicalizationPresentation>(
            std::move(presentation))};
  };

  SECTION("lone tensor uses the canonicalizer for its label") {
    // per-label canonicalizer given in the options ...
    auto q = Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"});
    CHECK(q->as<Tensor>().canonicalize(
              deduce_named_indices,
              with({.tensor_canonicalizers = {
                        {L"Q",
                         NullTensorCanonicalizer::instance()}}})) == nullptr);
    CHECK(*q == *Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"}));
    // the default canonicalizer sorts the bra and the ket of Q
    {
      const auto bp = q->as<Tensor>().canonicalize(deduce_named_indices);
      REQUIRE(bp);
      CHECK(bp->as<Constant>().value<int>() == -1);
      CHECK(*q == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
    }
    // a Q-labeled tensor alone and in a product are canonicalized alike, so
    // t and 2t combine
    {
      auto sum = Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"}) +
                 ex<Constant>(2) * Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"});
      canonicalize(sum, deduce_named_indices,
                   with({.tensor_canonicalizers = {
                             {L"Q", NullTensorCanonicalizer::instance()}}}));
      rapid_simplify(sum);
      CHECK(sum->is<Product>());
    }
  }

  SECTION("the canonicalizer for the empty label replaces the default") {
    auto q = Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"});
    CHECK(q->as<Tensor>().canonicalize(
              deduce_named_indices,
              with({.tensor_canonicalizers = {
                        {L"",
                         NullTensorCanonicalizer::instance()}}})) == nullptr);
    CHECK(*q == *Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"}));
  }

  SECTION("TN canonicalization uses the canonicalizer for the label") {
    auto canonicalized = [&Q](const CanonicalizeOptions& opts) {
      TensorNetwork tn(ExprPtrList{Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"})});
      const auto byproduct = tn.canonicalize(default_cardinal_tensor_labels(),
                                             deduce_named_indices, opts);
      const int phase = byproduct ? byproduct->as<Constant>().value<int>() : 1;
      return std::make_pair(std::dynamic_pointer_cast<Expr>(tn.tensors().at(0)),
                            phase);
    };
    {
      // a canonicalizer for label Q, which leaves Q as is and counts its uses
      struct CountingCanonicalizer : NullTensorCanonicalizer {
        std::shared_ptr<int> calls = std::make_shared<int>(0);
        ExprPtr apply(AbstractTensor&) const override {
          ++*calls;
          return {};
        }
      };
      const auto counting = std::make_shared<CountingCanonicalizer>();
      const auto [tensor, phase] =
          canonicalized(with({.tensor_canonicalizers = {{L"Q", counting}}}));
      CHECK(*counting->calls > 0);
      // the graph-based pass alone places the named indices in label order
      CHECK(*tensor == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
      CHECK(phase == -1);
    }
    const auto [tensor, phase] = canonicalized(CanonicalizeOptions{});
    CHECK(*tensor == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
    CHECK(phase == -1);
  }

  SECTION("index comparer") {
    const auto default_cmp = TensorCanonicalizer::default_index_comparer();
    const auto reversed = std::make_shared<const tensor_index_comparer_t>(
        [default_cmp](const Index& idx1, const Index& idx2) {
          return default_cmp(idx2, idx1);
        });
    auto q = Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"});
    q->as<Tensor>().canonicalize(deduce_named_indices,
                                 with({.index_comparer = reversed}));
    CHECK(*q == *Q({L"i_2", L"i_1"}, {L"a_2", L"a_1"}));
    // without presentation settings the default order is restored
    CHECK(q->as<Tensor>().canonicalize(deduce_named_indices) == nullptr);
    CHECK(*q == *Q({L"i_1", L"i_2"}, {L"a_1", L"a_2"}));
  }

  SECTION("cardinal tensor labels") {
    auto labels_of = [](const ExprPtr& product) {
      REQUIRE(product->is<Product>());
      std::vector<std::wstring> labels;
      for (const auto& factor : product->as<Product>().factors())
        labels.emplace_back(factor->as<Tensor>().label());
      return labels;
    };
    auto make = [] {
      return ex<Tensor>(L"Z", bra{L"i_1"}, ket{L"a_1"}) *
             ex<Tensor>(L"Y", bra{L"a_1"}, ket{L"i_1"});
    };
    const std::vector<std::wstring> y_first{L"Y", L"Z"};
    const std::vector<std::wstring> z_first{L"Z", L"Y"};

    auto internal = make();
    canonicalize(internal, deduce_named_indices, CanonicalizeOptions{});
    CHECK(labels_of(internal) == y_first);
    CHECK(internal->is_canonical(deduce_named_indices, CanonicalizeOptions{}));

    const auto z_opts =
        with({.cardinal_tensor_labels = container::vector<std::wstring>{L"Z"}});
    // a form with presentation settings is not recorded as canonical ...
    auto presented = make();
    canonicalize(presented, deduce_named_indices, z_opts);
    CHECK(labels_of(presented) == z_first);
    CHECK(!presented->is_canonical(deduce_named_indices, z_opts));
    CHECK(
        !presented->is_canonical(deduce_named_indices, CanonicalizeOptions{}));
    // ... and a canonical mark does not keep a presentation from being applied
    auto clone = internal->clone();
    canonicalize(clone, deduce_named_indices, z_opts);
    CHECK(labels_of(clone) == z_first);
  }

  SECTION("simplify accepts only the library's form") {
    auto expr = Q({L"i_2", L"i_1"}, {L"a_1", L"a_2"});
    CHECK_THROWS_AS(simplify(expr, deduce_named_indices,
                             with({.cardinal_tensor_labels =
                                       container::vector<std::wstring>{L"Q"}})),
                    Exception);
  }

  SECTION("options compare presentations by value") {
    const auto a =
        with({.cardinal_tensor_labels = container::vector<std::wstring>{L"Z"}});
    const auto b =
        with({.cardinal_tensor_labels = container::vector<std::wstring>{L"Z"}});
    const auto c =
        with({.cardinal_tensor_labels = container::vector<std::wstring>{L"Y"}});
    CHECK(a == b);
    CHECK(!(a == c));
    CHECK(!(a == CanonicalizeOptions{}));
    CHECK(
        a.copy_and_set(std::shared_ptr<const CanonicalizationPresentation>{}) ==
        CanonicalizeOptions{});
  }
}

TEST_CASE("canonicalization_zero_by_symmetry", "[algorithms][canonicalize]") {
  using namespace sequant;

  auto ibr = sequant::mbpt::make_legacy_spaces();
  auto ctx = get_default_context();
  ctx.set(ibr);
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
    canonicalize(zero, deduce_named_indices);
    REQUIRE(zero->is_zero());
    // a zero product keeps no factors, so zeros of different inputs agree
    REQUIRE(zero->as<Product>().factors().empty());
    {
      ExprPtr other_zero =
          ex<Tensor>(L"t", bra{L"a_2", L"a_1"}, ket{L"i_2", L"i_1"},
                     Symmetry::Symm) *
          ex<FNOperator>(cre({L"a_2", L"a_1"}), ann({L"i_1", L"i_2"}));
      canonicalize(other_zero, deduce_named_indices);
      REQUIRE(*other_zero == *zero);
      REQUIRE(other_zero->hash_value() == zero->hash_value());
    }
    simplify(zero, deduce_named_indices);
    REQUIRE(is_zero(zero));

    auto nonzero = make(Symmetry::Antisymm);
    simplify(nonzero, deduce_named_indices);
    REQUIRE(!nonzero->is_zero());
  }

  SECTION("symmetric tensor contracted with antisymmetric tensor") {
    // dummy indices only
    {
      auto expr = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"h", bra{L"p_3", L"p_4"}, ket{L"p_1", L"p_2"},
                             Symmetry::Symm);
      simplify(expr, deduce_named_indices);
      REQUIRE(is_zero(expr));
    }
    // named indices are not permuted, the contracted pair is
    {
      auto make = [] {
        return ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                          Symmetry::Antisymm) *
               ex<Tensor>(L"w", bra{L"p_5", L"p_6"}, ket{L"p_1", L"p_2"},
                          Symmetry::Symm);
      };
      auto expr = make();
      simplify(expr, deduce_named_indices);
      REQUIRE(is_zero(expr));

      // in a sum named index labels are meaningful
      auto sum = make() + ex<Tensor>(L"f", bra{L"p_3", L"p_4"},
                                     ket{L"p_5", L"p_6"}, Symmetry::Antisymm);
      simplify(sum, deduce_named_indices);
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
      simplify(expr, deduce_named_indices);
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
      simplify(expr, deduce_named_indices);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo("g{i3,i4;a3,a4}:A t{a3,a4;i1,i2}:A "
                                      "t{a1,a2;i3,i4}:A"));
    }
    // a lone antisymmetric tensor: swapping its named indices is not a
    // symmetry of the term
    {
      auto expr = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                             Symmetry::Antisymm);
      simplify(expr, deduce_named_indices);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo("g{p1,p2;p3,p4}:A"));
    }
    // automorphisms permuting antisymmetric bundles of unequal or odd size
    // in pairs are even
    for (std::wstring spec :
         {L"t{a1,a2;i1}:A u{;a1,a2}:A", L"X{a1,a2,a3;i1}:A Y{;a1,a2,a3}:A"}) {
      auto expr = deserialize(spec);
      simplify(expr, deduce_named_indices);
      REQUIRE(!expr->is_zero());
      REQUIRE_THAT(expr, EquivalentTo(spec));
    }
    {
      auto sum = deserialize(L"t{a1,a2;i1}:A u{;a1,a2}:A + w{;i1}:A");
      simplify(sum, deduce_named_indices);
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
    simplify(zero, deduce_named_indices);
    REQUIRE(is_zero(zero));

    // w antisymmetric too: the exchange is even
    auto nonzero = make(Symmetry::Antisymm);
    simplify(nonzero, deduce_named_indices);
    REQUIRE(!nonzero->is_zero());
  }

  // only topological canonicalization looks for automorphisms
  SECTION("rapid canonicalization does not detect a zero") {
    auto expr = ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                           Symmetry::Symm) *
                ex<FNOperator>(cre({L"a_1", L"a_2"}), ann({L"i_1", L"i_2"}));
    expr->canonicalize(
        deduce_named_indices,
        CanonicalizeOptions{}.copy_and_set(CanonicalizationMethod::Rapid));
    REQUIRE(!expr->is_zero());
    canonicalize(expr, deduce_named_indices);
    REQUIRE(expr->is_zero());
  }

  // odd-size bundles: antisymmetric in a1,a2,a3 against symmetric
  SECTION("odd-size antisymmetric bundle contracted with symmetric one") {
    auto expr = deserialize(L"X{a1,a2,a3;i1}:A Y{;a1,a2,a3}:S");
    simplify(expr, deduce_named_indices);
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
      canonicalize(e, deduce_named_indices, opts);
      return e;
    };
    const auto cx = canonicalized(x);
    CAPTURE(x, y, serialize(cx));
    REQUIRE(cx == canonicalized(y));
    auto again = cx->clone();
    canonicalize(again, deduce_named_indices, opts);
    REQUIRE(again == cx);
  }
}

TEST_CASE("canonicalize_options_equality", "[algorithms]") {
  using namespace sequant;
  const CanonicalizeOptions opts{.method = CanonicalizationMethod::Complete};
  REQUIRE(opts == CanonicalizeOptions{opts});
  // every member takes part, not only the method
  REQUIRE(!(opts == opts.copy_and_set(CanonicalizationMethod::Rapid)));
  REQUIRE(!(opts ==
            opts.copy_and_set(CanonicalizeOptions::IgnoreNamedIndexLabel::No)));
}

TEST_CASE("named_index_spec", "[algorithms]") {
  using namespace sequant;
  const NamedIndexSpec deduced = deduce_named_indices;
  const auto none = sequant::tests::explicit_named({});
  const auto i1 = sequant::tests::explicit_named({Index{L"i_1"}});
  REQUIRE(deduced.deduce());
  REQUIRE(deduced.indices_ptr() == nullptr);
  REQUIRE_THROWS_AS(deduced.indices(), Exception);
  REQUIRE(!none.deduce());
  REQUIRE(none.indices().empty());
  REQUIRE(i1.indices() == container::set<Index>{Index{L"i_1"}});
  REQUIRE(!(deduced == none));
  REQUIRE(!(none == i1));
  // the explicit named indices of a result are its indices
  const auto result = deserialize<ResultExpr>(L"R{a_1;i_1} = t{a_1;i_1}");
  REQUIRE(named_indices(result) ==
          sequant::tests::explicit_named({Index{L"a_1"}, Index{L"i_1"}}));
}

// the hooks of a lone tensor with a repeated index (a tensor network) default
// to the context's options; with no index named every index is a dummy, so
// the once-occurring a_5 is relabeled, which the deduced named indices (a_5
// among them) would not do
TEST_CASE("canonicalize_hook_defaults", "[algorithms]") {
  using namespace sequant;
  auto make_tensor = [] {
    return ex<Tensor>(L"t", bra{L"i_3"}, ket{L"i_3", L"a_5"});
  };
  const NamedIndexSpec all_dummy =
      sequant::tests::explicit_named(container::set<Index>{});

  SECTION("Tensor::canonicalize") {
    auto expected = make_tensor();
    expected->as<Tensor>().canonicalize(all_dummy);
    auto deduced = make_tensor();
    deduced->as<Tensor>().canonicalize(deduce_named_indices);
    REQUIRE(*expected != *deduced);
  }

  SECTION("Expr::rapid_canonicalize forwards its named indices") {
    auto expected = make_tensor();
    expected->canonicalize(all_dummy, CanonicalizeOptions{}.copy_and_set(
                                          CanonicalizationMethod::Rapid));
    auto deduced = make_tensor();
    deduced->canonicalize(
        deduce_named_indices,
        CanonicalizeOptions{}.copy_and_set(CanonicalizationMethod::Rapid));
    REQUIRE(*expected != *deduced);

    auto t = make_tensor();
    t->rapid_canonicalize(all_dummy);
    REQUIRE(*t == *expected);
  }
}

// a named index that the network does not contain (e.g. one named for a
// whole sum) does not affect
// the canonicalization of the network; the lexicographic pass relabels the
// anonymous indices in the order of the tensors they connect
TEST_CASE("canonicalize_named_index_absent", "[algorithms]") {
  using namespace sequant;
  const auto input = L"t{p_1;m_2} u{m_2;p_3} v{p_4;m_1} w{m_1;p_5}";
  const auto named = container::set<Index>{Index{L"p_1"}, Index{L"p_3"},
                                           Index{L"p_4"}, Index{L"p_5"}};
  auto canonicalized = [&](container::set<Index> named_indices) {
    auto e = deserialize(input);
    canonicalize(e, sequant::tests::explicit_named(std::move(named_indices)),
                 {.method = CanonicalizationMethod::Lexicographic});
    return e;
  };
  const auto expected =
      deserialize(L"t{p_1;m_1} u{m_1;p_3} v{p_4;m_2} w{m_2;p_5}");
  REQUIRE(canonicalized(named) == expected);
  auto with_absent = named;
  with_absent.emplace(L"p_2");
  REQUIRE(canonicalized(with_absent) == expected);
}

TEST_CASE("current_contexts_version", "[algorithms]") {
  using namespace sequant;
  const auto version = current_contexts_version();

  // any change to the configuration canonicalization depends on (the index
  // basis registry and the single-particle basis) of the context in effect for
  // any Statistics changes the version, for as long as it is in effect
  auto changed_by = [&](auto&& modify, Statistics s = Statistics::Arbitrary) {
    Context ctx(get_default_context(s));
    modify(ctx);
    auto resetter = set_scoped_default_context({{s, ctx}});
    return current_contexts_version() != version;
  };
  CHECK(changed_by([](Context& ctx) {
    ctx.set(ctx.spbasis() == SPBasis::Spinor ? SPBasis::Spinfree
                                             : SPBasis::Spinor);
  }));
  for (auto s : {Statistics::FermiDirac, Statistics::BoseEinstein,
                 Statistics::Arbitrary}) {
    // a copy of the registry is a configuration of its own
    CHECK(changed_by(
        [](Context& ctx) {
          ctx.set(IndexBasisRegistry(*ctx.index_basis_registry()));
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
  changed.set(IndexBasisRegistry(*bose_einstein.index_basis_registry()));
  set_default_context(changed, Statistics::BoseEinstein);
  CHECK(current_contexts_version() != version);
  set_default_context(bose_einstein, Statistics::BoseEinstein);
  CHECK(current_contexts_version() == version);
}

TEST_CASE("canonicalize_canonical", "[algorithms]") {
  using namespace sequant;

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
  const auto opts = CanonicalizeOptions{};

  SECTION("canonicalize marks its result") {
    for (auto e : {make_product(), make_sum()}) {
      canonicalize(e, deduce_named_indices);
      REQUIRE(e->is_canonical(deduce_named_indices));
      REQUIRE(e->is_canonical(deduce_named_indices, opts));
    }
    // also when the byproduct is absorbed into a new Product
    auto t = ex<Tensor>(L"g", bra{L"i_2", L"i_1"}, ket{L"a_1", L"a_2"},
                        Symmetry::Antisymm);
    canonicalize(t, deduce_named_indices);
    REQUIRE(t.is<Product>());
    REQUIRE(t->is_canonical(deduce_named_indices));
    REQUIRE(t->as<Product>().factor(0)->is_canonical(deduce_named_indices));
  }

  SECTION("rapid canonicalization does not mark") {
    auto e = make_sum();
    {
      const auto rapid = opts.copy_and_set(CanonicalizationMethod::Rapid);
      canonicalize(e, deduce_named_indices, rapid);
      REQUIRE(!e->is_canonical(deduce_named_indices, rapid));
      REQUIRE(!e->is_canonical(deduce_named_indices));
    }
    e->rapid_canonicalize(deduce_named_indices);
    REQUIRE(!e->is_canonical(deduce_named_indices));
    // and invalidates the mark
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
    e->rapid_canonicalize(deduce_named_indices);
    REQUIRE(!e->is_canonical(deduce_named_indices));
  }

  SECTION("mutation invalidates, canonicalize marks again") {
    auto e = make_product();
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
    e->as<Product>().append(1, ex<Tensor>(L"f", bra{L"i_3"}, ket{L"i_4"}));
    REQUIRE(!e->is_canonical(deduce_named_indices));
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
    e->as<Product>().scale(2);
    REQUIRE(!e->is_canonical(deduce_named_indices));
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
    // in-place mutation of a factor leaves the memoized hash of the Product
    // stale, but must neither throw nor leave the mark valid
    e->hash_value();
    auto& t = e->as<Product>().factor(1)->as<Tensor>();
    t.transform_indices(
        container::map<Index, Index>{{Index{L"i_1"}, Index{L"i_5"}}});
    t.reset_tags();
    REQUIRE(!e->is_canonical(deduce_named_indices));
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
  }

  SECTION("in-place mutation of a summand's factor invalidates") {
    auto e = make_sum();
    canonicalize(e, deduce_named_indices);
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
    REQUIRE(!e->is_canonical(deduce_named_indices));
    REQUIRE_NOTHROW(canonicalize(e, deduce_named_indices));
    REQUIRE(e->is_canonical(deduce_named_indices));
  }

  SECTION("canonicalizing a canonical expression is a no-op") {
    for (auto e : {make_product(), make_sum()}) {
      REQUIRE(count_product_canonicalizations(
                  [&] { canonicalize(e, deduce_named_indices); }) > 0);
      const auto* ptr = e.get();
      const auto latex = to_latex(e);
      const auto hash = e->hash_value();
      REQUIRE(count_product_canonicalizations(
                  [&] { canonicalize(e, deduce_named_indices); }) == 0);
      REQUIRE(count_product_canonicalizations(
                  [&] { e->canonicalize(deduce_named_indices); }) == 0);
      REQUIRE(e.get() == ptr);
      REQUIRE(to_latex(e) == latex);
      REQUIRE(e->hash_value() == hash);
    }
    // simplify() canonicalizes, so a simplified expression is not
    // canonicalized again
    auto e = make_sum();
    simplify(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
    REQUIRE(count_product_canonicalizations(
                [&] { simplify(e, deduce_named_indices); }) == 0);
  }

  SECTION("Topological canonicalization alone does not mark") {
    const CanonicalizeOptions topological{
        .method = CanonicalizationMethod::Topological};
    // a lone tensor and the same tensor scaled are spelled differently by
    // Topological canonicalization alone, so a summand it leaves behind must
    // still merge with a fresh copy
    ExprPtr z = deserialize(L"X{i_1;a_1} - X{i_6,i_3,i_5;a_1,a_2,a_4}:A-C-S");
    simplify(z, deduce_named_indices, topological);
    REQUIRE(!z->is_canonical(deduce_named_indices, topological));
    REQUIRE(!z->is_canonical(deduce_named_indices));
    ExprPtr d = z - deserialize(serialize(z));
    REQUIRE(simplify(d, deduce_named_indices, topological) == ex<Constant>(0));
  }

  SECTION("canonicalizing a Sum leaves its canonical summands alone") {
    auto e = make_sum();
    canonicalize(e, deduce_named_indices);
    const auto extra = ex<Tensor>(L"h", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
    e->as<Sum>().append(extra);
    REQUIRE(!e->is_canonical(deduce_named_indices));
    // only the new summand is canonicalized, in the rapid and in the full pass
    REQUIRE(count_product_canonicalizations(
                [&] { canonicalize(e, deduce_named_indices); }) == 2);
    REQUIRE(e->is_canonical(deduce_named_indices));
    // the result is that of canonicalizing from scratch
    auto reference = make_sum() + extra->clone();
    canonicalize(reference, deduce_named_indices);
    REQUIRE(to_latex(e) == to_latex(reference));
  }

  SECTION("a change of global state invalidates") {
    auto e = make_sum();
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));
    // a copy of the registry is a configuration of its own
    {
      auto ctx = get_default_context();
      ctx.set(IndexBasisRegistry(*ctx.index_basis_registry()));
      auto _ = set_scoped_default_context(ctx);
      REQUIRE(!e->is_canonical(deduce_named_indices));
      REQUIRE(count_product_canonicalizations(
                  [&] { canonicalize(e, deduce_named_indices); }) > 0);
      REQUIRE(e->is_canonical(deduce_named_indices));
    }
    REQUIRE(!e->is_canonical(deduce_named_indices));
    canonicalize(e, deduce_named_indices);
    REQUIRE(e->is_canonical(deduce_named_indices));

    {
      auto ctx = get_default_context();
      ctx.set(ctx.spbasis() == SPBasis::Spinor ? SPBasis::Spinfree
                                               : SPBasis::Spinor);
      auto _ = set_scoped_default_context(ctx);
      REQUIRE(!e->is_canonical(deduce_named_indices));
      canonicalize(e, deduce_named_indices);
      REQUIRE(e->is_canonical(deduce_named_indices));
    }
    REQUIRE(!e->is_canonical(deduce_named_indices));

    // a setting canonicalization does not read does not invalidate
    canonicalize(e, deduce_named_indices);
    {
      auto ctx = get_default_context();
      ctx.set(ctx.vacuum() == Vacuum::Physical ? Vacuum::SingleProduct
                                               : Vacuum::Physical);
      auto _ = set_scoped_default_context(ctx);
      REQUIRE(e->is_canonical(deduce_named_indices));
    }
  }

  SECTION("different options invalidate") {
    auto e = make_sum();
    canonicalize(e, deduce_named_indices);
    const NamedIndexSpec other =
        sequant::tests::explicit_named(container::set<Index>{Index{L"i_1"}});
    REQUIRE(!e->is_canonical(other));
    REQUIRE(count_product_canonicalizations([&] { canonicalize(e, other); }) >
            0);
    REQUIRE(e->is_canonical(other));
    REQUIRE(!e->is_canonical(deduce_named_indices));
  }

  SECTION("a clone of a canonical expression is canonical") {
    for (auto e : {make_product(), make_sum()}) {
      canonicalize(e, deduce_named_indices);
      auto c = e->clone();
      REQUIRE(c->is_canonical(deduce_named_indices));
      REQUIRE(count_product_canonicalizations(
                  [&] { canonicalize(c, deduce_named_indices); }) == 0);
      REQUIRE(c == e);
    }
  }
}

TEST_CASE("result_expr_canonicalization", "[algorithms][canonicalize]") {
  using namespace sequant;
  auto count = [](const ResultExpr& r) {
    return r.expression()->is<Sum>() ? r.expression()->size() : 1;
  };

  // the inner sum of S * (sum)
  auto inner = [](const ResultExpr& r) -> ExprPtr {
    REQUIRE(r.expression()->is<Product>());
    for (const auto& factor : r.expression()->as<Product>().factors())
      if (!factor->is<Tensor>() ||
          (factor->as<Tensor>().label() != reserved::symm_label() &&
           factor->as<Tensor>().label() != reserved::antisymm_label()))
        return factor;
    return nullptr;
  };
  auto size = [](const ExprPtr& e) { return e->is<Sum>() ? e->size() : 1; };

  SECTION("terms related by the result's symmetry are combined") {
    // under the symmetrizer of the result, the second term is the first with
    // the columns of the result swapped
    auto r = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-N-S = Ŝ{i_1,i_2;a_1,a_2}:N-N-S * "
        L"(X{a_1,a_2;i_1,i_2}:N-N-N + X{a_2,a_1;i_2,i_1}:N-N-N)");
    canonicalize(r);
    CHECK(size(inner(r)) == 1);
    // ... and so is a term that the symmetrizer multiplies on its own
    auto lone = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-N-S = Ŝ{i_1,i_2;a_1,a_2}:N-N-S * "
        L"X{a_1,a_2;i_1,i_2}:N-N-N + Ŝ{i_1,i_2;a_1,a_2}:N-N-S * "
        L"X{a_2,a_1;i_2,i_1}:N-N-N");
    simplify(lone);
    CHECK(count(lone) == 1);
    // ... but not without that symmetry
    auto nonsymm = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-N-N = Ŝ{i_1,i_2;a_1,a_2}:N-N-S * "
        L"(X{a_1,a_2;i_1,i_2}:N-N-N + X{a_2,a_1;i_2,i_1}:N-N-N)");
    canonicalize(nonsymm);
    CHECK(size(inner(nonsymm)) == 2);
    // antisymmetry absorbs the phase
    auto antisymm = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:A-N-S = Â{i_1,i_2;a_1,a_2}:A-N-S * "
        L"(X{a_1,a_2;i_1,i_2}:N-N-N - X{a_2,a_1;i_1,i_2}:N-N-N)");
    canonicalize(antisymm);
    const auto antisymm_inner = inner(antisymm);
    CHECK(size(antisymm_inner) == 1);
    const auto representative = antisymm_inner->is<Sum>()
                                    ? antisymm_inner->as<Sum>().summand(0)
                                    : antisymm_inner;
    REQUIRE(representative->is<Product>());
    // the representative's sign depends on which permutation it is
    CHECK(abs(representative->as<Product>().scalar()) == 2);
    // the phase of a term that includes the antisymmetrizer comes from
    // permuting the antisymmetrizer: A X - A X' = A X + A X = 2 A X , X' being
    // X with the virtual indices swapped
    auto antisymm_lone = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:A-N-S = Â{i_1,i_2;a_1,a_2}:A-N-S * "
        L"X{a_1,a_2;i_1,i_2}:N-N-N - Â{i_1,i_2;a_1,a_2}:A-N-S * "
        L"X{a_2,a_1;i_1,i_2}:N-N-N");
    simplify(antisymm_lone);
    REQUIRE(antisymm_lone.expression()->is<Product>());
    CHECK(abs(antisymm_lone.expression()->as<Product>().scalar()) == 2);
  }

  SECTION("terms are not combined without a symmetrizer") {
    // that the result has the symmetry does not make its terms have it: the
    // sum of the two is symmetric, either one alone is not
    auto r = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-N-S = X{a_1,a_2;i_1,i_2}:N-N-N + "
        L"X{a_2,a_1;i_2,i_1}:N-N-N");
    simplify(r);
    CHECK(count(r) == 2);
  }

  SECTION("only permutations within the block of the result") {
    // swapping the columns would exchange a virtual and an occupied index,
    // i.e. relate different blocks of R
    auto r = deserialize<ResultExpr>(
        L"R{a_1,i_3;i_1,i_2}:N-N-S = Ŝ{i_1,i_2;a_1,i_3}:N-N-S * "
        L"(X{a_1,i_3;i_1,i_2}:N-N-N + X{i_3,a_1;i_2,i_1}:N-N-N)");
    canonicalize(r);
    CHECK(size(inner(r)) == 2);
  }

  SECTION("external indices that occur more than once keep their labels") {
    // an element-wise (Hadamard) product: the externals occur twice, which
    // deduction would take for a contraction
    const auto input = L"R{a_1;i_1} = X{a_1;i_1} * Y{a_1;i_1}";
    auto r = deserialize<ResultExpr>(input);
    REQUIRE_NOTHROW(simplify(r));
    CHECK(r.expression() == deserialize(L"X{a_1;i_1} * Y{a_1;i_1}"));
    auto deduced = deserialize<ResultExpr>(input).expression();
    CHECK_THROWS_AS(simplify(deduced, deduce_named_indices), Exception);
  }
}

TEST_CASE("canonicalize_tensor_with_repeated_index", "[algorithms]") {
  using namespace sequant;

  for (const auto method : {CanonicalizationMethod::Topological,
                            CanonicalizationMethod::Complete}) {
    const CanonicalizeOptions opts{.method = method};

    // two spellings of one tensor network, differing in dummy names
    auto require_equal = [&opts](const ExprPtr& x, const ExprPtr& y) {
      REQUIRE(simplify(x - y, deduce_named_indices, opts) == ex<Constant>(0));
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
      REQUIRE(
          simplify(x() - ex<Constant>(2) * x(), deduce_named_indices, opts) ==
          simplify(ex<Constant>(-1) * x(), deduce_named_indices, opts));
    }
    // a lone and a scaled spelling merge
    REQUIRE(
        simplify(
            ex<Tensor>(L"h", bra{L"i_7"}, ket{L"i_7"}) -
                ex<Constant>(2) * ex<Tensor>(L"h", bra{L"i_1"}, ket{L"i_1"}),
            deduce_named_indices, opts) ==
        simplify(ex<Constant>(-1) * ex<Tensor>(L"h", bra{L"i_1"}, ket{L"i_1"}),
                 deduce_named_indices, opts));
  }
}

TEST_CASE("canonicalize_lone_term", "[algorithms]") {
  using namespace sequant;

  for (const auto method : {CanonicalizationMethod::Topological,
                            CanonicalizationMethod::Complete}) {
    const CanonicalizeOptions opts{.method = method};

    auto x = [](std::initializer_list<const wchar_t*> bra_idxs) {
      return ex<Tensor>(L"X", bra(bra_idxs), ket{L"a_1", L"a_2", L"a_4"},
                        Symmetry::Antisymm, BraKetSymmetry::Conjugate,
                        ColumnSymmetry::Symm);
    };
    // the summands of a Sum that merge into one term, and a Product of one
    // factor with a unit scalar, are canonicalized as that term
    for (auto input :
         {x({L"i_6", L"i_3", L"i_5"}) +
              ex<Constant>(2) * x({L"i_3", L"i_6", L"i_5"}),
          ex<Product>(1, ExprPtrList{x({L"i_6", L"i_3", L"i_5"})}),
          ex<Product>(-1, ExprPtrList{x({L"i_6", L"i_3", L"i_5"})})}) {
      simplify(input, deduce_named_indices, opts);
      REQUIRE(!input->is<Sum>());
      const auto once = serialize(input);
      simplify(input, deduce_named_indices, opts);
      REQUIRE(serialize(input) == once);
    }
  }
}

TEST_CASE("tensor_power_canonicalization", "[algorithms]") {
  using namespace sequant;
  auto ibr = mbpt::make_legacy_spaces();
  auto ctx = get_default_context();
  ctx.set(ibr);
  auto ctx_resetter = set_scoped_default_context(ctx);
  const TensorSymmetries dsym{.perm = Symmetry::Symm,
                              .braket = BraKetSymmetry::Symm,
                              .column = ColumnSymmetry::Symm};
  auto D = ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, dsym);
  auto Dp = ex<Tensor>(L"D", bra{L"i_2", L"i_1"}, ket{L"a_1", L"a_2"}, dsym);
  auto M = ex<Tensor>(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      TensorSymmetries{.perm = Symmetry::Antisymm});
  auto Mp = ex<Tensor>(L"M", bra{L"i_2", L"i_1"}, ket{L"a_1", L"a_2"},
                       TensorSymmetries{.perm = Symmetry::Antisymm});
  const NamedIndexSpec named =
      sequant::tests::explicit_named(container::set<Index>{
          Index{L"i_1"}, Index{L"i_2"}, Index{L"a_1"}, Index{L"a_2"}});

  // every product below has its externals in two slots of the same kind, so
  // without named indices they are deduced as dummies and three covariant
  // slots trip the strict bra-ket rule in TensorNetworkV3; canonicalize with
  // the externals named, as a ResultExpr head would
  auto canon = [&](ExprPtr e) {
    e->canonicalize(named);
    return e;
  };
  SECTION("Hadamard product with named externals canonicalizes") {
    auto e1 = canon(M * pow(D, rational{-1}));
    // i_1<->i_2 in both factors: M antisymm, D symm => e2 == -e1
    auto e2 = canon(Mp * pow(Dp, rational{-1}));
    REQUIRE(e1->is<Product>());
    REQUIRE(e2->is<Product>());
    // Product absorbs the phase into its scalar
    CHECK(e2->as<Product>().scalar() == -e1->as<Product>().scalar());
    CHECK(ranges::equal(e1->as<Product>().factors(),
                        e2->as<Product>().factors(),
                        [](auto const& x, auto const& y) { return *x == *y; }));
  }
  SECTION("power and base of the same label order uniquely") {
    auto a = canon(pow(D, rational{-1}) * pow(D, rational{2}) * M);
    auto b = canon(pow(D, rational{2}) * pow(D, rational{-1}) * M);
    CHECK(*a == *b);
    auto c = canon(D * pow(D, rational{-1}) * M);
    auto d = canon(pow(D, rational{-1}) * D * M);
    CHECK(*c == *d);
    CHECK(
        (D->as<AbstractTensor>() <
         pow(D, rational{-1})->as<AbstractTensor>()) !=
        (pow(D, rational{-1})->as<AbstractTensor>() < D->as<AbstractTensor>()));
  }
  SECTION("a power is not identified with its base by the network") {
    auto e = canon(D * pow(D, rational{-1}));
    REQUIRE(e->is<Product>());
    CHECK(e->as<Product>().factors().size() == 2);
    CHECK(canon(D * D)->hash_value() != e->hash_value());
  }
  SECTION("a lone power canonicalizes like a tensor") {
    auto p = pow(Dp, rational{-1});
    auto ref = pow(D, rational{-1});
    CHECK(!p->canonicalize(deduce_named_indices));
    CHECK(!ref->canonicalize(deduce_named_indices));
    CHECK(*p == *ref);

    auto s = pow(Dp, rational{-1}) - pow(D, rational{-1});
    simplify(s, deduce_named_indices);
    CHECK(s->is_zero());

    auto t = ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                        Symmetry::Antisymm);
    auto tp = ex<Tensor>(L"t", bra{L"i_2", L"i_1"}, ket{L"a_1", L"a_2"},
                         Symmetry::Antisymm);
    auto canon_with_phase = [](ExprPtr e) {
      auto bp = e->canonicalize(deduce_named_indices);
      return bp ? bp * e : e;
    };
    // even exponent: the phase squares to +1
    auto even = canon_with_phase(pow(tp, rational{2}));
    CHECK(*even == *canon_with_phase(pow(t, rational{2})));
    // odd exponent: the phase survives
    auto odd = canon_with_phase(pow(tp, rational{3}));
    CHECK(*odd == *(ex<Constant>(-1) * canon_with_phase(pow(t, rational{3}))));

    auto s2 = pow(tp, rational{3}) + pow(t, rational{3});
    simplify(s2, deduce_named_indices);
    CHECK(s2->is_zero());
  }
  SECTION("lexicographic canonicalization tells a power from its base") {
    // no cardinal labels, so that the order of the D factors is decided by
    // the tensor comparison
    auto opts =
        CanonicalizeOptions{}
            .copy_and_set(CanonicalizationMethod::Lexicographic)
            .copy_and_set(std::make_shared<const CanonicalizationPresentation>(
                CanonicalizationPresentation{
                    .cardinal_tensor_labels =
                        container::vector<std::wstring>{}}));
    auto X =
        pow(ex<Tensor>(L"D", bra{L"i_3"}, ket{L"a_3"}, dsym), rational{-1});
    auto Y = ex<Tensor>(L"D", bra{L"i_4"}, ket{L"a_4"}, dsym);
    auto Z = ex<Tensor>(L"t", bra{L"a_3", L"a_4"}, ket{L"i_3", L"i_4"});
    auto canon_lex = [&](ExprPtr e) {
      auto bp = e->canonicalize(deduce_named_indices, opts);
      return bp ? bp * e : e;
    };
    auto A = canon_lex(X->clone() * Y->clone() * Z->clone());
    auto B = canon_lex(Y->clone() * X->clone() * Z->clone());
    CHECK(*A == *B);
    // rapid canonicalization of a Sum canonicalizes its summands
    // lexicographically and merges them
    auto diff = X->clone() * Y->clone() * Z->clone() -
                Y->clone() * X->clone() * Z->clone();
    diff->rapid_canonicalize(deduce_named_indices);
    REQUIRE(diff->is<Sum>());
    CHECK(diff->as<Sum>().summands().empty());

    // two powers of one base differ only in the exponent
    auto Y2 =
        pow(ex<Tensor>(L"D", bra{L"i_4"}, ket{L"a_4"}, dsym), rational{2});
    auto A2 = canon_lex(X->clone() * Y2->clone() * Z->clone());
    auto B2 = canon_lex(Y2->clone() * X->clone() * Z->clone());
    CHECK(*A2 == *B2);
    auto diff2 = X->clone() * Y2->clone() * Z->clone() -
                 Y2->clone() * X->clone() * Z->clone();
    diff2->rapid_canonicalize(deduce_named_indices);
    REQUIRE(diff2->is<Sum>());
    CHECK(diff2->as<Sum>().summands().empty());
  }
  SECTION("relabelled duplicates merge in a sum") {
    // == M D^-1 - M D^-1 == 0
    auto s = M * pow(D, rational{-1}) + Mp * pow(Dp, rational{-1});
    simplify(s, deduce_named_indices);
    CHECK(s->is_zero());
  }
}

TEST_CASE("nested products", "[algorithms][canonicalize]") {
  using namespace sequant;

  // g (t (f t)): a product parenthesized for evaluation
  auto make = [](const wchar_t* a, const wchar_t* b, const wchar_t* i) {
    auto inner = ex<Product>(ExprPtrList{ex<Tensor>(L"f", bra{i}, ket{b}),
                                         ex<Tensor>(L"t", bra{b}, ket{L"i_1"})},
                             Product::Flatten::No);
    auto middle =
        ex<Product>(ExprPtrList{ex<Tensor>(L"t", bra{a}, ket{i}), inner},
                    Product::Flatten::No);
    return ex<Product>(
        ExprPtrList{ex<Tensor>(L"g", bra{L"i_2"}, ket{a}), middle},
        Product::Flatten::No);
  };
  auto depth = [](const ExprPtr& expr) {
    std::size_t result = 0;
    for (const ExprPtr* current = &expr;;) {
      const auto& factors = (*current)->as<Product>().factors();
      const auto it = std::ranges::find_if(
          factors, [](const ExprPtr& factor) { return factor->is<Product>(); });
      if (it == factors.end()) return result;
      ++result;
      current = &*it;
    }
  };

  const auto input = make(L"a_3", L"a_4", L"i_3");
  auto expr = input->clone();
  canonicalize(expr, deduce_named_indices);
  // the parenthesization is preserved
  CHECK(depth(expr) == 2);
  // the dummy indices of the nested products are named consistently: flat,
  // it is the same tensor network
  CHECK_THAT(flatten(expr->clone()), EquivalentTo(flatten(input->clone())));
  // the canonical form does not depend on the dummy labels
  auto other = make(L"a_9", L"a_8", L"i_7");
  canonicalize(other, deduce_named_indices);
  CHECK(*other == *expr);
  // and is marked as such
  CHECK(expr->is_canonical(deduce_named_indices));
  const auto copy = expr->clone();
  canonicalize(expr, deduce_named_indices);
  CHECK(*expr == *copy);

  SECTION("a symmetry of the network that the parenthesization breaks") {
    // X X T2 is symmetric under exchanging the X (T2 being column-symmetric),
    // X (X T2) is not: the two spellings are the same tree, which must not
    // depend on which X the canonicalization of the network puts first
    auto make_tree = [](const wchar_t* outer, const wchar_t* inner) {
      auto t2 = ex<Tensor>(L"T2", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                           Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                           ColumnSymmetry::Symm);
      auto x = [](const wchar_t* n) {
        return ex<Tensor>(L"X", bra{std::wstring(L"i_") + n},
                          ket{std::wstring(L"a_") + n});
      };
      return ex<Product>(
          ExprPtrList{x(outer), ex<Product>(ExprPtrList{x(inner), t2},
                                            Product::Flatten::No)},
          Product::Flatten::No);
    };
    auto first = make_tree(L"1", L"2");
    auto second = make_tree(L"2", L"1");
    canonicalize(first, deduce_named_indices);
    canonicalize(second, deduce_named_indices);
    CHECK(*first == *second);
    CHECK(depth(first) == 1);
  }
}

TEST_CASE("prettify", "[algorithms][canonicalize]") {
  using namespace sequant;

  auto labels_of = [](const ExprPtr& product) {
    REQUIRE(product->is<Product>());
    std::vector<std::wstring> labels;
    for (const auto& factor : product->as<Product>().factors())
      labels.emplace_back(factor->as<Tensor>().label());
    return labels;
  };
  auto make = [] {
    return ex<Tensor>(L"Z", bra{L"i_1"}, ket{L"a_1"}) *
           ex<Tensor>(L"Y", bra{L"a_1"}, ket{L"i_1"});
  };
  const std::vector<std::wstring> y_first{L"Y", L"Z"};
  const std::vector<std::wstring> z_first{L"Z", L"Y"};
  const PrettifyOptions z_opts{
      .presentation = {.cardinal_tensor_labels =
                           container::vector<std::wstring>{L"Z"}}};

  SECTION("sort only") {
    const auto input = make();
    const auto pretty = prettify(input, deduce_named_indices,
                                 PrettifyOptions{.recanonicalize = false});
    CHECK(labels_of(pretty) == y_first);
    // the input is not modified
    CHECK(labels_of(input) == z_first);
    CHECK(labels_of(
              prettify(input, deduce_named_indices,
                       PrettifyOptions{.recanonicalize = false,
                                       .presentation = z_opts.presentation})) ==
          z_first);
    // only sorted: the dummy indices keep their labels
    const auto sum = deserialize(L"t{a_3;i_3} + f{a_7;i_7}");
    CHECK(*prettify(sum, deduce_named_indices,
                    PrettifyOptions{.recanonicalize = false}) ==
          *deserialize(L"f{a_7;i_7} + t{a_3;i_3}"));
  }

  SECTION("full") {
    // with the default options the order of the library's canonical form
    CHECK(labels_of(prettify(make(), deduce_named_indices)) == y_first);
    CHECK(labels_of(prettify(make(), deduce_named_indices, z_opts)) == z_first);
    // an index comparer of the presentation settings
    const auto default_cmp = TensorCanonicalizer::default_index_comparer();
    const PrettifyOptions reversed{
        .presentation = {
            .index_comparer = std::make_shared<const tensor_index_comparer_t>(
                [default_cmp](const Index& idx1, const Index& idx2) {
                  return default_cmp(idx2, idx1);
                })}};
    const auto q = ex<Tensor>(L"Q", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                              Symmetry::Antisymm);
    const auto pretty = prettify(q, deduce_named_indices, reversed);
    CHECK(*pretty == *ex<Tensor>(L"Q", bra{L"i_2", L"i_1"}, ket{L"a_2", L"a_1"},
                                 Symmetry::Antisymm));
    CHECK_THAT(pretty, EquivalentTo(q));
  }

  SECTION("equivalence and idempotence") {
    for (const auto* input :
         {L"1/2 g{i_1,i_2;a_3,a_4}:A * t{a_3,a_4;i_1,i_2}:A + f{i_1;a_1} * "
          L"t{a_1;i_1}",
          L"R{i_3;a_3} * t{a_3;i_3} + 2 f{i_1;a_1} * t{a_1;i_1}"}) {
      CAPTURE(toUtf8(input));
      const auto expr = deserialize(input);
      for (const auto& opts :
           {PrettifyOptions{}, z_opts, mbpt::prettify_options()}) {
        const auto pretty = prettify(expr, deduce_named_indices, opts);
        CHECK_THAT(pretty, EquivalentTo(expr));
        CHECK(prettify(pretty, deduce_named_indices, opts) == pretty);
      }
    }
  }

  SECTION("named indices and the result are preserved") {
    const auto result = deserialize<ResultExpr>(
        L"R{a_1;i_1} = f{a_2;i_1} * t{a_1;i_2} * X{i_2;a_2} + "
        L"t{a_1;i_1} * f{i_2;a_2} * t{a_2;i_2}");
    const auto pretty = prettify(result, mbpt::prettify_options());
    CHECK(pretty.has_label());
    CHECK(pretty.label() == result.label());
    CHECK(pretty.bra() == result.bra());
    CHECK(pretty.ket() == result.ket());
    // the external indices keep their labels
    auto unique_indices = [](const ExprPtr& expr) {
      const auto groups = get_unique_indices(expr);
      container::set<Index> result;
      for (const auto* group : {&groups.bra, &groups.ket, &groups.aux})
        result.insert(group->begin(), group->end());
      return result;
    };
    CHECK(unique_indices(pretty.expression()) ==
          unique_indices(result.expression()));
    CHECK(unique_indices(pretty.expression()).size() == 2);
    CHECK_THAT(pretty, EquivalentTo(result));
  }

  SECTION("the mbpt preset") {
    // t precedes R among the mbpt labels, R precedes t lexicographically
    const auto expr = ex<Tensor>(L"R", bra{L"i_1"}, ket{L"a_1"}) *
                      ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
    CHECK(labels_of(prettify(expr, deduce_named_indices)) ==
          std::vector<std::wstring>{L"R", L"t"});
    CHECK(labels_of(
              prettify(expr, deduce_named_indices, mbpt::prettify_options())) ==
          std::vector<std::wstring>{L"t", L"R"});
  }

  SECTION("the parenthesization is preserved") {
    // (f t) (f t), the subproducts having distinct dummy indices that
    // canonicalizing each subproduct on its own would name alike
    auto make_nested = [](const wchar_t* a, const wchar_t* b,
                          const wchar_t* i) {
      auto left =
          ex<Product>(ExprPtrList{ex<Tensor>(L"f", bra{i}, ket{a}),
                                  ex<Tensor>(L"t", bra{a}, ket{L"i_1"})},
                      Product::Flatten::No);
      auto right =
          ex<Product>(ExprPtrList{ex<Tensor>(L"f", bra{L"i_4"}, ket{b}),
                                  ex<Tensor>(L"t", bra{b}, ket{i})},
                      Product::Flatten::No);
      return ex<Product>(ExprPtrList{left, right}, Product::Flatten::No);
    };
    auto is_nested = [](const ExprPtr& expr) {
      REQUIRE(expr->is<Product>());
      REQUIRE(expr->as<Product>().factors().size() == 2);
      for (const auto& factor : expr->as<Product>().factors()) {
        REQUIRE(factor->is<Product>());
        REQUIRE(factor->as<Product>().factors().size() == 2);
      }
      return true;
    };
    const auto input = make_nested(L"a_3", L"a_4", L"i_3");
    for (const auto& opts : {PrettifyOptions{}, mbpt::prettify_options()}) {
      const auto pretty = prettify(input, deduce_named_indices, opts);
      CHECK(is_nested(pretty));
      // flat, it is the same tensor network
      CHECK_THAT(flatten(pretty->clone()),
                 EquivalentTo(flatten(input->clone())));
      CHECK(*prettify(pretty, deduce_named_indices, opts) == *pretty);
      // the dummy labels do not matter
      CHECK(*prettify(make_nested(L"a_7", L"a_5", L"i_9"), deduce_named_indices,
                      opts) == *pretty);
    }

    // a product with a sum as a factor is only sorted
    const auto with_sum =
        ex<Product>(ExprPtrList{ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"}),
                                ex<Tensor>(L"g", bra{L"i_1"}, ket{L"a_3"}) +
                                    ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_3"})},
                    Product::Flatten::No);
    const auto pretty = prettify(with_sum, deduce_named_indices);
    REQUIRE(pretty->is<Product>());
    REQUIRE(pretty->as<Product>().factors().size() == 2);
    const auto& factors = pretty->as<Product>().factors();
    const auto sum_it = std::ranges::find_if(
        factors, [](const ExprPtr& factor) { return factor->is<Sum>(); });
    REQUIRE(sum_it != factors.end());
    // the summands are sorted, the dummy indices keep their labels
    CHECK(**sum_it == *(ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_3"}) +
                        ex<Tensor>(L"g", bra{L"i_1"}, ket{L"a_3"})));
  }

  SECTION("terms related by the symmetry of the result") {
    // the two terms are images of each other under the permutation of the
    // columns of the result; but a term is not replaced by its image unless it
    // is under a symmetrizer: that the result has the symmetry does not make
    // its terms have it
    const auto first = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-C-S = X{a_1;i_1} * Y{a_2;i_2}");
    const auto second = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-C-S = X{a_2;i_2} * Y{a_1;i_1}");
    for (const auto& opts : {PrettifyOptions{}, mbpt::prettify_options()}) {
      CHECK(*prettify(first, opts).expression() == *first.expression());
      CHECK(*prettify(second, opts).expression() == *second.expression());
    }
    // under a symmetrizer each summand is presented by its own image
    const auto symmetrized = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-C-S = Ŝ{i_1,i_2;a_1,a_2}:N-N-S * "
        L"(X{a_1;i_1} * Y{a_2;i_2} + Z{a_1;i_1} * W{a_2;i_2})");
    const auto symmetrized_image = deserialize<ResultExpr>(
        L"R{a_1,a_2;i_1,i_2}:N-C-S = Ŝ{i_1,i_2;a_1,a_2}:N-N-S * "
        L"(X{a_1;i_1} * Y{a_2;i_2} + Z{a_2;i_2} * W{a_1;i_1})");
    REQUIRE(symmetrized.expression()->is<Product>());
    for (const auto& opts : {PrettifyOptions{}, mbpt::prettify_options()}) {
      const auto pretty = prettify(symmetrized, opts);
      CHECK(*pretty.expression() ==
            *prettify(symmetrized_image, opts).expression());
      // the structure is kept
      CHECK(pretty.expression()->is<Product>());
    }
  }

  SECTION("the result is not canonical") {
    auto expr = make();
    canonicalize(expr, deduce_named_indices);
    REQUIRE(expr->is_canonical(deduce_named_indices));
    const auto pretty = prettify(expr, deduce_named_indices, z_opts);
    CHECK(labels_of(pretty) == z_first);
    CHECK(!pretty->is_canonical(deduce_named_indices));
    // the input keeps its mark
    CHECK(expr->is_canonical(deduce_named_indices));
  }
}
