#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/v1.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <SeQuant/core/bliss.hpp>

#include <memory>
#include <string>
#include <type_traits>

#include <range/v3/all.hpp>

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
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_5"}, ket{L"a_1"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_5", L"a_2"},
                     Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo("S{a1,a2;i1,i2} f{a3;i3} t{i3;a2} t{i1,i2;a1,a3}"));
    }
    {
      auto input =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo(
              "S{a_1,a_2;i_1,i_2} f{a_3;i_3} t{i_2;a_3} t{i_1,i_3;a_1,a_2}"));
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
            parse_expr(L"1/2 t{a3,a1,a2;i4,i5,i2}:N-C-S g{i4,i5;i3,i1}:N-C-S");
        //      auto input1 = parse_expr(L"1/2 t{a1,a2,a3;i5,i2,i4}:N-C-S
        //      g{i4,i5;i3,i1}:N-C-S");
        auto input2 =
            parse_expr(L"1/2 t{a1,a3,a2;i5,i4,i2}:N-C-S g{i5,i4;i1,i3}:N-C-S");
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
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Nonsymm) *
          q2 * ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Variable>(L"p") *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) *
          ex<Variable>(L"q1") *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::Nonsymm);
      canonicalize(input);
      REQUIRE_THAT(input,
                   SimplifiesTo("p q1 q2^* S{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "t{i_2;a_3} t{i_1,i_3;a_1,a_2}"));
    }
    {  // Product containing adjoint of a Tensor
      auto f2 = ex<Tensor>(L"f", bra{L"a_1", L"a_2"}, ket{L"i_5", L"i_2"},
                           Symmetry::Nonsymm, BraKetSymmetry::Nonsymm);
      f2->adjoint();
      auto input1 =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) * f2;
      canonicalize(input1);
      REQUIRE_THAT(input1,
                   SimplifiesTo("S{a_1,a_2;i_1,i_2} f{a_3;i_3} "
                                "f⁺{i_1,i_3;a_1,a_2}:N-N-S t{i_2;a_3}"));
      auto input2 =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::Nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::Nonsymm) * f2 *
          ex<Variable>(L"w") * ex<Constant>(rational{1, 2});
      canonicalize(input2);
      REQUIRE_THAT(input2,
                   SimplifiesTo("1/2 w S{a_1,a_2;i_1,i_2} f{a_3;i_3} "
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
      Context ctx = get_default_context();
      ctx.set(BraKetSymmetry::Symm);
      auto resetter = set_scoped_default_context(ctx);
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
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_2", L"i_4"}, Symmetry::Nonsymm) +
            ex<Constant>(-4) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_4", L"i_3"}, Symmetry::Nonsymm);
        canonicalize(input);
        REQUIRE_THAT(
            input,
            EquivalentTo(
                "-8 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2}"));
      }

      {
        auto term1 =
            ex<Constant>(-4) *
            ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_3", L"i_2", L"i_4"}, Symmetry::Nonsymm);
        auto term2 =
            ex<Constant>(-4) *
            ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_2", L"i_4", L"i_3"}, Symmetry::Nonsymm);
        canonicalize(term1);
        canonicalize(term2);
        REQUIRE_THAT(term1,
                     EquivalentTo("-4 S{i_1,i_2,i_3;a_1,a_2,a_3} f{i_4;i_3} "
                                  "t{a_1,a_2,a_3;i_1,i_4,i_2}"));
        REQUIRE_THAT(term2,
                     EquivalentTo("-4 S{i_1,i_2,i_3;a_1,a_2,a_3} f{i_4;i_3} "
                                  "t{a_1,a_2,a_3;i_1,i_4,i_2}"));
        auto sum_of_terms = term1 + term2;
        simplify(sum_of_terms);
        REQUIRE_THAT(
            sum_of_terms,
            EquivalentTo(
                "-8 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2}"));
      }

      {  // Terms 2 and 4 from spin-traced result
        auto input =
            ex<Constant>(2) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_4", L"i_2"}, Symmetry::Nonsymm) +
            ex<Constant>(2) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Nonsymm);
        canonicalize(input);
        REQUIRE_THAT(
            input, EquivalentTo(
                       "4 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i4,i1,i2}"));
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
