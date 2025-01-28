#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <SeQuant/core/bliss.hpp>

#include <memory>
#include <string>
#include <type_traits>

#include <range/v3/all.hpp>

TEST_CASE("Canonicalizer", "[algorithms]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  auto ctx_resetter = set_scoped_default_context(
      Context(sequant::mbpt::make_legacy_spaces(), Vacuum::SingleProduct));

  SECTION("tests") {
    // Case 7: with protoindices
    {
      auto& l = Logger::instance();
      //      l.tensor_network = l.canonicalize = l.canonicalize_dot =
      //          l.canonicalize_input_graph = true;

      if (true) {
        auto input1 = L"g{i3,i4;a3<i1,i4>,a4<i2>} * s{a1<i1,i2>;a5<i3>}";
        auto input2 = L"g{i3,i4;a3<i1,i3>,a4<i2>} * s{a2<i1,i2>;a6<i4>}";
        auto input3 = L"g{i3,i4;a3<i1,i4>,a4<i2>} * s{a2<i1,i2>;a6<i3>}";
        auto input4 = L"g{i3,i4;a3<i1,i3>,a4<i2>} * s{a2<i1,i2>;a6<i4>}";

        // input1 with transposed tensors
        auto input1_ref = L"s{a1<i1,i2>;a5<i3>} * g{i3,i4;a3<i1,i4>,a4<i2>}";

        std::wcout << "============== " << input1_ref
                   << " ===============" << std::endl;
        TensorNetwork tn1(parse_expr(input1_ref).as<Product>());
        TensorNetwork::CanonicalizationMetadata cbp1;
        auto tn1_canon_factor = tn1.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(), false, nullptr,
            /* rename_named_indices = */ true, &cbp1);

        for (auto&& [orig_idx, new_idx] : cbp1.named_index_replacements) {
          std::wcout << "named index " << orig_idx.to_latex()
                     << " should be renamed to " << new_idx.to_latex()
                     << " in canonical form" << std::endl;
        }

        for (auto& input : {input1, input2, input3, input4}) {
          std::wcout << "============== " << input
                     << " ===============" << std::endl;
          TensorNetwork tn(parse_expr(input).as<Product>());
          TensorNetwork::CanonicalizationMetadata cbp;
          auto tn_canon_factor = tn.canonicalize(
              TensorCanonicalizer::cardinal_tensor_labels(), false, nullptr,
              /* rename_named_indices = */ true, &cbp);

          for (auto&& [orig_idx, new_idx] : cbp.named_index_replacements) {
            std::wcout << "named index " << orig_idx.to_latex()
                       << " should be renamed to " << new_idx.to_latex()
                       << " in canonical form" << std::endl;
          }

          std::wcout << "graph(" << input1_ref << ") <=> graph(" << input
                     << "): " << cbp1.graph->cmp(*cbp.graph) << std::endl;
        }
      }

      {
        auto input1 =
            L"s{a2<i1,i2>;a6<i2,i4>} * g{i3,i4;a3<i2,i4>,a4<i1,i3>} * "
            L"t{a3<i2,i4>,a6<i2,i4>;i4,i2}";
        auto input2 =
            L"g{i3,i4;a3<i1,i4>,a4<i2,i3>} * t{a3<i1,i4>,a5<i1,i4>;i4,i1} * "
            L"s{a1<i1,i2>;a5<i1,i4>}";

        std::wcout << "============== " << input1
                   << " ===============" << std::endl;
        auto ex1 = parse_expr(input1);
        TensorNetwork tn1(ex1.as<Product>());
        TensorNetwork::CanonicalizationMetadata cbp1;
        auto tn1_canon_factor = tn1.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(), false, nullptr,
            /* rename_named_indices = */ true, &cbp1);
        for (auto&& [orig_idx, new_idx] : cbp1.named_index_replacements) {
          std::wcout << "named index " << orig_idx.to_latex()
                     << " should be renamed to " << new_idx.to_latex()
                     << " in canonical form" << std::endl;
        }

        std::wcout << "============== " << input2
                   << " ===============" << std::endl;
        auto ex2 = parse_expr(input2);
        TensorNetwork tn2(ex2.as<Product>());
        TensorNetwork::CanonicalizationMetadata cbp2;
        auto tn2_canon_factor = tn2.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(), false, nullptr,
            /* rename_named_indices = */ true, &cbp2);

        for (auto&& [orig_idx, new_idx] : cbp2.named_index_replacements) {
          std::wcout << "named index " << orig_idx.to_latex()
                     << " should be renamed to " << new_idx.to_latex()
                     << " in canonical form" << std::endl;
        }

        std::wcout << "graph(" << input1 << ") <=> graph(" << input2
                   << "): " << cbp1.graph->cmp(*cbp2.graph) << std::endl;

        REQUIRE(cbp1.graph->cmp(*cbp2.graph) == 0);

        //        std::wcout << canonicalize(ex1).to_latex() << " should be
        //        equal " << canonicalize(ex2).to_latex() << std::endl;
      }

      //      l.tensor_network = l.canonicalize = l.canonicalize_dot =
      //          l.canonicalize_input_graph = false;
    }
  }

  SECTION("Tensors") {
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                           Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_3", L"p_4"},
                           Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p4,p3}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                           Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::symm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("g{p1,p2;p3,p4}:S"));
    }
    {
      auto op = ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           Symmetry::antisymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("-g{p1,p2;p3,p4}:A"));
    }

    // aux indices
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_2"}, aux{L"p_3"},
                           Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1;p2;p3}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1,p2;p4,p3;p5}"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::symm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("B{p1,p2;p3,p4;p5}:S"));
    }
    {
      auto op = ex<Tensor>(L"B", bra{L"p_1", L"p_2"}, ket{L"p_4", L"p_3"},
                           aux{L"p_5"}, Symmetry::antisymm);
      canonicalize(op);
      REQUIRE_THAT(op, SimplifiesTo("-B{p1,p2;p3,p4;p5}:A"));
    }
  }

  SECTION("Products") {
    {
      auto input =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", bra{L"i_5"}, ket{L"a_1"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1", L"i_2"}, ket{L"a_5", L"a_2"},
                     Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo("S{a1,a2;i1,i2} f{a3;i3} t{i3;a2} t{i1,i2;a1,a3}"));
    }
    {
      auto input =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo("S{a1,a3;i1,i2} f{a2;i3} t{i_2;a_2} t{i1,i3;a1,a3}"));
    }
    {  // Product containing Variables
      auto q2 = ex<Variable>(L"q2");
      q2->adjoint();
      auto input =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::nonsymm) *
          q2 * ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::nonsymm) *
          ex<Variable>(L"p") *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::nonsymm) *
          ex<Variable>(L"q1") *
          ex<Tensor>(L"t", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo(
              "p q1 q2^* S{a1,a3;i1,i2} f{a2;i3} t {i2;a2} t{i1,i3;a1,a3}"));
    }
    {  // Product containing adjoint of a Tensor
      auto f2 = ex<Tensor>(L"f", bra{L"i_5", L"i_2"}, ket{L"a_1", L"a_2"},
                           Symmetry::nonsymm, BraKetSymmetry::nonsymm);
      f2->adjoint();
      auto input1 =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::nonsymm) * f2;
      canonicalize(input1);
      REQUIRE_THAT(
          input1,
          SimplifiesTo("S{a1,a2;i1,i3} f{a3;i2} f⁺{a1,a2;i1,i2} t{i3;a3}"));
      auto input2 =
          ex<Tensor>(L"S", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"f", bra{L"a_5"}, ket{L"i_5"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_5"}, Symmetry::nonsymm) * f2 *
          ex<Variable>(L"w") * ex<Constant>(rational{1, 2});
      canonicalize(input2);
      REQUIRE_THAT(
          input2,
          SimplifiesTo(
              "1/2 w S{a1,a2;i1,i3} f{a3;i2} f⁺{a1,a2;i1,i2} t{i3;a3}"));
    }
  }
  {
    auto input = ex<Constant>(rational{1, 2}) *
                 ex<Tensor>(L"B", bra{L"p_2"}, ket{L"p_4"}, aux{L"p_5"},
                            Symmetry::nonsymm) *
                 ex<Tensor>(L"B", bra{L"p_1"}, ket{L"p_3"}, aux{L"p_5"},
                            Symmetry::nonsymm) *
                 ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm) *
                 ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm);
    canonicalize(input);
    REQUIRE_THAT(input,
                 SimplifiesTo("1/2 t{p1;p3} t{p2;p4} B{p3;p1;p5} B{p4;p2;p5}"));
  }

  SECTION("sum of products") {
    {
      // CASE 1: Non-symmetric tensors
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm);
      simplify(input);
      canonicalize(input);
      REQUIRE_THAT(input, SimplifiesTo("g{p3,p4;p1,p2} t{p1;p3} t{p2;p4}"));
    }

    // CASE 2: Symmetric tensors
    {
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::symm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::symm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE_THAT(input, SimplifiesTo("g{p2,p3;p1,p4}:S t{p1;p2} t{p4;p3}"));
    }

    // Case 3: Anti-symmetric tensors
    {
      auto input =
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm) +
          ex<Constant>(rational{1, 2}) *
              ex<Tensor>(L"g", bra{L"p_2", L"p_1"}, ket{L"p_4", L"p_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", bra{L"p_3"}, ket{L"p_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"p_4"}, ket{L"p_2"}, Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE_THAT(input, SimplifiesTo("g{p2,p3;p1,p4}:A t{p1;p2} t{p4;p3}"));
    }

    // Case 4: permuted indices
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"i_1"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::antisymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         Symmetry::antisymm);
      canonicalize(input);
      REQUIRE(input->size() == 1);
      REQUIRE_THAT(input,
                   SimplifiesTo("g{i3,i4;i1,a3}:A t{a2;i3} t{a1,a3;i2,i4}:A"));
    }

    // Case 4: permuted indices from CCSD R2 biorthogonal configuration
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"i_1"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::nonsymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         Symmetry::nonsymm);

      canonicalize(input);
      REQUIRE(input->size() == 1);
      REQUIRE_THAT(input,
                   SimplifiesTo("g{i3,i4;i1,a3} t{a2;i4} t{a1,a3;i3,i2}"));
    }

    {  // Case 5: CCSDT R3: S3 * F * T3

      {  // Terms 1 and 6 from spin-traced result
        auto input =
            ex<Constant>(-4) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_2", L"i_4"}, Symmetry::nonsymm) +
            ex<Constant>(-4) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_4", L"i_3"}, Symmetry::nonsymm);
        canonicalize(input);
        REQUIRE_THAT(
            input,
            SimplifiesTo(
                "-8 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2}"));
      }

      {
        auto term1 =
            ex<Constant>(-4) *
            ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_3", L"i_2", L"i_4"}, Symmetry::nonsymm);
        auto term2 =
            ex<Constant>(-4) *
            ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
            ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_2", L"i_4", L"i_3"}, Symmetry::nonsymm);
        canonicalize(term1);
        canonicalize(term2);
        REQUIRE_THAT(
            term1,
            SimplifiesTo(
                "-4 S{i1,i3,i4;a1,a2,a3} f{i2;i4} t{a1,a2,a3;i1,i2,i3}"));
        REQUIRE_THAT(
            term2,
            SimplifiesTo(
                "-4 S{i1,i3,i4;a1,a2,a3} f{i2;i4} t{a1,a2,a3;i1,i2,i3}"));
        auto sum_of_terms = term1 + term2;
        simplify(sum_of_terms);
        REQUIRE_THAT(
            sum_of_terms,
            SimplifiesTo(
                "-8 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2}"));
      }

      {  // Terms 2 and 4 from spin-traced result
        auto input =
            ex<Constant>(2) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_3", L"i_4", L"i_2"}, Symmetry::nonsymm) +
            ex<Constant>(2) *
                ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                           ket{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
                ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                           ket{L"i_2", L"i_3", L"i_4"}, Symmetry::nonsymm);
        canonicalize(input);
        REQUIRE_THAT(
            input, SimplifiesTo(
                       "4 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i4,i1,i2}"));
      }
    }

    // Case 6: Case 4 w/ aux indices
    {
      auto input =
          ex<Constant>(rational{4, 3}) *
              ex<Tensor>(L"B", bra{L"i_3"}, ket{L"a_3"}, aux{L"p_5"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"B", bra{L"i_4"}, ket{L"i_1"}, aux{L"p_5"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_4", L"i_2"},
                         Symmetry::nonsymm) -
          ex<Constant>(rational{1, 3}) *
              ex<Tensor>(L"B", bra{L"i_3"}, ket{L"i_1"}, aux{L"p_5"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"B", bra{L"i_4"}, ket{L"a_3"}, aux{L"p_5"},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_3", L"i_2"},
                         Symmetry::nonsymm);

      canonicalize(input);
      simplify(input);
      REQUIRE_THAT(
          input,
          SimplifiesTo("t{a2;i4} t{a1,a3;i3,i2} B{i3;i1;p5} B{i4;a3;p5}"));
    }
  }
}
