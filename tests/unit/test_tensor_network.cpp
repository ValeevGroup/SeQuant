//
// Created by Eduard Valeyev on 3/23/18.
//

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cmath>
#include <codecvt>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <SeQuant/core/timer.hpp>
#include <range/v3/all.hpp>

using namespace sequant;

TEMPLATE_TEST_CASE("tensor_network_shared", "[elements]", TensorNetwork,
                   TensorNetworkV2) {
  using TN = TestType;

  SECTION("canonicalize_slots") {
    SECTION("phase_difference") {
      const Product prod1 =
          parse_expr(L"g{i_2,i_3;a_2,a_3}:A-C-S * t{a_2;i_2}:A-C-S")
              ->as<Product>();
      const Product prod2 =
          parse_expr(L"g{i_2,i_3;a_2,a_3}:A-C-S * t{a_2;i_3}:A-C-S")
              ->as<Product>();

      TN tn1(prod1.factors());
      TN tn2(prod2.factors());

      const auto& canon1 =
          tn1.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels());
      const auto& canon2 =
          tn2.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels());

      REQUIRE(canon1.hash_value() == canon2.hash_value());
      REQUIRE(canon1.phase != canon2.phase);
    }
  }

  SECTION("canonicalize") {
    SECTION("need particle reorder?") {
      {
        Index::reset_tmp_index();
        auto input1 = parse_expr(
            L"t{p1,p2;p3,p4}:N-C-S t{p4,p5;p6,p7}:N-C-S t{p7,p8;p9,p1}:N-C-S");
        // N.B. renaming external index changes local canonical order produced
        // by TNV1 (due to the use of DefaultTensorCanonicalizer)
        auto input2 = parse_expr(
            L"t{p1,p2;p3,p4}:N-C-S t{p4,p5;p11,p7}:N-C-S t{p7,p8;p9,p1}:N-C-S");
        std::wcout << "input1 = " << to_latex(input1) << std::endl;
        std::wcout << "input2 = " << to_latex(input2) << std::endl;

        {
          TN tn1(*input1);
          tn1.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                           false);
          TN tn2(*input2);
          tn2.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                           false);

          // std::wcout << "tn1[0] = " <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn1.tensors()[0])) <<
          // "\n"; std::wcout << "tn1[1] = " <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn1.tensors()[1])) <<
          // "\n"; std::wcout << "tn1[2] = " <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn1.tensors()[2])) <<
          // "\n"; std::wcout << "tn2[0] = " <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn2.tensors()[0])) <<
          // "\n"; std::wcout << "tn2[1] = " <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn2.tensors()[1])) <<
          // "\n"; std::wcout << "tn2[2] = " <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn2.tensors()[2])) <<
          // "\n";

          // TNv1 fails to canonicalize this correctly
          if constexpr (std::is_same_v<TN, TensorNetworkV2>) {
            // input2 obtained from input1 by i6 -> i11, which "frees" i6 for
            // dummy renamings so canonical(input2) is obtained from
            // canonical(input1) by i6 -> i11 and i7 -> i6
            REQUIRE(tn1.tensors()[0]->_to_latex() ==
                    tn2.tensors()[0]->_to_latex());
            REQUIRE(ranges::equal(tn1.tensors()[1]->_bra(),
                                  tn2.tensors()[1]->_bra()));
            REQUIRE(ranges::equal(tn1.tensors()[2]->_ket(),
                                  tn2.tensors()[2]->_ket()));
          }
        }
      }
    }
  }
}

TEST_CASE("tensor_network", "[elements]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using sequant::Context;
  namespace t = sequant::mbpt::tensor;
  namespace o = sequant::mbpt::op;

  SECTION("constructors") {
    {  // with Tensors
      auto t1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"});
      auto t2 = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"i_2"});
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));

      auto t1_x_t2_p_t2 = t1 * (t2 + t2);  // can only use a flat tensor product
      REQUIRE_THROWS_AS(TensorNetwork(*t1_x_t2_p_t2), std::logic_error);
    }

    {  // with NormalOperators
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<FNOperator>(cre({L"i_1"}), ann({L"i_2"}), V);
      auto t2 = ex<FNOperator>(cre({L"i_2"}), ann({L"i_1"}), V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));
    }

    {  // with Tensors and NormalOperators
      auto tmp = t::A(nₚ(-2)) * t::H_(2) * t::T_(2) * t::T_(2);
      REQUIRE_NOTHROW(TensorNetwork(tmp->as<Product>().factors()));
    }

  }  // SECTION("constructors")

  SECTION("accessors") {
    {
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"});
      auto t2 = ex<FNOperator>(cre({L"i_1"}), ann({L"i_3"}), V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));
      TensorNetwork tn(*t1_x_t2);

      // edges
      auto edges = tn.edges();
      REQUIRE(edges.size() == 3);

      // ext indices
      auto ext_indices = tn.ext_indices();
      REQUIRE(ext_indices.size() == 2);

      // tensors
      auto tensors = tn.tensors();
      REQUIRE(size(tensors) == 2);
      REQUIRE(std::dynamic_pointer_cast<Expr>(tensors[0]));
      REQUIRE(std::dynamic_pointer_cast<Expr>(tensors[1]));
      REQUIRE(*std::dynamic_pointer_cast<Expr>(tensors[0]) == *t1);
      REQUIRE(*std::dynamic_pointer_cast<Expr>(tensors[1]) == *t2);

      // index replacements performed by canonicalize() ... since canonicalize()
      // not invoked this is empty
      auto idxrepl = tn.idxrepl();
      REQUIRE(idxrepl.size() == 0);
    }
  }  // SECTION("accessors")

  SECTION("canonicalizer") {
    {
      {  // with no external indices, hence no named indices whatsoever
        Index::reset_tmp_index();
        constexpr const auto V = Vacuum::SingleProduct;
        auto t1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"});
        auto t2 = ex<FNOperator>(cre({L"i_1"}), ann({L"i_2"}), V);
        auto t1_x_t2 = t1 * t2;
        TensorNetwork tn(*t1_x_t2);
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

        REQUIRE(size(tn.tensors()) == 2);
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
        //        std::wcout <<
        //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) <<
        //        std::endl; std::wcout <<
        //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) <<
        //        std::endl;
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                L"{F^{{i_2}}_{{i_1}}}");
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                L"{\\tilde{a}^{{i_1}}_{{i_2}}}");
        REQUIRE(tn.idxrepl().size() == 0);
      }

      {
        Index::reset_tmp_index();
        constexpr const auto V = Vacuum::SingleProduct;
        auto t1 = ex<Tensor>(L"F", bra{L"i_2"}, ket{L"i_17"});
        auto t2 = ex<FNOperator>(cre({L"i_2"}), ann({L"i_3"}), V);
        auto t1_x_t2 = t1 * t2;

        // with all external named indices
        {
          TensorNetwork tn(*t1_x_t2);
          tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

          REQUIRE(size(tn.tensors()) == 2);
          REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
          REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
          // std::wcout <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) <<
          // std::endl; std::wcout <<
          // to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) <<
          // std::endl;
          REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                  L"{\\tilde{a}^{{i_1}}_{{i_3}}}");
          REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                  L"{F^{{i_{17}}}_{{i_1}}}");
        }

        // with explicit named indices
        {
          Index::reset_tmp_index();
          TensorNetwork tn(*t1_x_t2);

          using named_indices_t = TensorNetwork::named_indices_t;
          named_indices_t indices{Index{L"i_17"}};
          tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false,
                          &indices);

          REQUIRE(size(tn.tensors()) == 2);
          REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
          REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
          //        std::wcout <<
          //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]))
          //        << std::endl; std::wcout <<
          //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]))
          //        << std::endl;
          REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                  L"{\\tilde{a}^{{i_1}}_{{i_3}}}");
          REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                  L"{F^{{i_{17}}}_{{i_1}}}");
        }
      }
    }
  }  // SECTION("canonicalizer")

  SECTION("bliss graph") {
    auto ctx_resetter = set_scoped_default_context(
        Context(sequant::mbpt::make_legacy_spaces(), Vacuum::SingleProduct));
    Index::reset_tmp_index();
    // to generate expressions in specified (i.e., platform-independent) manner
    // can't use operator expression (due to unspecified order of evaluation of
    // function arguments), must use initializer list
    auto tmp = ex<Product, std::initializer_list<ExprPtr>>(
        {t::A(nₚ(-2)), t::H_(2), t::T_(2), t::T_(2), t::T_(2)});
    // canonicalize to avoid dependence on the implementation details of
    // mbpt::sr::make_op
    canonicalize(tmp);
    // std::wcout << "A2*H2*T2*T2*T2 = " << to_latex(tmp) << std::endl;
    TensorNetwork tn(tmp->as<Product>().factors());

    // make graph, with all labels
    REQUIRE_NOTHROW(
        tn.make_bliss_graph({.make_labels = true, .make_texlabels = true}));
    auto gdata = tn.make_bliss_graph();
    const auto& [graph, vlabels, vtexlabels, vcolors, vtypes] = gdata;

    // test dot representation
    // N.B. cluster tensor vertices only
    std::basic_ostringstream<wchar_t> oss;
    REQUIRE_NOTHROW(graph->write_dot(
        oss, vlabels, vtexlabels,
        {.vertex_to_subgraph = [&](std::size_t vertex_ordinal) {
          return gdata.vertex_to_tensor_cluster(vertex_ordinal);
        }}));
    // std::wcout << "oss.str() = " << std::endl << oss.str() << std::endl;
    const std::wstring actual = oss.str();
    // clang-format off
    const std::wstring expected = L"graph g {\n"
L"node [ style=filled, penwidth=2, margin=0];\n"
L"v0 [ label=\"i_1\", texlbl=\"${i_1}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v0 -- v21\n"
L"v0 -- v42\n"
L"v1 [ label=\"i_2\", texlbl=\"${i_2}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v1 -- v21\n"
L"v1 -- v42\n"
L"v2 [ label=\"i_3\", texlbl=\"${i_3}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v2 -- v30\n"
L"v2 -- v57\n"
L"v3 [ label=\"i_4\", texlbl=\"${i_4}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v3 -- v30\n"
L"v3 -- v57\n"
L"v4 [ label=\"i_5\", texlbl=\"${i_5}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v4 -- v34\n"
L"v4 -- v53\n"
L"v5 [ label=\"i_6\", texlbl=\"${i_6}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v5 -- v34\n"
L"v5 -- v53\n"
L"v6 [ label=\"i_7\", texlbl=\"${i_7}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v6 -- v38\n"
L"v6 -- v49\n"
L"v7 [ label=\"i_8\", texlbl=\"${i_8}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
L"v7 -- v38\n"
L"v7 -- v49\n"
L"v8 [ label=\"a_1\", texlbl=\"${a_1}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v8 -- v22\n"
L"v8 -- v41\n"
L"v9 [ label=\"a_2\", texlbl=\"${a_2}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v9 -- v22\n"
L"v9 -- v41\n"
L"v10 [ label=\"a_3\", texlbl=\"${a_3}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v10 -- v29\n"
L"v10 -- v58\n"
L"v11 [ label=\"a_4\", texlbl=\"${a_4}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v11 -- v29\n"
L"v11 -- v58\n"
L"v12 [ label=\"a_5\", texlbl=\"${a_5}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v12 -- v33\n"
L"v12 -- v54\n"
L"v13 [ label=\"a_6\", texlbl=\"${a_6}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v13 -- v33\n"
L"v13 -- v54\n"
L"v14 [ label=\"a_7\", texlbl=\"${a_7}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v14 -- v37\n"
L"v14 -- v50\n"
L"v15 [ label=\"a_8\", texlbl=\"${a_8}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
L"v15 -- v37\n"
L"v15 -- v50\n"
L"v16 [ label=\"κ_1\", texlbl=\"${\\kappa_1}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
L"v16 -- v25\n"
L"v16 -- v46\n"
L"v17 [ label=\"κ_2\", texlbl=\"${\\kappa_2}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
L"v17 -- v25\n"
L"v17 -- v46\n"
L"v18 [ label=\"κ_3\", texlbl=\"${\\kappa_3}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
L"v18 -- v26\n"
L"v18 -- v45\n"
L"v19 [ label=\"κ_4\", texlbl=\"${\\kappa_4}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
L"v19 -- v26\n"
L"v19 -- v45\n"
L"subgraph cluster0 {\n"
L"v20 [ label=\"A\", texlbl=\"$A$\", color=\"#257a61\", fillcolor=\"#94f4c2\" ];\n"
L"v20 -- v23\n"
L"v21 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v21 -- v23\n"
L"v22 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v22 -- v23\n"
L"v23 [ label=\"bka\", color=\"#257a61\", fillcolor=\"#94f4c2\" ];\n"
L"}\n"
L"subgraph cluster1 {\n"
L"v24 [ label=\"g\", texlbl=\"$g$\", color=\"#300a49\", fillcolor=\"#c05092\" ];\n"
L"v24 -- v27\n"
L"v25 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v25 -- v27\n"
L"v26 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v26 -- v27\n"
L"v27 [ label=\"bka\", color=\"#300a49\", fillcolor=\"#c05092\" ];\n"
L"}\n"
L"subgraph cluster2 {\n"
L"v28 [ label=\"t\", texlbl=\"$t$\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
L"v28 -- v31\n"
L"v29 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v29 -- v31\n"
L"v30 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v30 -- v31\n"
L"v31 [ label=\"bka\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
L"}\n"
L"subgraph cluster3 {\n"
L"v32 [ label=\"t\", texlbl=\"$t$\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
L"v32 -- v35\n"
L"v33 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v33 -- v35\n"
L"v34 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v34 -- v35\n"
L"v35 [ label=\"bka\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
L"}\n"
L"subgraph cluster4 {\n"
L"v36 [ label=\"t\", texlbl=\"$t$\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
L"v36 -- v39\n"
L"v37 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v37 -- v39\n"
L"v38 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v38 -- v39\n"
L"v39 [ label=\"bka\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
L"}\n"
L"subgraph cluster5 {\n"
L"v40 [ label=\"ã\", texlbl=\"$\\tilde{a}$\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"v40 -- v43\n"
L"v41 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v41 -- v43\n"
L"v42 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v42 -- v43\n"
L"v43 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"}\n"
L"subgraph cluster6 {\n"
L"v44 [ label=\"ã\", texlbl=\"$\\tilde{a}$\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"v44 -- v47\n"
L"v45 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v45 -- v47\n"
L"v46 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v46 -- v47\n"
L"v47 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"}\n"
L"subgraph cluster7 {\n"
L"v48 [ label=\"ã\", texlbl=\"$\\tilde{a}$\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"v48 -- v51\n"
L"v49 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v49 -- v51\n"
L"v50 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v50 -- v51\n"
L"v51 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"}\n"
L"subgraph cluster8 {\n"
L"v52 [ label=\"ã\", texlbl=\"$\\tilde{a}$\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"v52 -- v55\n"
L"v53 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v53 -- v55\n"
L"v54 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v54 -- v55\n"
L"v55 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"}\n"
L"subgraph cluster9 {\n"
L"v56 [ label=\"ã\", texlbl=\"$\\tilde{a}$\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"v56 -- v59\n"
L"v57 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
L"v57 -- v59\n"
L"v58 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
L"v58 -- v59\n"
L"v59 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
L"}\n"
L"}\n";
    // clang-format on

    REQUIRE(actual == expected);

    // make sure can generate without labels also
    REQUIRE_NOTHROW(
        tn.make_bliss_graph({.make_labels = true, .make_texlabels = false})
            .vertex_texlabels.size() == 0);
    REQUIRE_NOTHROW(
        tn.make_bliss_graph({.make_labels = false, .make_texlabels = true})
            .vertex_labels.size() == 0);
    {
      auto g =
          tn.make_bliss_graph({.make_labels = false, .make_texlabels = false});
      REQUIRE_NOTHROW(g.vertex_labels.size() == 0 &&
                      g.vertex_texlabels.size() == 0);
    }

    // compute automorphism group
    {
      bliss::Stats stats;
      graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

      std::vector<std::vector<unsigned int>> aut_generators;
      auto save_aut = [&aut_generators](const unsigned int n,
                                        const unsigned int* aut) {
        aut_generators.emplace_back(aut, aut + n);
      };
      graph->find_automorphisms(stats, &bliss::aut_hook<decltype(save_aut)>,
                                &save_aut);
      std::basic_ostringstream<wchar_t> oss;
      bliss::print_auts(aut_generators, oss, decltype(vlabels){});
      const std::wstring actual = oss.str();
      const std::wstring expected =
          L"(18,19)\n"
          L"(16,17)\n"
          L"(0,1)\n"
          L"(8,9)\n"
          L"(6,7)\n"
          L"(14,15)\n"
          L"(4,5)\n"
          L"(2,3)\n"
          L"(12,13)\n"
          L"(10,11)\n"
          L"(2,4)(3,5)(10,12)(11,13)(28,32)(29,33)(30,34)(31,35)(52,56)(53,57)("
          L"54,58)(55,59)\n"
          L"(4,6)(5,7)(12,14)(13,15)(32,36)(33,37)(34,38)(35,39)(48,52)(49,53)("
          L"50,54)(51,55)\n";

      REQUIRE(actual == expected);

      // change to 1 to user vertex labels rather than indices
      if (0) {
        std::basic_ostringstream<wchar_t> oss2;
        bliss::print_auts(aut_generators, oss2, vlabels);
        std::wcout << oss2.str() << std::endl;
      }
    }

  }  // SECTION("bliss graph")

  SECTION("misc1") {
    if (false) {
      Index::reset_tmp_index();
      // TN1 from manuscript
      auto g = ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                          Symmetry::antisymm);
      auto ta = ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_1", L"i_2"},
                           Symmetry::antisymm);
      auto tb = ex<Tensor>(L"t", bra{L"a_2", L"a_4"}, ket{L"i_3", L"i_4"},
                           Symmetry::antisymm);

      auto tmp = g * ta * tb;
      // std::wcout << "TN1 = " << to_latex(tmp) << std::endl;
      TensorNetwork tn(tmp->as<Product>().factors());

      // make graph
      // N.B. treat all indices as dummy so that the automorphism ignores the
      using named_indices_t = TensorNetwork::named_indices_t;
      named_indices_t indices{};
      REQUIRE_NOTHROW(tn.make_bliss_graph({.named_indices = &indices}));
      auto [graph, vlabels, vtexlabels, vcolors, vtypes] =
          tn.make_bliss_graph({.named_indices = &indices});

      // create dot
      {
        std::basic_ostringstream<wchar_t> oss;
        REQUIRE_NOTHROW(graph->write_dot(oss, vlabels));
        // std::wcout << "oss.str() = " << std::endl << oss.str() << std::endl;
      }

      bliss::Stats stats;
      graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

      std::vector<std::vector<unsigned int>> aut_generators;
      auto save_aut = [&aut_generators](const unsigned int n,
                                        const unsigned int* aut) {
        aut_generators.emplace_back(aut, aut + n);
      };
      graph->find_automorphisms(stats, &bliss::aut_hook<decltype(save_aut)>,
                                &save_aut);
      CHECK(aut_generators.size() ==
            2);  // there are 2 generators, i1<->i2, i3<->i4

      std::basic_ostringstream<wchar_t> oss;
      bliss::print_auts(aut_generators, oss, vlabels);
      CHECK(oss.str() == L"({i_3},{i_4})\n({i_1},{i_2})\n");
      // std::wcout << oss.str() << std::endl;
    }
  }  // SECTION("misc1")
}

template <typename Container>
std::vector<sequant::ExprPtr> to_tensors(const Container& cont) {
  std::vector<sequant::ExprPtr> tensors;

  std::transform(cont.begin(), cont.end(), std::back_inserter(tensors),
                 [](const auto& tensor) {
                   auto casted =
                       std::dynamic_pointer_cast<sequant::Expr>(tensor);
                   REQUIRE(casted != nullptr);
                   return casted;
                 });
  return tensors;
}

template <typename Container>
sequant::ExprPtr to_product(const Container& container) {
  return sequant::ex<sequant::Product>(to_tensors(container));
}

namespace sequant {
class TensorNetworkV2Accessor {
 public:
  auto get_canonical_bliss_graph(
      sequant::TensorNetworkV2 tn,
      const sequant::TensorNetwork::named_indices_t* named_indices = nullptr) {
    tn.canonicalize_graph(named_indices ? *named_indices : tn.ext_indices_);
    tn.init_edges();
    auto graph = tn.create_graph(
        {.named_indices = named_indices,
         .make_labels = Logger::instance().canonicalize_dot,
         .make_texlabels = Logger::instance().canonicalize_dot});
    return std::make_pair(std::move(graph.bliss_graph), graph.vertex_labels);
  }
};
}  // namespace sequant

TEST_CASE("tensor_network_v2", "[elements]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using sequant::Context;
  namespace t = sequant::mbpt::tensor;
  namespace o = sequant::mbpt::op;

  auto ctx_resetter = sequant::set_scoped_default_context(Context(
      mbpt::make_sr_spaces(), Vacuum::SingleProduct, IndexSpaceMetric::Unit,
      BraKetSymmetry::conjugate, SPBasis::spinorbital));

  SECTION("Edges") {
    using Vertex = TensorNetworkV2::Vertex;
    using Edge = TensorNetworkV2::Edge;
    using Origin = TensorNetworkV2::Origin;

    Vertex v1(Origin::Bra, 0, 1, Symmetry::antisymm);
    Vertex v2(Origin::Bra, 0, 0, Symmetry::antisymm);
    Vertex v3(Origin::Ket, 1, 0, Symmetry::symm);
    Vertex v4(Origin::Ket, 1, 3, Symmetry::symm);
    Vertex v5(Origin::Bra, 3, 0, Symmetry::nonsymm);
    Vertex v6(Origin::Bra, 3, 2, Symmetry::nonsymm);
    Vertex v7(Origin::Ket, 3, 1, Symmetry::nonsymm);
    Vertex v8(Origin::Ket, 5, 0, Symmetry::symm);

    const Index dummy(L"a_1");

    Edge e1(v1, dummy);
    e1.connect_to(v4);
    Edge e2(v2, dummy);
    e2.connect_to(v3);
    Edge e3(v3, dummy);
    e3.connect_to(v5);
    Edge e4(v4, dummy);
    e4.connect_to(v6);

    Edge e5(v8, dummy);
    e5.connect_to(v6);
    Edge e6(v8, dummy);
    REQUIRE_THROWS_AS(e6.connect_to(v7), std::logic_error);

    // Due to tensor symmetries, these edges are considered equal
    REQUIRE(e1 == e2);
    REQUIRE(!(e1 < e2));
    REQUIRE(!(e2 < e1));

    // Smallest terminal index wins
    REQUIRE(!(e1 == e3));
    REQUIRE(e1 < e3);
    REQUIRE(!(e3 < e1));

    // For non-symmetric tensors the connection slot is taken into account
    REQUIRE(!(e3 == e4));
    REQUIRE(e3 < e4);
    REQUIRE(!(e4 < e3));

    // Unconnected edges always come before fully connected ones
    REQUIRE(!(e6 == e1));
    REQUIRE(e6 < e1);
    REQUIRE(!(e1 < e6));
  }

  SECTION("constructors") {
    {  // with Tensors
      auto t1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"});
      auto t2 = ex<Tensor>(L"t", bra{L"i_2"}, ket{L"i_1"});
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetworkV2(*t1_x_t2));

      auto t1_x_t2_p_t2 = t1 * (t2 + t2);  // can only use a flat tensor product
      REQUIRE_THROWS_AS(TensorNetworkV2(*t1_x_t2_p_t2), std::logic_error);

      // must be covariant: no bra to bra or ket to ket
      t2->adjoint();
      auto t1_x_t2_adjoint = t1 * t2;
      REQUIRE_THROWS_AS(TensorNetworkV2(t1_x_t2_adjoint), std::logic_error);
    }

    {  // with NormalOperators
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<FNOperator>(cre({L"i_1"}), ann({L"i_2"}), V);
      auto t2 = ex<FNOperator>(cre({L"i_2"}), ann({L"i_1"}), V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetworkV2(*t1_x_t2));
    }

    {  // with Tensors and NormalOperators
      auto tmp = t::A(nₚ(-2)) * t::H_(2) * t::T_(2) * t::T_(2);
      REQUIRE_NOTHROW(TensorNetworkV2(tmp->as<Product>().factors()));
    }

  }  // SECTION("constructors")

  SECTION("accessors") {
    {
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"});
      auto t2 = ex<FNOperator>(cre({L"i_1"}), ann({L"i_3"}), V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetworkV2(*t1_x_t2));
      TensorNetworkV2 tn(*t1_x_t2);

      // edges
      auto edges = tn.edges();
      REQUIRE(edges.size() == 3);

      // ext indices
      auto ext_indices = tn.ext_indices();
      REQUIRE(ext_indices.size() == 2);

      // tensors
      auto tensors = tn.tensors();
      REQUIRE(size(tensors) == 2);
      REQUIRE(std::dynamic_pointer_cast<Expr>(tensors[0]));
      REQUIRE(std::dynamic_pointer_cast<Expr>(tensors[1]));
      REQUIRE(*std::dynamic_pointer_cast<Expr>(tensors[0]) == *t1);
      REQUIRE(*std::dynamic_pointer_cast<Expr>(tensors[1]) == *t2);
    }
  }  // SECTION("accessors")

  SECTION("canonicalizer") {
    {  // with no external indices, hence no named indices whatsoever
      Index::reset_tmp_index();
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<Tensor>(L"F", bra{L"i_1"}, ket{L"i_2"});
      auto t2 = ex<FNOperator>(cre({L"i_1"}), ann({L"i_2"}), V);
      auto t1_x_t2 = t1 * t2;
      TensorNetworkV2 tn(*t1_x_t2);
      tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

      REQUIRE(size(tn.tensors()) == 2);
      REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
      REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
      //        std::wcout <<
      //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) <<
      //        std::endl; std::wcout <<
      //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) <<
      //        std::endl;
      REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
              L"{F^{{i_2}}_{{i_1}}}");
      REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
              L"{\\tilde{a}^{{i_1}}_{{i_2}}}");
    }

    {
      Index::reset_tmp_index();
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<Tensor>(L"F", bra{L"i_2"}, ket{L"i_17"});
      auto t2 = ex<FNOperator>(cre({L"i_2"}), ann({L"i_3"}), V);
      auto t1_x_t2 = t1 * t2;

      // with all external named indices
      SECTION("implicit") {
        TensorNetworkV2 tn(*t1_x_t2);
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

        REQUIRE(size(tn.tensors()) == 2);
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
        // std::wcout <<
        // to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) <<
        // std::endl; std::wcout <<
        // to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) <<
        // std::endl;
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                L"{\\tilde{a}^{{i_1}}_{{i_3}}}");
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                L"{F^{{i_{17}}}_{{i_1}}}");
      }

      // with explicit named indices
      SECTION("explicit") {
        Index::reset_tmp_index();
        TensorNetworkV2 tn(*t1_x_t2);

        using named_indices_t = TensorNetworkV2::NamedIndexSet;
        named_indices_t indices{Index{L"i_17"}};
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false,
                        &indices);

        REQUIRE(size(tn.tensors()) == 2);
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]));
        REQUIRE(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]));
        //        std::wcout <<
        //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0]))
        //        << std::endl; std::wcout <<
        //        to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1]))
        //        << std::endl;
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                L"{\\tilde{a}^{{i_2}}_{{i_1}}}");
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                L"{F^{{i_{17}}}_{{i_2}}}");
      }
    }

    SECTION("particle non-conserving") {
      const auto input1 = parse_expr(L"P{;a1,a3}");
      const auto input2 = parse_expr(L"P{a1,a3;}");
      const std::wstring expected1 = L"{{P^{{a_1}{a_3}}_{}}}";
      const std::wstring expected2 = L"{{P^{}_{{a_1}{a_3}}}}";

      for (int variant : {1, 2}) {
        for (bool fast : {true, false}) {
          TensorNetworkV2 tn(
              std::vector<ExprPtr>{variant == 1 ? input1 : input2});
          tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), fast);
          REQUIRE(tn.tensors().size() == 1);
          auto result = ex<Product>(to_tensors(tn.tensors()));
          REQUIRE(to_latex(result) == (variant == 1 ? expected1 : expected2));
        }
      }
    }

    SECTION("non-symmetric") {
      const auto input =
          parse_expr(L"A{i9,i12;i7,i3}:A I1{i7,i3;;x5}:N I2{;i9,i12;x5}:N")
              .as<Product>()
              .factors();
      const std::wstring expected =
          L"A{i_1,i_2;i_3,i_4}:A * I1{i_3,i_4;;x_1}:N * I2{;i_1,i_2;x_1}:N";

      for (bool fast : {true, false}) {
        TensorNetworkV2 tn(input);
        tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), fast);
        const auto result = ex<Product>(to_tensors(tn.tensors()));
        REQUIRE_THAT(result, SimplifiesTo(expected));
      }
    }

    SECTION("particle-1,2-symmetry") {
      const std::vector<std::pair<std::wstring, std::wstring>> pairs = {
          {L"S{i_1,i_2,i_3;a_1,a_2,a_3}:N * f{i_4;i_2}:N * "
           L"t{a_1,a_2,a_3;i_4,i_3,i_1}:N",
           L"S{i_1,i_2,i_3;a_1,a_2,a_3}:N * f{i_4;i_1}:N * "
           L"t{a_1,a_2,a_3;i_2,i_3,i_4}:N"},
          {L"Γ{o_2,o_4;o_1,o_3}:N * g{i_1,o_1;o_2,e_1}:N * "
           L"t{o_3,e_1;o_4,i_1}:N",
           L"Γ{o_2,o_4;o_1,o_3}:N * g{i_1,o_3;o_4,e_1}:N * "
           L"t{o_1,e_1;o_2,i_1}:N"}};
      for (const auto& pair : pairs) {
        const auto first = parse_expr(pair.first).as<Product>().factors();
        const auto second = parse_expr(pair.second).as<Product>().factors();

        TensorNetworkV2Accessor accessor;
        auto [first_graph, first_labels] =
            accessor.get_canonical_bliss_graph(TensorNetworkV2(first));
        auto [second_graph, second_labels] =
            accessor.get_canonical_bliss_graph(TensorNetworkV2(second));
        if (first_graph->cmp(*second_graph) != 0) {
          std::wstringstream stream;
          stream << "First graph:\n";
          first_graph->write_dot(stream, first_labels, {},
                                 {.display_colors = true});
          stream << "Second graph:\n";
          second_graph->write_dot(stream, second_labels, {},
                                  {.display_colors = true});
          stream << "TN graph:\n";
          auto [wick_graph, labels, texlabels, d1, d2] =
              TensorNetwork(first).make_bliss_graph();
          wick_graph->write_dot(stream, labels, texlabels,
                                {.display_colors = true});

          FAIL(to_string(stream.str()));
        }

        TensorNetworkV2 tn1(first);
        TensorNetworkV2 tn2(second);

        tn1.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);
        tn2.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

        REQUIRE(tn1.tensors().size() == tn2.tensors().size());
        for (std::size_t i = 0; i < tn1.tensors().size(); ++i) {
          auto t1 = std::dynamic_pointer_cast<Expr>(tn1.tensors()[i]);
          auto t2 = std::dynamic_pointer_cast<Expr>(tn2.tensors()[i]);
          REQUIRE(t1);
          REQUIRE(t2);
          REQUIRE(to_latex(t1) == to_latex(t2));
        }
      }
    }

    SECTION("miscellaneous") {
      const std::vector<std::pair<std::wstring, std::wstring>> inputs = {
          {L"g{i_1,a_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A",
           L"g{i_1,a_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A"},
          {L"g{a_1,i_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A",
           L"-1 g{i_1,a_1;i_2,i_3}:A * I{i_2,i_3;i_1,a_1}:A"},

          {L"g{i_1,a_1;i_2,i_3}:N * I{i_2,i_3;i_1,a_1}:N",
           L"g{i_1,a_1;i_2,i_3}:N * I{i_2,i_3;i_1,a_1}:N"},
          {L"g{a_1,i_1;i_2,i_3}:N * I{i_2,i_3;i_1,a_1}:N",
           L"g{i_1,a_1;i_2,i_3}:N * I{i_3,i_2;i_1,a_1}:N"},
      };

      for (const auto& [input, expected] : inputs) {
        const auto input_tensors = parse_expr(input).as<Product>().factors();

        TensorNetworkV2 tn(input_tensors);
        ExprPtr factor = tn.canonicalize(
            TensorCanonicalizer::cardinal_tensor_labels(), true);

        ExprPtr prod = to_product(tn.tensors());
        if (factor) {
          prod = ex<Product>(
              prod.as<Product>().scale(factor.as<Constant>().value()));
        }

        REQUIRE_THAT(prod, SimplifiesTo(expected));
      }
    }

    SECTION("special") {
      auto factors =
          parse_expr(
              L"S{i_1;a_1<i_1>}:N-C-S g{i_2,a_1<i_1>;a_2<i_2>,i_1}:N-C-S "
              L"t{a_2<i_2>;i_2}:N-C-S")
              ->as<Product>()
              .factors();

      TensorNetworkV2 tn(factors);

      ExprPtr factor =
          tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);
      ExprPtr result = to_product(tn.tensors());
      if (factor) {
        result *= factor;
      }

      REQUIRE(result);
    }

#ifndef SEQUANT_SKIP_LONG_TESTS
    SECTION("Exhaustive SRCC example") {
      // Note: the exact canonical form written here is implementation-defined
      // and doesn't actually matter What does, is that all equivalent ways of
      // writing it down, canonicalizes to the same exact form
      const Product expectedExpr =
          parse_expr(
              L"A{i1,i2;a1,a2} g{i3,i4;a3,a4} t{a1,a3;i3,i4} t{a2,a4;i1,i2}",
              Symmetry::antisymm)
              .as<Product>();

      const auto expected = expectedExpr.factors();

      TensorNetworkV2Accessor accessor;
      const auto [canonical_graph, canonical_graph_labels] =
          accessor.get_canonical_bliss_graph(TensorNetworkV2(expected));

      //      std::wcout << "Canonical graph:\n";
      //      canonical_graph->write_dot(std::wcout, canonical_graph_labels);
      //      std::wcout << std::endl;

      std::vector<Index> indices;
      for (std::size_t i = 0; i < expected.size(); ++i) {
        const Tensor& tensor = expected[i].as<Tensor>();
        for (const Index& idx : tensor.indices()) {
          if (std::find(indices.begin(), indices.end(), idx) == indices.end()) {
            indices.push_back(idx);
          }
        }
      }
      std::sort(indices.begin(), indices.end());

      const auto original_indices = indices;

      // Make sure to clone all expressions in order to not accidentally
      // modify the ones in expected (even though they are const... the
      // pointer-like semantics of expressions messes with const semantics)
      std::remove_const_t<decltype(expected)> factors;
      for (const auto& factor : expected) {
        factors.push_back(factor.clone());
      }
      std::sort(factors.begin(), factors.end());

      const auto is_occ = [](const Index& idx) {
        return idx.space() == Index(L"i_1").space();
      };

      // Iterate over all tensor permutations and all permutations of possible
      // index name swaps
      REQUIRE(std::is_sorted(factors.begin(), factors.end()));
      REQUIRE(std::is_sorted(indices.begin(), indices.end()));
      REQUIRE(std::is_partitioned(indices.begin(), indices.end(), is_occ));
      REQUIRE(std::partition_point(indices.begin(), indices.end(), is_occ) ==
              indices.begin() + 4);
      std::size_t total_variations = 0;
      do {
        do {
          do {
            total_variations++;

            // Compute index replacements
            container::map<Index, Index> idxrepl;
            for (std::size_t i = 0; i < indices.size(); ++i) {
              REQUIRE(original_indices[i].space() == indices[i].space());

              idxrepl.insert(
                  std::make_pair(original_indices.at(i), indices.at(i)));
            }

            // Apply index replacements to a copy of the current tensor
            // permutation
            auto copy = factors;
            for (ExprPtr& expr : copy) {
              expr.as<Tensor>().transform_indices(idxrepl);
              reset_tags(expr.as<Tensor>());
            }

            TensorNetworkV2 tn(copy);

            // At the heart of our canonicalization lies the fact that we can
            // always create the uniquely defined canonical graph for a given
            // network
            const auto [current_graph, current_graph_labels] =
                accessor.get_canonical_bliss_graph(tn);
            if (current_graph->cmp(*canonical_graph) != 0) {
              std::wcout << "Canonical graph for " << deparse(ex<Product>(copy))
                         << ":\n";
              current_graph->write_dot(std::wcout, current_graph_labels);
              std::wcout << std::endl;
            }
            REQUIRE(current_graph->cmp(*canonical_graph) == 0);

            tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                            false);

            std::vector<ExprPtr> actual;
            std::transform(tn.tensors().begin(), tn.tensors().end(),
                           std::back_inserter(actual), [](const auto& t) {
                             assert(std::dynamic_pointer_cast<Expr>(t));
                             return std::dynamic_pointer_cast<Expr>(t);
                           });

            // The canonical graph must not change due to the other
            // canonicalization steps we perform
            REQUIRE(accessor.get_canonical_bliss_graph(TensorNetworkV2(actual))
                        .first->cmp(*canonical_graph) == 0);

            REQUIRE(actual.size() == expected.size());

            if (!std::equal(expected.begin(), expected.end(), actual.begin())) {
              std::wostringstream sstream;
              sstream
                  << "Expected all tensors to be equal (actual == expected), "
                     "but got:\n";
              for (std::size_t i = 0; i < expected.size(); ++i) {
                std::wstring equality =
                    actual[i] == expected[i] ? L" == " : L" != ";

                sstream << deparse(actual[i]) << equality
                        << deparse(expected[i]) << "\n";
              }
              sstream << "\nInput was " << deparse(ex<Product>(factors))
                      << "\n";
              FAIL(to_string(sstream.str()));
            }
          } while (std::next_permutation(indices.begin() + 4, indices.end()));
        } while (std::next_permutation(indices.begin(), indices.begin() + 4));
      } while (std::next_permutation(factors.begin(), factors.end()));

      // 4! (tensors) * 4! (internal indices) * 4! (external indices)
      REQUIRE(total_variations == 24 * 24 * 24);
    }
#endif

    SECTION("idempotency") {
      const std::vector<std::wstring> inputs = {
          L"F{i1;i8} g{i8,i9;i1,i7}",
          L"A{i9,i12;i7,i3}:A I1{i7,i3;;x5}:N I2{;i9,i12;x5}:N",
          L"f{i4;i1}:N t{a1,a2,a3;i2,i3,i4}:N S{i1,i2,i3;a1,a2,a3}:N",
          L"P{a1,a3;} k{i8;i2}",
          L"L{x6;;x2} P{;a1,a3}",
      };

      for (const std::wstring& current : inputs) {
        auto factors1 = parse_expr(current).as<Product>().factors();
        auto factors2 = parse_expr(current).as<Product>().factors();

        TensorNetworkV2 reference_tn(factors1);
        reference_tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                                  false);

        TensorNetworkV2 check_tn(factors2);
        check_tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                              false);

        REQUIRE(to_latex(to_product(reference_tn.tensors())) ==
                to_latex(to_product(check_tn.tensors())));

        for (bool fast : {true, false, true, true, false, false, true}) {
          reference_tn.canonicalize(
              TensorCanonicalizer::cardinal_tensor_labels(), fast);

          REQUIRE(to_latex(to_product(reference_tn.tensors())) ==
                  to_latex(to_product(check_tn.tensors())));
        }
      }
    }  // SECTION("idempotency")

  }  // SECTION("canonicalizer")

  SECTION("misc1") {
    if (false) {
      Index::reset_tmp_index();
      // TN1 from manuscript
      auto g = ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                          Symmetry::antisymm);
      auto ta = ex<Tensor>(L"t", bra{L"a_1", L"a_3"}, ket{L"i_1", L"i_2"},
                           Symmetry::antisymm);
      auto tb = ex<Tensor>(L"t", bra{L"a_2", L"a_4"}, ket{L"i_3", L"i_4"},
                           Symmetry::antisymm);

      auto tmp = g * ta * tb;
      // std::wcout << "TN1 = " << to_latex(tmp) << std::endl;
      TensorNetworkV2 tn(tmp->as<Product>().factors());

      // make graph
      // N.B. treat all indices as dummy so that the automorphism ignores the
      using named_indices_t = TensorNetworkV2::NamedIndexSet;
      named_indices_t indices{};
      REQUIRE_NOTHROW(tn.create_graph({.named_indices = &indices}));
      TensorNetworkV2::Graph graph =
          tn.create_graph({.named_indices = &indices});

      // can disable label production
      {
        // make sure can generate without labels also
        REQUIRE_NOTHROW(
            tn.create_graph({.make_labels = true, .make_texlabels = false})
                .vertex_texlabels.size() == 0);
        REQUIRE_NOTHROW(
            tn.create_graph({.make_labels = false, .make_texlabels = true})
                .vertex_labels.size() == 0);
        {
          auto g =
              tn.create_graph({.make_labels = false, .make_texlabels = false});
          REQUIRE_NOTHROW(g.vertex_labels.size() == 0 &&
                          g.vertex_texlabels.size() == 0);
        }
      }

      // create dot
      {
        std::basic_ostringstream<wchar_t> oss;
        REQUIRE_NOTHROW(graph.bliss_graph->write_dot(oss, graph.vertex_labels));
        // std::wcout << "oss.str() = " << std::endl << oss.str() <<
        // std::endl;
      }

      bliss::Stats stats;
      graph.bliss_graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

      std::vector<std::vector<unsigned int>> aut_generators;
      auto save_aut = [&aut_generators](const unsigned int n,
                                        const unsigned int* aut) {
        aut_generators.emplace_back(aut, aut + n);
      };
      graph.bliss_graph->find_automorphisms(
          stats, &bliss::aut_hook<decltype(save_aut)>, &save_aut);
      CHECK(aut_generators.size() ==
            2);  // there are 2 generators, i1<->i2, i3<->i4

      std::basic_ostringstream<wchar_t> oss;
      bliss::print_auts(aut_generators, oss, graph.vertex_labels);
      CHECK(oss.str() == L"({i_3},{i_4})\n({i_1},{i_2})\n");
      // std::wcout << oss.str() << std::endl;
    }
  }
}
