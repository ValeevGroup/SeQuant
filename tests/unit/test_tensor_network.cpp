//
// Created by Eduard Valeyev on 3/23/18.
//

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

TEST_CASE("TensorNetwork", "[elements]") {
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

    // make graph
    REQUIRE_NOTHROW(tn.make_bliss_graph());
    auto [graph, vlabels, vtexlabels, vcolors, vtypes] = tn.make_bliss_graph();

    // create dot
    std::basic_ostringstream<wchar_t> oss;
    REQUIRE_NOTHROW(graph->write_dot(oss, vlabels, vtexlabels));
    // std::wcout << "oss.str() = " << std::endl << oss.str() << std::endl;
    const std::wstring actual = oss.str();
    // clang-format off
    const std::wstring expected =
        L"graph g {\n"
"node [ style=filled, penwidth=2, margin=0];\n"
"v0 [ label=\"a_1\", texlbl=\"${a_1}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v0 -- v29\n"
"v0 -- v58\n"
"v1 [ label=\"a_2\", texlbl=\"${a_2}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v1 -- v29\n"
"v1 -- v58\n"
"v2 [ label=\"a_3\", texlbl=\"${a_3}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v2 -- v33\n"
"v2 -- v54\n"
"v3 [ label=\"a_4\", texlbl=\"${a_4}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v3 -- v33\n"
"v3 -- v54\n"
"v4 [ label=\"a_5\", texlbl=\"${a_5}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v4 -- v37\n"
"v4 -- v50\n"
"v5 [ label=\"a_6\", texlbl=\"${a_6}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v5 -- v37\n"
"v5 -- v50\n"
"v6 [ label=\"a_7\", texlbl=\"${a_7}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v6 -- v22\n"
"v6 -- v41\n"
"v7 [ label=\"a_8\", texlbl=\"${a_8}$\", color=\"#13db00\", fillcolor=\"#98db04\" ];\n"
"v7 -- v22\n"
"v7 -- v41\n"
"v8 [ label=\"i_1\", texlbl=\"${i_1}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v8 -- v30\n"
"v8 -- v57\n"
"v9 [ label=\"i_2\", texlbl=\"${i_2}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v9 -- v30\n"
"v9 -- v57\n"
"v10 [ label=\"i_3\", texlbl=\"${i_3}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v10 -- v34\n"
"v10 -- v53\n"
"v11 [ label=\"i_4\", texlbl=\"${i_4}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v11 -- v34\n"
"v11 -- v53\n"
"v12 [ label=\"i_5\", texlbl=\"${i_5}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v12 -- v38\n"
"v12 -- v49\n"
"v13 [ label=\"i_6\", texlbl=\"${i_6}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v13 -- v38\n"
"v13 -- v49\n"
"v14 [ label=\"i_7\", texlbl=\"${i_7}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v14 -- v21\n"
"v14 -- v42\n"
"v15 [ label=\"i_8\", texlbl=\"${i_8}$\", color=\"#3e2d55\", fillcolor=\"#f8b4aa\" ];\n"
"v15 -- v21\n"
"v15 -- v42\n"
"v16 [ label=\"κ_1\", texlbl=\"${\\kappa_1}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
"v16 -- v25\n"
"v16 -- v46\n"
"v17 [ label=\"κ_2\", texlbl=\"${\\kappa_2}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
"v17 -- v25\n"
"v17 -- v46\n"
"v18 [ label=\"κ_3\", texlbl=\"${\\kappa_3}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
"v18 -- v26\n"
"v18 -- v45\n"
"v19 [ label=\"κ_4\", texlbl=\"${\\kappa_4}$\", color=\"#e174c5\", fillcolor=\"#e1e8c5\" ];\n"
"v19 -- v26\n"
"v19 -- v45\n"
"v20 [ label=\"A\", color=\"#257a61\", fillcolor=\"#94f4c2\" ];\n"
"v20 -- v23\n"
"v21 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v21 -- v23\n"
"v22 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v22 -- v23\n"
"v23 [ label=\"bka\", color=\"#257a61\", fillcolor=\"#94f4c2\" ];\n"
"v24 [ label=\"g\", color=\"#300a49\", fillcolor=\"#c05092\" ];\n"
"v24 -- v27\n"
"v25 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v25 -- v27\n"
"v26 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v26 -- v27\n"
"v27 [ label=\"bka\", color=\"#300a49\", fillcolor=\"#c05092\" ];\n"
"v28 [ label=\"t\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
"v28 -- v31\n"
"v29 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v29 -- v31\n"
"v30 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v30 -- v31\n"
"v31 [ label=\"bka\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
"v32 [ label=\"t\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
"v32 -- v35\n"
"v33 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v33 -- v35\n"
"v34 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v34 -- v35\n"
"v35 [ label=\"bka\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
"v36 [ label=\"t\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
"v36 -- v39\n"
"v37 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v37 -- v39\n"
"v38 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v38 -- v39\n"
"v39 [ label=\"bka\", color=\"#e812d9\", fillcolor=\"#e890d9\" ];\n"
"v40 [ label=\"ã\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v40 -- v43\n"
"v41 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v41 -- v43\n"
"v42 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v42 -- v43\n"
"v43 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v44 [ label=\"ã\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v44 -- v47\n"
"v45 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v45 -- v47\n"
"v46 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v46 -- v47\n"
"v47 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v48 [ label=\"ã\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v48 -- v51\n"
"v49 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v49 -- v51\n"
"v50 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v50 -- v51\n"
"v51 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v52 [ label=\"ã\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v52 -- v55\n"
"v53 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v53 -- v55\n"
"v54 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v54 -- v55\n"
"v55 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v56 [ label=\"ã\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"v56 -- v59\n"
"v57 [ label=\"bra2a\", color=\"#06d223\", fillcolor=\"#30d28c\" ];\n"
"v57 -- v59\n"
"v58 [ label=\"ket2a\", color=\"#c849cb\", fillcolor=\"#c892cb\" ];\n"
"v58 -- v59\n"
"v59 [ label=\"bka\", color=\"#2a13ee\", fillcolor=\"#a898ee\" ];\n"
"}\n";
    // clang-format on

    REQUIRE(actual == expected);

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
          L"(14,15)\n"
          L"(6,7)\n"
          L"(12,13)\n"
          L"(4,5)\n"
          L"(10,11)\n"
          L"(8,9)\n"
          L"(2,3)\n"
          L"(0,1)\n"
          L"(0,2)(1,3)(8,10)(9,11)(28,32)(29,33)(30,34)(31,35)(52,56)(53,57)("
          L"54,58)(55,59)\n"
          L"(2,4)(3,5)(10,12)(11,13)(32,36)(33,37)(34,38)(35,39)(48,52)(49,53)("
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
      REQUIRE_NOTHROW(tn.make_bliss_graph(&indices));
      auto [graph, vlabels, vtexlabels, vcolors, vtypes] =
          tn.make_bliss_graph(&indices);

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

    // profile canonicalization for synthetic tests in
    // DOI 10.1016/j.cpc.2018.02.014
    if (false) {
      for (auto testcase : {0, 1, 2, 3}) {
        // - testcase=0,2 are "equivalent" and correspond to the "frustrated"
        //   case in Section 5.3 of DOI 10.1016/j.cpc.2018.02.014
        // - testcase=1 corresponds to the "frustrated" case in Section 5.4 of
        //   DOI 10.1016/j.cpc.2018.02.014
        // - testcase=3 corresponds to the "No symmetry dummy"
        //   case in Section 5.1 of DOI 10.1016/j.cpc.2018.02.014
        if (testcase == 0)
          std::wcout
              << "canonicalizing network with 1 totally-symmetric tensor with "
                 "N indices and 1 asymmetric tensor with N indices\n";
        else if (testcase == 3)
          std::wcout << "canonicalizing network with 1 asymmetric tensor with "
                        "N indices and 1 asymmetric tensor with N indices\n";
        else
          std::wcout << "canonicalizing network with n equivalent asymmetric "
                        "tensors with N/n indices each and 1 asymmetric tensor "
                        "with N indices\n";

        std::wcout << "N,n,min_time,geommean_time,max_time\n";

        for (auto N :
             {1, 2, 4, 8, 16, 32, 64, 128, 256}) {  // total number of indices

          int n;
          switch (testcase) {
            case 0:
              n = 1;
              break;
            case 1:
              n = N / 2;
              break;
            case 2:
              n = N;
              break;
            case 3:
              n = 1;
              break;
            default:
              abort();
          }
          if (n == 0 || n > N) continue;

          auto ctx_resetter = set_scoped_default_context(
              (static_cast<std::size_t>(N) > Index::min_tmp_index())
                  ? Context(get_default_context())
                        .set_first_dummy_index_ordinal(N + 1)
                  : get_default_context());

          // make list of indices
          std::vector<Index> indices;
          for (auto i = 0; i != N; ++i) {
            std::wostringstream oss;
            oss << "i_" << i;
            indices.emplace_back(oss.str());
          }
          std::random_device rd;

          // randomly sample connectivity between bra and ket tensors
          const auto S = 10;  // how many samples to take

          auto product_time =
              1.;  // product of all times, need to get geometric mean
          auto min_time =
              std::numeric_limits<double>::max();  // total time for all samples
          auto max_time =
              std::numeric_limits<double>::min();  // total time for all samples
          for (auto s = 0; s != S; ++s) {
            // make tensors of independently (and randomly) permuted
            // contravariant and covariant indices
            auto contravariant_indices = indices;
            auto covariant_indices = indices;

            std::shuffle(contravariant_indices.begin(),
                         contravariant_indices.end(), std::mt19937{rd()});
            std::shuffle(covariant_indices.begin(), covariant_indices.end(),
                         std::mt19937{rd()});

            auto utensors =
                covariant_indices | ranges::views::chunk(N / n) |
                ranges::views::transform([&](const auto& idxs) {
                  return ex<Tensor>(
                      L"u", bra(idxs), ket{},
                      (testcase == 3
                           ? Symmetry::nonsymm
                           : ((n == 1) ? Symmetry::symm : Symmetry::nonsymm)));
                }) |
                ranges::to_vector;
            CHECK(utensors.size() == static_cast<std::size_t>(n));
            auto dtensors =
                contravariant_indices | ranges::views::chunk(N) |
                ranges::views::transform([&](const auto& idxs) {
                  return ex<Tensor>(L"d", bra{}, ket(idxs), Symmetry::nonsymm);
                }) |
                ranges::to_vector;
            CHECK(dtensors.size() == 1);

            ExprPtr expr;
            for (auto g = 0; g != n; ++g) {
              if (g == 0)
                expr = utensors[0] * dtensors[0];
              else
                expr = expr * utensors[g];
            }

            TensorNetwork tn(expr->as<Product>().factors());

            // produce misc data for publication
            if (false && s == 0) {
              std::wcout << "N=" << N << " n=" << n << " expr:\n"
                         << expr->to_latex() << std::endl;

              // make graph
              REQUIRE_NOTHROW(tn.make_bliss_graph());
              auto [graph, vlabels, vtexlabels, vcolors, vtypes] =
                  tn.make_bliss_graph();

              // create dot
              std::basic_ostringstream<wchar_t> oss;
              REQUIRE_NOTHROW(graph->write_dot(oss, vlabels));
              std::wcout << "bliss graph:" << std::endl
                         << oss.str() << std::endl;
            }

            sequant::TimerPool<> timer;
            timer.start();
            tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                            false);
            timer.stop();
            const auto elapsed_seconds = timer.read();
            product_time *= elapsed_seconds;
            min_time = std::min(min_time, elapsed_seconds);
            max_time = std::max(max_time, elapsed_seconds);
          }

          const auto geommean_time = std::pow(product_time, 1. / S);
          std::wcout << N << "," << n << "," << min_time << "," << geommean_time
                     << "," << max_time << "\n";
        }
      }
    }

  }  // SECTION("misc1")

}  // TEST_CASE("TensorNetwork")

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
    auto graph = tn.create_graph(named_indices);
    return std::make_pair(std::move(graph.bliss_graph), graph.vertex_labels);
  }
};
}  // namespace sequant

TEST_CASE("TensorNetworkV2", "[elements]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using sequant::Context;
  namespace t = sequant::mbpt::tensor;
  namespace o = sequant::mbpt::op;

  sequant::set_default_context(Context(
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
           L"g{i_1,a_1;i_2,i_3}:N * I{i_2,i_3;a_1,i_1}:N"},
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
      REQUIRE_NOTHROW(tn.create_graph(&indices));
      TensorNetworkV2::Graph graph = tn.create_graph(&indices);

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

    // profile canonicalization for synthetic tests in
    // DOI 10.1016/j.cpc.2018.02.014
    if (false) {
      for (auto testcase : {0, 1, 2, 3}) {
        // - testcase=0,2 are "equivalent" and correspond to the "frustrated"
        //   case in Section 5.3 of DOI 10.1016/j.cpc.2018.02.014
        // - testcase=1 corresponds to the "frustrated" case in Section 5.4 of
        //   DOI 10.1016/j.cpc.2018.02.014
        // - testcase=3 corresponds to the "No symmetry dummy"
        //   case in Section 5.1 of DOI 10.1016/j.cpc.2018.02.014
        if (testcase == 0)
          std::wcout << "canonicalizing network with 1 totally-symmetric "
                        "tensor with "
                        "N indices and 1 asymmetric tensor with N indices\n";
        else if (testcase == 3)
          std::wcout << "canonicalizing network with 1 asymmetric tensor with "
                        "N indices and 1 asymmetric tensor with N indices\n";
        else
          std::wcout << "canonicalizing network with n equivalent asymmetric "
                        "tensors with N/n indices each and 1 asymmetric tensor "
                        "with N indices\n";

        std::wcout << "N,n,min_time,geommean_time,max_time\n";

        for (auto N :
             {1, 2, 4, 8, 16, 32, 64, 128, 256}) {  // total number of indices

          int n;
          switch (testcase) {
            case 0:
              n = 1;
              break;
            case 1:
              n = N / 2;
              break;
            case 2:
              n = N;
              break;
            case 3:
              n = 1;
              break;
            default:
              abort();
          }
          if (n == 0 || n > N) continue;

          auto ctx_resetter = set_scoped_default_context(
              (static_cast<std::size_t>(N) > Index::min_tmp_index())
                  ? Context(get_default_context())
                        .set_first_dummy_index_ordinal(N + 1)
                  : get_default_context());

          // make list of indices
          std::vector<Index> indices;
          for (auto i = 0; i != N; ++i) {
            std::wostringstream oss;
            oss << "i_" << i;
            indices.emplace_back(oss.str());
          }
          std::random_device rd;

          // randomly sample connectivity between bra and ket tensors
          const auto S = 10;  // how many samples to take

          auto product_time =
              1.;  // product of all times, need to get geometric mean
          auto min_time = std::numeric_limits<double>::max();  // total time for
                                                               // all samples
          auto max_time = std::numeric_limits<double>::min();  // total time for
                                                               // all samples
          for (auto s = 0; s != S; ++s) {
            // make tensors of independently (and randomly) permuted
            // contravariant and covariant indices
            auto contravariant_indices = indices;
            auto covariant_indices = indices;

            std::shuffle(contravariant_indices.begin(),
                         contravariant_indices.end(), std::mt19937{rd()});
            std::shuffle(covariant_indices.begin(), covariant_indices.end(),
                         std::mt19937{rd()});

            auto utensors =
                covariant_indices | ranges::views::chunk(N / n) |
                ranges::views::transform([&](const auto& idxs) {
                  return ex<Tensor>(
                      L"u", bra(idxs), ket{},
                      (testcase == 3
                           ? Symmetry::nonsymm
                           : ((n == 1) ? Symmetry::symm : Symmetry::nonsymm)));
                }) |
                ranges::to_vector;
            CHECK(utensors.size() == static_cast<std::size_t>(n));
            auto dtensors =
                contravariant_indices | ranges::views::chunk(N) |
                ranges::views::transform([&](const auto& idxs) {
                  return ex<Tensor>(L"d", bra{}, ket(idxs), Symmetry::nonsymm);
                }) |
                ranges::to_vector;
            CHECK(dtensors.size() == 1);

            ExprPtr expr;
            for (auto g = 0; g != n; ++g) {
              if (g == 0)
                expr = utensors[0] * dtensors[0];
              else
                expr = expr * utensors[g];
            }

            TensorNetworkV2 tn(expr->as<Product>().factors());

            // produce misc data for publication
            if (false && s == 0) {
              std::wcout << "N=" << N << " n=" << n << " expr:\n"
                         << expr->to_latex() << std::endl;

              // make graph
              REQUIRE_NOTHROW(tn.create_graph());
              TensorNetworkV2::Graph graph = tn.create_graph();

              // create dot
              std::basic_ostringstream<wchar_t> oss;
              REQUIRE_NOTHROW(
                  graph.bliss_graph->write_dot(oss, graph.vertex_labels));
              std::wcout << "bliss graph:" << std::endl
                         << oss.str() << std::endl;
            }

            sequant::TimerPool<> timer;
            timer.start();
            tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                            false);
            timer.stop();
            const auto elapsed_seconds = timer.read();
            product_time *= elapsed_seconds;
            min_time = std::min(min_time, elapsed_seconds);
            max_time = std::max(max_time, elapsed_seconds);
          }

          const auto geommean_time = std::pow(product_time, 1. / S);
          std::wcout << N << "," << n << "," << min_time << "," << geommean_time
                     << "," << max_time << "\n";
        }
      }
    }

  }  // SECTION("misc1")

}  // TEST_CASE("TensorNetworkV2")
