//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/mbpt/sr.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <range/v3/all.hpp>

// TODO: Add test cases with auxiliary indices

TEST_CASE("TensorNetwork", "[elements]") {
  using namespace sequant;

  using namespace sequant::mbpt::sr;

  SECTION("constructors") {
    {  // with Tensors
      auto t1 =
          ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t2 =
          ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));

      auto t1_x_t2_p_t2 = t1 * (t2 + t2);  // can only use a flat tensor product
      REQUIRE_THROWS_AS(TensorNetwork(*t1_x_t2_p_t2), std::logic_error);
    }

    {  // with NormalOperators
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_2"}, V);
      auto t2 = ex<FNOperator>(WstrList{L"i_2"}, WstrList{L"i_1"}, V);
      auto t1_x_t2 = t1 * t2;
      REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));
    }

    {  // with Tensors and NormalOperators
      auto tmp = A(-2) * H_(2) * T_(2) * T_(2);
      REQUIRE_NOTHROW(TensorNetwork(tmp->as<Product>().factors()));
    }

  }  // SECTION("constructors")

  SECTION("accessors") {
    {
      constexpr const auto V = Vacuum::SingleProduct;
      auto t1 =
          ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
      auto t2 = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_3"}, V);
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
        auto t1 =
            ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}, WstrList{});
        auto t2 = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_2"}, V);
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
                L"{F^{{i_1}}_{{i_2}}}");
        REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[1])) ==
                L"{\\tilde{a}^{{i_2}}_{{i_1}}}");
        REQUIRE(tn.idxrepl().size() == 2);
      }

      {
        Index::reset_tmp_index();
        constexpr const auto V = Vacuum::SingleProduct;
        auto t1 =
            ex<Tensor>(L"F", WstrList{L"i_2"}, WstrList{L"i_17"}, WstrList{});
        auto t2 = ex<FNOperator>(WstrList{L"i_2"}, WstrList{L"i_3"}, V);
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
                  L"{\\tilde{a}^{{i_2}}_{{i_1}}}");
          REQUIRE(to_latex(std::dynamic_pointer_cast<Expr>(tn.tensors()[0])) ==
                  L"{F^{{i_{17}}}_{{i_2}}}");
        }
      }
    }
  }  // SECTION("canonicalizer")

  SECTION("bliss graph") {
    Index::reset_tmp_index();
    // to generate expressions in specified (i.e., platform-independent) manner
    // can't use operator expression (due to unspecified order of evaluation of
    // function arguments), must use initializer list
    auto tmp = ex<Product, std::initializer_list<ExprPtr>>(
        {A(-2), H_(2), T_(2), T_(2), T_(2)});
    // canonicalize to avoid dependence on the implementation details of
    // mbpt::sr::make_op
    canonicalize(tmp);
    // std::wcout << "A2*H2*T2*T2*T2 = " << to_latex(tmp) << std::endl;
    TensorNetwork tn(tmp->as<Product>().factors());

    // make graph
    REQUIRE_NOTHROW(tn.make_bliss_graph());
    auto [graph, vlabels, vcolors, vtypes] = tn.make_bliss_graph();

    // create dot
    std::basic_ostringstream<wchar_t> oss;
    REQUIRE_NOTHROW(graph->write_dot(oss, vlabels));
    std::wcout << "oss.str() = " << std::endl << oss.str() << std::endl;
    REQUIRE(oss.str() ==
            L"graph g {\n"
            "v0 [label=\"{a_1}\"; color=\"#9e3,ba0\"];\n"
            "v0 -- v29\n"
            "v0 -- v58\n"
            "v1 [label=\"{a_2}\"; color=\"#9e3,ba0\"];\n"
            "v1 -- v29\n"
            "v1 -- v58\n"
            "v2 [label=\"{a_3}\"; color=\"#9e3,ba0\"];\n"
            "v2 -- v33\n"
            "v2 -- v54\n"
            "v3 [label=\"{a_4}\"; color=\"#9e3,ba0\"];\n"
            "v3 -- v33\n"
            "v3 -- v54\n"
            "v4 [label=\"{a_5}\"; color=\"#9e3,ba0\"];\n"
            "v4 -- v37\n"
            "v4 -- v50\n"
            "v5 [label=\"{a_6}\"; color=\"#9e3,ba0\"];\n"
            "v5 -- v37\n"
            "v5 -- v50\n"
            "v6 [label=\"{a_7}\"; color=\"#9e3,ba0\"];\n"
            "v6 -- v22\n"
            "v6 -- v41\n"
            "v7 [label=\"{a_8}\"; color=\"#9e3,ba0\"];\n"
            "v7 -- v22\n"
            "v7 -- v41\n"
            "v8 [label=\"{i_1}\"; color=\"#a78,ee8\"];\n"
            "v8 -- v30\n"
            "v8 -- v57\n"
            "v9 [label=\"{i_2}\"; color=\"#a78,ee8\"];\n"
            "v9 -- v30\n"
            "v9 -- v57\n"
            "v10 [label=\"{i_3}\"; color=\"#a78,ee8\"];\n"
            "v10 -- v34\n"
            "v10 -- v53\n"
            "v11 [label=\"{i_4}\"; color=\"#a78,ee8\"];\n"
            "v11 -- v34\n"
            "v11 -- v53\n"
            "v12 [label=\"{i_5}\"; color=\"#a78,ee8\"];\n"
            "v12 -- v38\n"
            "v12 -- v49\n"
            "v13 [label=\"{i_6}\"; color=\"#a78,ee8\"];\n"
            "v13 -- v38\n"
            "v13 -- v49\n"
            "v14 [label=\"{i_7}\"; color=\"#a78,ee8\"];\n"
            "v14 -- v21\n"
            "v14 -- v42\n"
            "v15 [label=\"{i_8}\"; color=\"#a78,ee8\"];\n"
            "v15 -- v21\n"
            "v15 -- v42\n"
            "v16 [label=\"{\\kappa_1}\"; color=\"#703,062\"];\n"
            "v16 -- v25\n"
            "v16 -- v46\n"
            "v17 [label=\"{\\kappa_2}\"; color=\"#703,062\"];\n"
            "v17 -- v25\n"
            "v17 -- v46\n"
            "v18 [label=\"{\\kappa_3}\"; color=\"#703,062\"];\n"
            "v18 -- v26\n"
            "v18 -- v45\n"
            "v19 [label=\"{\\kappa_4}\"; color=\"#703,062\"];\n"
            "v19 -- v26\n"
            "v19 -- v45\n"
            "v20 [label=\"A\"; color=\"#518,020\"];\n"
            "v20 -- v23\n"
            "v21 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v21 -- v23\n"
            "v22 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v22 -- v23\n"
            "v23 [label=\"bka\"; color=\"#518,020\"];\n"
            "v24 [label=\"g\"; color=\"#2e0,351\"];\n"
            "v24 -- v27\n"
            "v25 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v25 -- v27\n"
            "v26 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v26 -- v27\n"
            "v27 [label=\"bka\"; color=\"#2e0,351\"];\n"
            "v28 [label=\"t\"; color=\"#43,e44\"];\n"
            "v28 -- v31\n"
            "v29 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v29 -- v31\n"
            "v30 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v30 -- v31\n"
            "v31 [label=\"bka\"; color=\"#43,e44\"];\n"
            "v32 [label=\"t\"; color=\"#43,e44\"];\n"
            "v32 -- v35\n"
            "v33 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v33 -- v35\n"
            "v34 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v34 -- v35\n"
            "v35 [label=\"bka\"; color=\"#43,e44\"];\n"
            "v36 [label=\"t\"; color=\"#43,e44\"];\n"
            "v36 -- v39\n"
            "v37 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v37 -- v39\n"
            "v38 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v38 -- v39\n"
            "v39 [label=\"bka\"; color=\"#43,e44\"];\n"
            "v40 [label=\"ã\"; color=\"#cbf,be5\"];\n"
            "v40 -- v43\n"
            "v41 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v41 -- v43\n"
            "v42 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v42 -- v43\n"
            "v43 [label=\"bka\"; color=\"#cbf,be5\"];\n"
            "v44 [label=\"ã\"; color=\"#cbf,be5\"];\n"
            "v44 -- v47\n"
            "v45 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v45 -- v47\n"
            "v46 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v46 -- v47\n"
            "v47 [label=\"bka\"; color=\"#cbf,be5\"];\n"
            "v48 [label=\"ã\"; color=\"#cbf,be5\"];\n"
            "v48 -- v51\n"
            "v49 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v49 -- v51\n"
            "v50 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v50 -- v51\n"
            "v51 [label=\"bka\"; color=\"#cbf,be5\"];\n"
            "v52 [label=\"ã\"; color=\"#cbf,be5\"];\n"
            "v52 -- v55\n"
            "v53 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v53 -- v55\n"
            "v54 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v54 -- v55\n"
            "v55 [label=\"bka\"; color=\"#cbf,be5\"];\n"
            "v56 [label=\"ã\"; color=\"#cbf,be5\"];\n"
            "v56 -- v59\n"
            "v57 [label=\"bra2a\"; color=\"#eaa,2ab\"];\n"
            "v57 -- v59\n"
            "v58 [label=\"ket2a\"; color=\"#5a8,fd3\"];\n"
            "v58 -- v59\n"
            "v59 [label=\"bka\"; color=\"#cbf,be5\"];\n"
            "}\n");

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
      REQUIRE(oss.str() ==
              L"(18,19)\n"
              "(16,17)\n"
              "(6,7)\n"
              "(14,15)\n"
              "(0,1)\n"
              "(8,9)\n"
              "(2,3)\n"
              "(4,5)\n"
              "(10,11)\n"
              "(12,13)\n"
              "(2,4)(3,5)(10,12)(11,13)(32,36)(33,37)(34,38)(35,39)(48,52)(49,"
              "53)(50,54)(51,55)\n"
              "(0,2)(1,3)(8,10)(9,11)(28,32)(29,33)(30,34)(31,35)(52,56)(53,57)"
              "(54,58)(55,59)\n");
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
      auto g =
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},
                     WstrList{}, Symmetry::antisymm);
      auto ta =
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_3"}, WstrList{L"i_1", L"i_2"},
                     WstrList{}, Symmetry::antisymm);
      auto tb =
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_4"}, WstrList{L"i_3", L"i_4"},
                     WstrList{}, Symmetry::antisymm);

      auto tmp = g * ta * tb;
      // std::wcout << "TN1 = " << to_latex(tmp) << std::endl;
      TensorNetwork tn(tmp->as<Product>().factors());

      // make graph
      // N.B. treat all indices as dummy so that the automorphism ignores the
      using named_indices_t = TensorNetwork::named_indices_t;
      named_indices_t indices{};
      REQUIRE_NOTHROW(tn.make_bliss_graph(&indices));
      auto [graph, vlabels, vcolors, vtypes] = tn.make_bliss_graph(&indices);

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
                      L"u", idxs, std::vector<Index>{}, std::vector<Index>{},
                      (testcase == 3
                           ? Symmetry::nonsymm
                           : ((n == 1) ? Symmetry::symm : Symmetry::nonsymm)));
                }) |
                ranges::to_vector;
            CHECK(utensors.size() == static_cast<std::size_t>(n));
            auto dtensors = contravariant_indices | ranges::views::chunk(N) |
                            ranges::views::transform([&](const auto& idxs) {
                              return ex<Tensor>(L"d", std::vector<Index>{}, std::vector<Index>{},
                                                idxs, Symmetry::nonsymm);
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
              auto [graph, vlabels, vcolors, vtypes] = tn.make_bliss_graph();

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

}  // TEST_CASE("Tensor")
