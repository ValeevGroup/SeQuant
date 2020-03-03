//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>

#include "SeQuant/core/bliss.hpp"
#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor_network.hpp"
#include "SeQuant/domain/mbpt/sr/sr.hpp"

TEST_CASE("TensorNetwork", "[elements]") {

  using namespace sequant;

  using namespace sequant::mbpt::sr::so;

  SECTION("constructors") {
    {  // with Tensors
      auto t1 = ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"});
      auto t2 = ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"i_1"});
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

    { // with Tensors and NormalOperators
      auto tmp = A(2) * H2() * T_(2) * T_(2);
      REQUIRE_NOTHROW(TensorNetwork(tmp->as<Product>().factors()));
    };

  }  // SECTION("constructors")

  SECTION("bliss graph") {
    Index::reset_tmp_index();
    auto tmp = A(2) * H2() * T_(2) * T_(2) * T_(2);
    // std::wcout << "A2*H2*T2*T2*T2 = " << to_latex(tmp) << std::endl;
    TensorNetwork tn(tmp->as<Product>().factors());

    // make graph
    REQUIRE_NOTHROW(tn.make_bliss_graph());
    auto [graph, vlabels, vcolors, vtypes] = tn.make_bliss_graph();

    // create dot
    std::basic_ostringstream<wchar_t> oss;
    REQUIRE_NOTHROW(graph->write_dot(oss, vlabels));
    //std::wcout << oss.str() << std::endl;
    // the write_dot output seems to be platform-dependent, e.g. g++-8 on Linux
    // TODO investigate the source of platform dependence
//    REQUIRE(oss.str() ==
//            L"graph g {\n"
//            "v0 [label=\"{a_{102}}\"; color=\"#4f1dd0\"];\n"
//            "v0 -- v22\n"
//            "v0 -- v25\n"
//            "v1 [label=\"{a_{103}}\"; color=\"#4f1dd0\"];\n"
//            "v1 -- v22\n"
//            "v1 -- v25\n"
//            "v2 [label=\"{a_{110}}\"; color=\"#4f1dd0\"];\n"
//            "v2 -- v37\n"
//            "v2 -- v42\n"
//            "v3 [label=\"{a_{111}}\"; color=\"#4f1dd0\"];\n"
//            "v3 -- v37\n"
//            "v3 -- v42\n"
//            "v4 [label=\"{a_{114}}\"; color=\"#4f1dd0\"];\n"
//            "v4 -- v45\n"
//            "v4 -- v50\n"
//            "v5 [label=\"{a_{115}}\"; color=\"#4f1dd0\"];\n"
//            "v5 -- v45\n"
//            "v5 -- v50\n"
//            "v6 [label=\"{a_{118}}\"; color=\"#4f1dd0\"];\n"
//            "v6 -- v53\n"
//            "v6 -- v58\n"
//            "v7 [label=\"{a_{119}}\"; color=\"#4f1dd0\"];\n"
//            "v7 -- v53\n"
//            "v7 -- v58\n"
//            "v8 [label=\"{i_{100}}\"; color=\"#a78ee8\"];\n"
//            "v8 -- v21\n"
//            "v8 -- v26\n"
//            "v9 [label=\"{i_{101}}\"; color=\"#a78ee8\"];\n"
//            "v9 -- v21\n"
//            "v9 -- v26\n"
//            "v10 [label=\"{i_{108}}\"; color=\"#a78ee8\"];\n"
//            "v10 -- v38\n"
//            "v10 -- v41\n"
//            "v11 [label=\"{i_{109}}\"; color=\"#a78ee8\"];\n"
//            "v11 -- v38\n"
//            "v11 -- v41\n"
//            "v12 [label=\"{i_{112}}\"; color=\"#a78ee8\"];\n"
//            "v12 -- v46\n"
//            "v12 -- v49\n"
//            "v13 [label=\"{i_{113}}\"; color=\"#a78ee8\"];\n"
//            "v13 -- v46\n"
//            "v13 -- v49\n"
//            "v14 [label=\"{i_{116}}\"; color=\"#a78ee8\"];\n"
//            "v14 -- v54\n"
//            "v14 -- v57\n"
//            "v15 [label=\"{i_{117}}\"; color=\"#a78ee8\"];\n"
//            "v15 -- v54\n"
//            "v15 -- v57\n"
//            "v16 [label=\"{κ_{104}}\"; color=\"#808f74\"];\n"
//            "v16 -- v29\n"
//            "v16 -- v34\n"
//            "v17 [label=\"{κ_{105}}\"; color=\"#808f74\"];\n"
//            "v17 -- v29\n"
//            "v17 -- v34\n"
//            "v18 [label=\"{κ_{106}}\"; color=\"#808f74\"];\n"
//            "v18 -- v30\n"
//            "v18 -- v33\n"
//            "v19 [label=\"{κ_{107}}\"; color=\"#808f74\"];\n"
//            "v19 -- v30\n"
//            "v19 -- v33\n"
//            "v20 [label=\"A\"; color=\"#5cd738\"];\n"
//            "v20 -- v23\n"
//            "v21 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v21 -- v23\n"
//            "v22 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v22 -- v23\n"
//            "v23 [label=\"bka\"; color=\"#5cd738\"];\n"
//            "v24 [label=\"ã\"; color=\"#8b0711\"];\n"
//            "v24 -- v27\n"
//            "v25 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v25 -- v27\n"
//            "v26 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v26 -- v27\n"
//            "v27 [label=\"bka\"; color=\"#8b0711\"];\n"
//            "v28 [label=\"g\"; color=\"#4c04e1\"];\n"
//            "v28 -- v31\n"
//            "v29 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v29 -- v31\n"
//            "v30 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v30 -- v31\n"
//            "v31 [label=\"bka\"; color=\"#4c04e1\"];\n"
//            "v32 [label=\"ã\"; color=\"#8b0711\"];\n"
//            "v32 -- v35\n"
//            "v33 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v33 -- v35\n"
//            "v34 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v34 -- v35\n"
//            "v35 [label=\"bka\"; color=\"#8b0711\"];\n"
//            "v36 [label=\"t\"; color=\"#b4d687\"];\n"
//            "v36 -- v39\n"
//            "v37 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v37 -- v39\n"
//            "v38 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v38 -- v39\n"
//            "v39 [label=\"bka\"; color=\"#b4d687\"];\n"
//            "v40 [label=\"ã\"; color=\"#8b0711\"];\n"
//            "v40 -- v43\n"
//            "v41 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v41 -- v43\n"
//            "v42 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v42 -- v43\n"
//            "v43 [label=\"bka\"; color=\"#8b0711\"];\n"
//            "v44 [label=\"t\"; color=\"#b4d687\"];\n"
//            "v44 -- v47\n"
//            "v45 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v45 -- v47\n"
//            "v46 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v46 -- v47\n"
//            "v47 [label=\"bka\"; color=\"#b4d687\"];\n"
//            "v48 [label=\"ã\"; color=\"#8b0711\"];\n"
//            "v48 -- v51\n"
//            "v49 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v49 -- v51\n"
//            "v50 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v50 -- v51\n"
//            "v51 [label=\"bka\"; color=\"#8b0711\"];\n"
//            "v52 [label=\"t\"; color=\"#b4d687\"];\n"
//            "v52 -- v55\n"
//            "v53 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v53 -- v55\n"
//            "v54 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v54 -- v55\n"
//            "v55 [label=\"bka\"; color=\"#b4d687\"];\n"
//            "v56 [label=\"ã\"; color=\"#8b0711\"];\n"
//            "v56 -- v59\n"
//            "v57 [label=\"bra2a\"; color=\"#eaa2ab\"];\n"
//            "v57 -- v59\n"
//            "v58 [label=\"ket2a\"; color=\"#5a8fd3\"];\n"
//            "v58 -- v59\n"
//            "v59 [label=\"bka\"; color=\"#8b0711\"];\n"
//            "}\n");

    // compute automorphism group
    {
      bliss::Stats stats;
      graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

      std::vector<std::vector<unsigned int>> aut_generators;
      auto save_aut = [&aut_generators](const unsigned int n,
                                        const unsigned int* aut) {
        aut_generators.emplace_back(aut, aut + n);
      };
      graph->find_automorphisms(stats, &bliss::aut_hook<decltype(save_aut)>, &save_aut);

      // print automorphisms
      auto print_auts = [&aut_generators](auto&& stream, auto&& vlabels,
                                          bool use_labels) {
        ranges::for_each(aut_generators, [&stream, &vlabels,
                                          &use_labels](auto&& gen) {

          // see bliss::print_permutation
          auto print = [&stream, &vlabels,
                        &use_labels](const std::vector<unsigned int>& perm) {
            const unsigned int offset = 0;
            const unsigned int N = perm.size();
            for (unsigned int i = 0; i < N; i++) {
              unsigned int j = perm[i];
              if (j == i) continue;
              bool is_first = true;
              while (j != i) {
                if (j < i) {
                  is_first = false;
                  break;
                }
                j = perm[j];
              }
              if (!is_first) continue;
              stream << "("
                     << (use_labels ? vlabels.at(i)
                                    : std::to_wstring(i + offset))
                     << ",";
              j = perm[i];
              while (j != i) {
                stream << (use_labels ? vlabels.at(j)
                                      : std::to_wstring(j + offset));
                j = perm[j];
                if (j != i) stream << ",";
              }
              stream << ")";
            }
          };

          print(gen);
          stream << std::endl;

        });
      };
      std::basic_ostringstream<wchar_t> oss;
      print_auts(oss, vlabels, false);
      //std::wcout << oss.str() << std::endl;
      // the print_auts output seems to be platform-dependent, e.g. g++-8 on Linux
      // TODO investigate the source of platform dependence
//      REQUIRE(oss.str() ==
//              L"(18,19)\n"
//              "(16,17)\n"
//              "(0,1)\n"
//              "(8,9)\n"
//              "(2,3)\n"
//              "(10,11)\n"
//              "(4,5)\n"
//              "(6,7)\n"
//              "(12,13)\n"
//              "(14,15)\n"
//              "(4,6)(5,7)(12,14)(13,15)(44,52)(45,53)(46,54)(47,55)(48,56)(49,57)(50,58)(51,59)\n"
//              "(2,4)(3,5)(10,12)(11,13)(36,44)(37,45)(38,46)(39,47)(40,48)(41,49)(42,50)(43,51)\n");
      // change to 1 to user vertex labels rather than indices
      if (0) {
        std::basic_ostringstream<wchar_t> oss2;
        print_auts(oss2, vlabels, true);
        std::wcout << oss2.str() << std::endl;
      }
    }

  }  // SECTION("bliss graph")

}  // TEST_CASE("Tensor")
