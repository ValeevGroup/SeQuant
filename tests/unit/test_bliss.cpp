//
// Created by Eduard Valeyev on 2/26/19.
//

#include "catch.hpp"

#include "../../external/bliss/graph.hh"
#include "../../external/bliss/utils.hh"

#include <cstdio>
#include <iostream>
#include <numeric>
#include <vector>

constexpr const bool use_colors = true;

/**
 * The hook function that prints the found automorphisms.
 * \a param must be a file descriptor (FILE *).
 */
static void report_aut(void* param, const unsigned int n,
                       const unsigned int* aut) {
  assert(param);
  fprintf((FILE*)param, "Generator: ");
  bliss::print_permutation((FILE*)param, n, aut, 1);
  fprintf((FILE*)param, "\n");
}

TEST_CASE("Bliss", "[elements]") {

  SECTION("basic operation") {
      const int n = 16;  // # of vertices

      bliss::Graph sg1(n), sg2(n);

      auto add_edges = [](auto& sg, int v1, std::initializer_list<int> v2s) {
        assert(v1 >= 0 && v1 < sg.get_nof_vertices());
        for (auto&& v2 : v2s) {
          assert(v2 >= 0 && v2 < sg.get_nof_vertices());
          sg.add_edge(v1, v2);
        }
      };

      auto assign_colors =
          [](auto& sg,
             std::initializer_list<std::initializer_list<int>> color_partitions) {
            int partition_cnt = 0;
            for (auto&& partition : color_partitions) {
              for (auto&& v : partition) {
                sg.change_color(v, partition_cnt);
              }
              ++partition_cnt;
            }
          };

      /* Now make the two graphs */

      // 1. tensor vertices: 4 antisymmetric tensors = 8 vertices, each tensor has
      // order 3 => each vertex has degree 3
      // 2. line vertices: 8 lines = 8 vertices, each vertex has degree 2
      // b1 -> 0, k1 -> 1, b2 -> 2, k2 -> 3, b3 -> 4, k3 -> 5, b4 -> 6, k4 -> 7, i1
      // -> 8, i2 -> 9, ... a1 -> 12, a2 -> 13, ...
      const char* labels[] = {"b1", "k1", "b2", "k2", "b3", "k3", "b4", "k4",
                              "i1", "i2", "i3", "i4", "a1", "a2", "a3", "a4"};

      // "make" tensors in both graphs
      for (auto g = 0; g != 2; ++g) {
        auto& sg = g == 0 ? sg1 : sg2;
        sg.add_edge(0, 1);
        sg.add_edge(2, 3);
        sg.add_edge(4, 5);
        sg.add_edge(6, 7);
      }

      // Graph1

      //  i1 -- { b1 k3 }
      add_edges(sg1, 8, {0, 5});
      //  i2 -- { b1 k3 }
      add_edges(sg1, 9, {0, 5});
      //  i3 -- { b2 k4 }
      add_edges(sg1, 10, {2, 7});
      //  i4 -- { b2 k4 }
      add_edges(sg1, 11, {2, 7});
      //  a1 -- { k1 b3 }
      add_edges(sg1, 12, {1, 4});
      //  a2 -- { k1 b4 }
      add_edges(sg1, 13, {1, 6});
      //  a3 -- { k2 b3 }
      add_edges(sg1, 14, {3, 4});
      //  a4 -- { k2 b4 }
      add_edges(sg1, 15, {3, 6});

      // Graph2

      //    i1 -- { b1 k4 }
      add_edges(sg2, 8, {0, 7});
      //    i2 -- { b1 k4 }
      add_edges(sg2, 9, {0, 7});
      //    i3 -- { b2 k3 }
      add_edges(sg2, 10, {2, 5});
      //    i4 -- { b2 k3 }
      add_edges(sg2, 11, {2, 5});
      //    a1 -- { k1 b3 }
      add_edges(sg2, 12, {1, 4});
      //    a2 -- { k1 b4 }
      add_edges(sg2, 13, {1, 6});
      //    a3 -- { k2 b3 }
      add_edges(sg2, 14, {3, 4});
      //    a4 -- { k2 b4 }
      add_edges(sg2, 15, {3, 6});

      // coloring does not matter in this example, but in general needs to be used
      // to reduce cost and avoid accidental false symmetries
      if (use_colors) {
        // colors are the same for both graphs since the vertices are identically
        // ordered
        assign_colors(
            sg1,
            {{0}, {1}, {2}, {3}, {4, 6}, {5, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}});
        assign_colors(
            sg2,
            {{0}, {1}, {2}, {3}, {4, 6}, {5, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}});
      }

      sg1.set_splitting_heuristic(bliss::Graph::shs_fsm);
      sg2.set_splitting_heuristic(bliss::Graph::shs_fsm);

      bliss::Graph *cg1, *cg2;

      for (auto g = 0; g != 2; ++g) {
        auto& sg = g == 0 ? sg1 : sg2;
        auto*& cg = g == 0 ? cg1 : cg2;
        bliss::Stats stats;
        //    const unsigned int* cl = sg.canonical_form(stats, &report_aut,
        //    stdout);
        const unsigned int* cl = sg.canonical_form(stats, nullptr, nullptr);

//        fprintf(stdout, "Canonical labeling for graph %d: ", g);
//        bliss::print_permutation(stdout, sg.get_nof_vertices(), cl, 1);
//        fprintf(stdout, "\n");

        cg = sg.permute(cl);
      }

      REQUIRE(cg1->cmp(*cg2) == 0);
    }

}  // TEST_CASE("Bliss")