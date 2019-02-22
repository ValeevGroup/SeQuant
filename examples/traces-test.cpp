//
// nautyex7.c modified by Eduard Valeyev on 2019-02-21.
//

/* This program demonstrates how an isomorphism is found between
   two diagrams. Uses Traces.
*/

extern "C" {
#include "traces.h"
}

#include <iostream>
#include <numeric>
#include <vector>

constexpr const bool use_colors = true;

int main(int argc, char* argv[]) {
  DYNALLSTAT(int, lab1, lab1_sz);
  DYNALLSTAT(int, lab2, lab2_sz);
  DYNALLSTAT(int, ptn1, ptn1_sz);
  DYNALLSTAT(int, ptn2, ptn2_sz);
  DYNALLSTAT(int, orbits, orbits_sz);
  DYNALLSTAT(int, map, map_sz);
  static DEFAULTOPTIONS_TRACES(options);
  TracesStats stats;
  /* Declare and initialize sparse graph structures */
  SG_DECL(sg1);
  SG_DECL(sg2);
  SG_DECL(cg1);
  SG_DECL(cg2);

  /* Select option for canonical labelling */

  options.getcanon = TRUE;

  // setup

  const int n = 16;   // # of vertices
  const int ne = 20;  // # of edges
  const int m = SETWORDSNEEDED(n);
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  DYNALLOC1(int, lab1, lab1_sz, n, "malloc");
  DYNALLOC1(int, lab2, lab2_sz, n, "malloc");
  DYNALLOC1(int, ptn1, ptn1_sz, n, "malloc");
  DYNALLOC1(int, ptn2, ptn2_sz, n, "malloc");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
  DYNALLOC1(int, map, map_sz, n, "malloc");

  /// makes a sequences of vertices of degrees given by @c vertex_degrees
  auto make_vertices = [](auto& sg, const std::vector<int>& vertex_degrees) {
    auto nvertices = vertex_degrees.size(); /* Number of vertices */
    auto nedges = std::accumulate(begin(vertex_degrees), end(vertex_degrees),
                                  0); /* Number of edges .. Traces only handles
                                         undirected graphs, so de=e */
    SG_ALLOC(sg, nvertices, nedges, "malloc");
    sg.nv = nvertices;
    sg.nde = nedges;

    size_t i = 0;
    for (auto&& deg : vertex_degrees) {
      assert(deg >= 0);
      sg.d[i] = deg;
      sg.v[i] = i == 0 ? 0 : sg.v[i - 1] + sg.d[i - 1];
      ++i;
    }

    // fill edges with -1 to make possible adding edges one at a time
    std::fill(sg.e, sg.e + nedges, -1);
  };

  /// attaches v1 to v2s
  auto add_edges = [](auto& sg, int v1, std::initializer_list<int> v2s) {
    assert(v1 >= 0 && v1 < sg.nv);
    const auto v1_edge_offset = sg.v[v1];
    const auto v1_degree = sg.d[v1];

    // find first free edge of v1
    int v1_edge = v1_edge_offset;
    while (sg.e[v1_edge] != -1) {
      // make sure we are not adding an edge already added
      assert(std::find(begin(v2s), end(v2s), sg.e[v1_edge]) == end(v2s));
      ++v1_edge;
    }

    for (auto&& v2 : v2s) {
      // make sure we have not run out of edges for v1
      assert(v1_edge < v1_edge_offset + v1_degree);

      assert(v2 >= 0 && v2 < sg.nv);

      const auto v2_edge_offset = sg.v[v2];
      const auto v2_degree = sg.d[v2];

      // find first free edge of v2
      int v2_edge = v2_edge_offset;
      while (sg.e[v2_edge] != -1) {
        // make sure we are not adding an edge already added
        assert(sg.e[v2_edge] != v1);
        ++v2_edge;
        // make sure we have not run out of edges
        assert(v2_edge < v2_edge_offset + v2_degree);
      }

      std::cout << "adding edge#" << (v1_edge - v1_edge_offset) << " of vertex "
                << v1 << " and edge#" << (v2_edge - v2_edge_offset)
                << " of vertex " << v2 << std::endl;

      sg.e[v1_edge] = v2;
      sg.e[v2_edge] = v1;
      ++v1_edge;
    }
  };

  /// attaches v to ve
  auto add_edge = [&](auto& sg, int v, int ve) { add_edges(sg, v, {ve}); };

  /// attaches v to ve
  auto validate_edges = [&](auto& sg) {
    auto it = std::find(sg.e, sg.e + sg.nde, -1);
    auto d = it - sg.e;
    assert(d == sg.nde);
  };

  auto assign_colors =
      [&](auto& sg, auto& lab, auto& ptn,
          std::initializer_list<std::initializer_list<int>> partitions) {
        size_t lab_cnt = 0;
        for (auto&& part : partitions) {
          const auto part_size = std::size(part);
          size_t cnt = 1;
          for (auto&& i : part) {
            assert(lab_cnt < sg.nv);
            lab[lab_cnt] = i;
            ptn[lab_cnt] = cnt == part_size ? 0 : 1;
            ++cnt;
            ++lab_cnt;
          }
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

  make_vertices(sg1, {3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2});
  make_vertices(sg2, {3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2});

  // "make" tensors in both graphs
  for (auto g = 0; g != 2; ++g) {
    auto& sg = g == 0 ? sg1 : sg2;
    add_edge(sg, 0, 1);
    add_edge(sg, 2, 3);
    add_edge(sg, 4, 5);
    add_edge(sg, 6, 7);
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

  validate_edges(sg1);
  validate_edges(sg2);

  // coloring does not matter in this example, but in general needs to be used
  // to reduce cost and avoid accidental false symmetries
  if (use_colors) {
    options.defaultptn = 0;
    // colors are the same for both graphs since the vertices are identically
    // ordered
    assign_colors(
        sg1, lab1, ptn1,
        {{0}, {1}, {2}, {3}, {4, 6}, {5, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}});
    assign_colors(
        sg2, lab2, ptn2,
        {{0}, {1}, {2}, {3}, {4, 6}, {5, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}});
  } else {
    options.defaultptn = 1;
  }

  /* Label sg1, result in cg1 and labelling in lab1; similarly sg2.
     It is not necessary to pre-allocate space in cg1 and cg2, but
     they have to be initialised as we did above.  */

  Traces(&sg1, lab1, ptn1, orbits, &options, &stats, &cg1);
  Traces(&sg2, lab2, ptn2, orbits, &options, &stats, &cg2);

  /* Compare canonically labeled graphs */

  if (aresame_sg(&cg1, &cg2)) {
    printf("Isomorphic.\n");
    if (n <= 1000) {
      /* Write the isomorphism.  For each i, vertex lab1[i]
         of sg1 maps onto vertex lab2[i] of sg2.  We compute
         the map in order of labelling because it looks better. */

      printf("graph 1 canonical labels:\n");
      for (auto i = 0; i < n; ++i) printf(" %s=>%d", labels[i], lab1[i]);
      printf("\ngraph 2 canonical labels:\n");
      for (auto i = 0; i < n; ++i) printf(" %s=>%d", labels[i], lab2[i]);

      for (auto i = 0; i < n; ++i) map[lab1[i]] = lab2[i];
      printf("\ngraph 1->2 automorphism:\n");
      for (auto i = 0; i < n; ++i) printf(" %s=>%s", labels[i], labels[map[i]]);
      printf("\n");
    }
  } else
    printf("Not isomorphic.\n");

  exit(0);
}
