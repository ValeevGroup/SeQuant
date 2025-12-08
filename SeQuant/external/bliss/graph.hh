#ifndef BLISS_GRAPH_HH
#define BLISS_GRAPH_HH

/*
  Copyright (c) 2003-2015 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * \namespace bliss
 * The namespace bliss contains all the classes and functions of the bliss
 * tool except for the C programming language API.
 */
namespace bliss {
class AbstractGraph;
}

#include <array>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <locale>
#include <map>
#include <optional>
#include <vector>

#include "bignum.hh"
#include "heap.hh"
#include "kqueue.hh"
#include "kstack.hh"
#include "orbit.hh"
#include "partition.hh"
#include "uintseqhash.hh"

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/iota.hpp>

#include <boost/container/small_vector.hpp>

namespace bliss {

/**
 * \brief Statistics returned by the bliss search algorithm.
 */
class Stats {
  friend class AbstractGraph;
  /** \internal The size of the automorphism group. */
  BigNum group_size;
  /** \internal An approximation (due to possible overflows) of
   * the size of the automorphism group. */
  long double group_size_approx;
  /** \internal The number of nodes in the search tree. */
  long unsigned int nof_nodes;
  /** \internal The number of leaf nodes in the search tree. */
  long unsigned int nof_leaf_nodes;
  /** \internal The number of bad nodes in the search tree. */
  long unsigned int nof_bad_nodes;
  /** \internal The number of canonical representative updates. */
  long unsigned int nof_canupdates;
  /** \internal The number of generator permutations. */
  long unsigned int nof_generators;
  /** \internal The maximal depth of the search tree. */
  unsigned long int max_level;
  /** */
  void reset() {
    group_size.assign(1);
    group_size_approx = 1.0;
    nof_nodes = 0;
    nof_leaf_nodes = 0;
    nof_bad_nodes = 0;
    nof_canupdates = 0;
    nof_generators = 0;
    max_level = 0;
  }

 public:
  Stats() { reset(); }
  /** Print the statistics. */
  size_t print(FILE* const fp) const {
    size_t r = 0;
    r += fprintf(fp, "Nodes:          %lu\n", nof_nodes);
    r += fprintf(fp, "Leaf nodes:     %lu\n", nof_leaf_nodes);
    r += fprintf(fp, "Bad nodes:      %lu\n", nof_bad_nodes);
    r += fprintf(fp, "Canrep updates: %lu\n", nof_canupdates);
    r += fprintf(fp, "Generators:     %lu\n", nof_generators);
    r += fprintf(fp, "Max level:      %lu\n", max_level);
    r += fprintf(fp, "|Aut|:          ") + group_size.print(fp) +
         fprintf(fp, "\n");
    fflush(fp);
    return r;
  }
  /** An approximation (due to possible overflows/rounding errors) of
   * the size of the automorphism group. */
  long double get_group_size_approx() const { return group_size_approx; }
  /** The number of nodes in the search tree. */
  long unsigned int get_nof_nodes() const { return nof_nodes; }
  /** The number of leaf nodes in the search tree. */
  long unsigned int get_nof_leaf_nodes() const { return nof_leaf_nodes; }
  /** The number of bad nodes in the search tree. */
  long unsigned int get_nof_bad_nodes() const { return nof_bad_nodes; }
  /** The number of canonical representative updates. */
  long unsigned int get_nof_canupdates() const { return nof_canupdates; }
  /** The number of generator permutations. */
  long unsigned int get_nof_generators() const { return nof_generators; }
  /** The maximal depth of the search tree. */
  unsigned long int get_max_level() const { return max_level; }
};

/**
 * \brief An abstract base class for different types of graphs.
 */
class AbstractGraph {
  friend class Partition;

 public:
  AbstractGraph() = default;
  virtual ~AbstractGraph();

  /**
   * Set the verbose output level for the algorithms.
   * \param level  the level of verbose output, 0 means no verbose output
   */
  void set_verbose_level(const unsigned int level);

  /**
   * Set the file stream for the verbose output.
   * \param fp  the file stream; if null, no verbose output is written
   */
  void set_verbose_file(FILE* const fp);

  /**
   * Add a new vertex with color \a color in the graph and return its index.
   */
  virtual unsigned int add_vertex(const unsigned int color = 0) = 0;

  /**
   * Add an edge between vertices \a source and \a target.
   * Duplicate edges between vertices are ignored but try to avoid introducing
   * them in the first place as they are not ignored immediately but will
   * consume memory and computation resources for a while.
   */
  virtual void add_edge(const unsigned int source,
                        const unsigned int target) = 0;

  /**
   * Change the color of the vertex \a vertex to \a color.
   */
  virtual void change_color(const unsigned int vertex,
                            const unsigned int color) = 0;

  /**
   * Check whether \a perm is an automorphism of this graph.
   * Unoptimized, mainly for debugging purposes.
   */
  virtual bool is_automorphism(const std::vector<unsigned int>& perm) const;

  /** Activate/deactivate failure recording.
   * May not be called during the search, i.e. from an automorphism reporting
   * hook function.
   * \param active  if true, activate failure recording, deactivate otherwise
   */
  void set_failure_recording(const bool active) {
    assert(!in_search);
    opt_use_failure_recording = active;
  }

  /** Activate/deactivate component recursion.
   * The choice affects the computed canonical labelings;
   * therefore, if you want to compare whether two graphs are isomorphic by
   * computing and comparing (for equality) their canonical versions,
   * be sure to use the same choice for both graphs.
   * May not be called during the search, i.e. from an automorphism reporting
   * hook function.
   * \param active  if true, activate component recursion, deactivate otherwise
   */
  void set_component_recursion(const bool active) {
    assert(!in_search);
    opt_use_comprec = active;
  }

  /**
   * Return the number of vertices in the graph.
   */
  virtual unsigned int get_nof_vertices() const = 0;

  /**
   * Return a new graph that is the result of applying the permutation \a perm
   * to this graph. This graph is not modified.
   * \a perm must contain N=this.get_nof_vertices() elements and be a bijection
   * on {0,1,...,N-1}, otherwise the result is undefined or a segfault.
   */
  virtual AbstractGraph* permute(const unsigned int* const perm) const = 0;
  virtual AbstractGraph* permute(
      const std::vector<unsigned int>& perm) const = 0;

  /**
   * Find a set of generators for the automorphism group of the graph.
   * The function \a hook (if non-null) is called each time a new generator
   * for the automorphism group is found.
   * The first argument \a user_param for the hook is the
   * \a hook_user_param given below,
   * the second argument \a n is the length of the automorphism (equal to
   * get_nof_vertices()) and
   * the third argument \a aut is the automorphism
   * (a bijection on {0,...,get_nof_vertices()-1}).
   * The memory for the automorphism \a aut will be invalidated immediately
   * after the return from the hook function;
   * if you want to use the automorphism later, you have to take a copy of it.
   * Do not call any member functions in the hook.
   * The search statistics are copied in \a stats.
   */
  void find_automorphisms(Stats& stats,
                          void (*hook)(void* user_param, unsigned int n,
                                       const unsigned int* aut),
                          void* hook_user_param);

  /**
   * Otherwise the same as find_automorphisms() except that
   * a canonical labeling of the graph (a bijection on
   * {0,...,get_nof_vertices()-1}) is returned.
   * The memory allocated for the returned canonical labeling will remain
   * valid only until the next call to a member function with the exception
   * that constant member functions (for example, bliss::Graph::permute()) can
   * be called without invalidating the labeling.
   * To compute the canonical version of an undirected graph, call this
   * function and then bliss::Graph::permute() with the returned canonical
   * labeling.
   * Note that the computed canonical version may depend on the applied version
   * of bliss as well as on some other options (for instance, the splitting
   * heuristic selected with bliss::Graph::set_splitting_heuristic()).
   */
  const unsigned int* canonical_form(Stats& stats,
                                     void (*hook)(void* user_param,
                                                  unsigned int n,
                                                  const unsigned int* aut),
                                     void* hook_user_param);

  /**
   * Write the graph to a file in a variant of the DIMACS format.
   * See the <A href="http://www.tcs.hut.fi/Software/bliss/">bliss website</A>
   * for the definition of the file format.
   * Note that in the DIMACS file the vertices are numbered from 1 to N while
   * in this C++ API they are from 0 to N-1.
   * Thus the vertex n in the file corresponds to the vertex n-1 in the API.
   * \param fp  the file stream where the graph is written
   */
  virtual void write_dimacs(FILE* const fp) = 0;

  /**
   * Write the graph to a file in the graphviz dotty format.
   * \param fp  the file stream where the graph is written
   */
  virtual void write_dot(FILE* const fp) = 0;

  /**
   * Write the graph in a file in the graphviz dotty format.
   * Do nothing if the file cannot be written.
   * \param file_name  the name of the file to which the graph is written
   */
  virtual void write_dot(const char* const file_name) = 0;

  /**
   * Get a 32-bit hash value for the graph.
   * \return  the hash value
   * \note this is used by all graph algorithms
   */
  virtual unsigned int get_hash() = 0;

  /**
 * Get a 64-bit hash value for the graph.
 * \return  the hash value
 * \note this can be used to a wider range of values than get_hash()
 */
  virtual std::uint64_t get_hash64() = 0;

  /**
   * Disable/enable the "long prune" method.
   * The choice affects the computed canonical labelings;
   * therefore, if you want to compare whether two graphs are isomorphic by
   * computing and comparing (for equality) their canonical versions,
   * be sure to use the same choice for both graphs.
   * May not be called during the search, i.e. from an automorphism reporting
   * hook function.
   * \param active  if true, activate "long prune", deactivate otherwise
   */
  void set_long_prune_activity(const bool active) {
    assert(!in_search);
    opt_use_long_prune = active;
  }

 protected:
  /** \internal
   * How much verbose output is produced (0 means none) */
  unsigned int verbose_level = 0;
  /** \internal
   * The output stream for verbose output. */
  FILE* verbstr = stdout;

 protected:
  /** \internal
   * The ordered partition used in the search algorithm. */
  Partition p;

  /** \internal
   * Whether the search for automorphisms and a canonical labeling is
   * in progress.
   */
  bool in_search = {};

  /** \internal
   * Is failure recording in use?
   */
  bool opt_use_failure_recording = true; /* Default value for using failure recording */
  /* The "tree-specific" invariant value for the point when current path
   * got different from the first path */
  unsigned int failure_recording_fp_deviation = {};

  /** \internal
   * Is component recursion in use?
   */
  bool opt_use_comprec = true; /* Default value for using component recursion */

  unsigned int refine_current_path_certificate_index = {};
  bool refine_compare_certificate = {};
  bool refine_equal_to_first = {};
  unsigned int refine_first_path_subcertificate_end = {};
  int refine_cmp_to_best = {};
  unsigned int refine_best_path_subcertificate_end = {};

  static const unsigned int CERT_SPLIT = 0;  // UINT_MAX;
  static const unsigned int CERT_EDGE = 1;   // UINT_MAX-1;
  /** \internal
   * Add a triple (v1,v2,v3) in the certificate.
   * May modify refine_equal_to_first and refine_cmp_to_best.
   * May also update eqref_hash and failure_recording_fp_deviation. */
  void cert_add(const unsigned int v1, const unsigned int v2,
                const unsigned int v3);

  /** \internal
   * Add a redundant triple (v1,v2,v3) in the certificate.
   * Can also just dicard the triple.
   * May modify refine_equal_to_first and refine_cmp_to_best.
   * May also update eqref_hash and failure_recording_fp_deviation. */
  void cert_add_redundant(const unsigned int x, const unsigned int y,
                          const unsigned int z);

  /**\internal
   * Is the long prune method in use?
   */
  bool opt_use_long_prune = true;  /* Default value for using "long prune" */
  /**\internal
   * Maximum amount of memory (in megabytes) available for
   * the long prune method
   */
  static const unsigned int long_prune_options_max_mem = 50;
  /**\internal
   * Maximum amount of automorphisms stored for the long prune method;
   * less than this is stored if the memory limit above is reached first
   */
  static const unsigned int long_prune_options_max_stored_auts = 100;

  unsigned int long_prune_max_stored_autss = {};
  std::vector<std::vector<bool>*> long_prune_fixed;
  std::vector<std::vector<bool>*> long_prune_mcrs;
  std::vector<bool> long_prune_temp;
  unsigned int long_prune_begin = {};
  unsigned int long_prune_end = {};
  /** \internal
   * Initialize the "long prune" data structures.
   */
  void long_prune_init();
  /** \internal
   * Release the memory allocated for "long prune" data structures.
   */
  void long_prune_deallocate();
  void long_prune_add_automorphism(const unsigned int* aut);
  std::vector<bool>& long_prune_get_fixed(const unsigned int index);
  std::vector<bool>& long_prune_allocget_fixed(const unsigned int index);
  std::vector<bool>& long_prune_get_mcrs(const unsigned int index);
  std::vector<bool>& long_prune_allocget_mcrs(const unsigned int index);
  /** \internal
   * Swap the i:th and j:th stored automorphism information;
   * i and j must be "in window, i.e. in [long_prune_begin,long_prune_end[
   */
  void long_prune_swap(const unsigned int i, const unsigned int j);

  /*
   * Data structures and routines for refining the partition p into equitable
   */
  Heap neighbour_heap;
  virtual bool split_neighbourhood_of_unit_cell(Partition::Cell*) = 0;
  virtual bool split_neighbourhood_of_cell(Partition::Cell* const) = 0;
  void refine_to_equitable();
  void refine_to_equitable(Partition::Cell* const unit_cell);
  void refine_to_equitable(Partition::Cell* const unit_cell1,
                           Partition::Cell* const unit_cell2);

  /** \internal
   * \return false if it was detected that the current certificate
   *         is different from the first and/or best (whether this is checked
   *         depends on in_search and refine_compare_certificate flags.
   */
  bool do_refine_to_equitable();

  unsigned int eqref_max_certificate_index;
  /** \internal
   * Whether eqref_hash is updated during equitable refinement process.
   */
  bool compute_eqref_hash = {};
  UintSeqHash eqref_hash;

  /** \internal
   * Check whether the current partition p is equitable.
   * Performance: very slow, use only for debugging purposes.
   */
  virtual bool is_equitable() const = 0;

  unsigned int* first_path_labeling = nullptr;
  unsigned int* first_path_labeling_inv = nullptr;
  Orbit first_path_orbits;
  unsigned int* first_path_automorphism = nullptr;

  unsigned int* best_path_labeling = nullptr;
  unsigned int* best_path_labeling_inv = nullptr;
  Orbit best_path_orbits;
  unsigned int* best_path_automorphism = nullptr;

  void update_labeling(unsigned int* const lab);
  void update_labeling_and_its_inverse(unsigned int* const lab,
                                       unsigned int* const lab_inv);
  void update_orbit_information(Orbit& o, const unsigned int* perm);

  void reset_permutation(unsigned int* perm);

  /* Mainly for debugging purposes */
  virtual bool is_automorphism(unsigned int* const perm);

  std::vector<unsigned int> certificate_current_path;
  std::vector<unsigned int> certificate_first_path;
  std::vector<unsigned int> certificate_best_path;

  unsigned int certificate_index = {};
  virtual void initialize_certificate() = 0;

  virtual void remove_duplicate_edges() = 0;
  virtual void make_initial_equitable_partition() = 0;
  virtual Partition::Cell* find_next_cell_to_be_splitted(
      Partition::Cell* cell) = 0;

  void search(const bool canonical, Stats& stats);

  void (*report_hook)(void* user_param, unsigned int n,
                      const unsigned int* aut) = nullptr;
  void* report_user_param = nullptr;

  /*
   *
   * Nonuniform component recursion (NUCR)
   *
   */

  /** The currently traversed component */
  unsigned int cr_level = {};

  /** \internal
   * The "Component End Point" data structure
   */
  class CR_CEP {
   public:
    /** At which level in the search was this CEP created */
    unsigned int creation_level;
    /** The current component has been fully traversed when the partition has
     * this many discrete cells left */
    unsigned int discrete_cell_limit;
    /** The component to be traversed after the current one */
    unsigned int next_cr_level;
    /** The next component end point */
    unsigned int next_cep_index;
    bool first_checked;
    bool best_checked;
  };
  /** \internal
   * A stack for storing Component End Points
   */
  std::vector<CR_CEP> cr_cep_stack;

  /** \internal
   * Find the first non-uniformity component at the component recursion
   * level \a level.
   * The component is stored in \a cr_component.
   * If no component is found, \a cr_component is empty.
   * Returns false if all the cells in the component recursion level \a level
   * were discrete.
   * Modifies the max_ival and max_ival_count fields of Partition:Cell
   * (assumes that they are 0 when called and
   *  quarantees that they are 0 when returned).
   */
  virtual bool nucr_find_first_component(const unsigned int level) = 0;
  virtual bool nucr_find_first_component(const unsigned int level,
                                         std::vector<unsigned int>& component,
                                         unsigned int& component_elements,
                                         Partition::Cell*& sh_return) = 0;
  /** \internal
   * The non-uniformity component found by nucr_find_first_component()
   * is stored here.
   */
  std::vector<unsigned int> cr_component;
  /** \internal
   * The number of vertices in the component \a cr_component
   */
  unsigned int cr_component_elements = {};
};

/**
 * \brief The class for undirected, vertex colored graphs.
 *
 * Multiple edges between vertices are not allowed (i.e., are ignored).
 */
class Graph : public AbstractGraph {
 public:
  friend struct ConstGraphCmp;

  /**
   * The possible splitting heuristics.
   * The selected splitting heuristics affects the computed canonical
   * labelings; therefore, if you want to compare whether two graphs
   * are isomorphic by computing and comparing (for equality) their
   * canonical versions, be sure to use the same splitting heuristics
   * for both graphs.
   */
  typedef enum {
    /** First non-unit cell.
     * Very fast but may result in large search spaces on difficult graphs.
     * Use for large but easy graphs. */
    shs_f = 0,
    /** First smallest non-unit cell.
     * Fast, should usually produce smaller search spaces than shs_f. */
    shs_fs,
    /** First largest non-unit cell.
     * Fast, should usually produce smaller search spaces than shs_f. */
    shs_fl,
    /** First maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_fm,
    /** First smallest maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_fsm,
    /** First largest maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_flm
  } SplittingHeuristic;

 protected:
  class Vertex {
   public:
    Vertex();
    ~Vertex();
    void add_edge(const unsigned int other_vertex);
    void remove_duplicate_edges(std::vector<bool>& tmp);
    void sort_edges();

    unsigned int color;
    boost::container::small_vector<unsigned int, 8> edges;
    unsigned int nof_edges() const { return edges.size(); }
  };
  std::vector<Vertex> vertices;
  void sort_edges();
  void remove_duplicate_edges();

  /** \internal
   * Partition independent invariant.
   * Returns the color of the vertex.
   * Time complexity: O(1).
   */
  static unsigned int vertex_color_invariant(const Graph* const g,
                                             const unsigned int v);
  /** \internal
   * Partition independent invariant.
   * Returns the degree of the vertex.
   * DUPLICATE EDGES MUST HAVE BEEN REMOVED BEFORE.
   * Time complexity: O(1).
   */
  static unsigned int degree_invariant(const Graph* const g,
                                       const unsigned int v);
  /** \internal
   * Partition independent invariant.
   * Returns 1 if there is an edge from the vertex to itself, 0 if not.
   * Time complexity: O(k), where k is the number of edges leaving the vertex.
   */
  static unsigned int selfloop_invariant(const Graph* const g,
                                         const unsigned int v);

  bool refine_according_to_invariant(unsigned int (*inv)(const Graph* const g,
                                                         const unsigned int v));

  /*
   * Routines needed when refining the partition p into equitable
   */
  bool split_neighbourhood_of_unit_cell(Partition::Cell*);
  bool split_neighbourhood_of_cell(Partition::Cell* const);

  /** \internal
   * \copydoc AbstractGraph::is_equitable() const
   */
  bool is_equitable() const;

  /* Splitting heuristics, documented in more detail in graph.cc */
  SplittingHeuristic sh;
  Partition::Cell* find_next_cell_to_be_splitted(Partition::Cell* cell);
  Partition::Cell* sh_first();
  Partition::Cell* sh_first_smallest();
  Partition::Cell* sh_first_largest();
  Partition::Cell* sh_first_max_neighbours();
  Partition::Cell* sh_first_smallest_max_neighbours();
  Partition::Cell* sh_first_largest_max_neighbours();

  void make_initial_equitable_partition();

  void initialize_certificate();

  bool is_automorphism(unsigned int* const perm);

  bool nucr_find_first_component(const unsigned int level);
  bool nucr_find_first_component(const unsigned int level,
                                 std::vector<unsigned int>& component,
                                 unsigned int& component_elements,
                                 Partition::Cell*& sh_return);

 public:
  /**
   * Create a new graph with \a N vertices and no edges.
   */
  Graph(const unsigned int N = 0);

  /**
   * Destroy the graph.
   */
  ~Graph();

  /**
   * Read the graph from the file \a fp in a variant of the DIMACS format.
   * See the <A href="http://www.tcs.hut.fi/Software/bliss/">bliss website</A>
   * for the definition of the file format.
   * Note that in the DIMACS file the vertices are numbered from 1 to N while
   * in this C++ API they are from 0 to N-1.
   * Thus the vertex n in the file corresponds to the vertex n-1 in the API.
   *
   * \param fp      the file stream for the graph file
   * \param errstr  if non-null, the possible error messages are printed
   *                in this file stream
   * \return        a new Graph object or 0 if reading failed for some
   *                reason
   */
  static Graph* read_dimacs(FILE* const fp, FILE* const errstr = stderr);

  /**
   * Write the graph to a file in a variant of the DIMACS format.
   * See the <A href="http://www.tcs.hut.fi/Software/bliss/">bliss website</A>
   * for the definition of the file format.
   */
  void write_dimacs(FILE* const fp);

  /**
   * \copydoc AbstractGraph::write_dot(FILE * const fp)
   */
  void write_dot(FILE* const fp);

  /**
   * \copydoc AbstractGraph::write_dot(const char * const file_name)
   */
  void write_dot(const char* const file_name);

  /// converts a color expressed as utin32_t to RGB color represented as 3 integers with 8 bits worth of information
  /// increase xsat (xsat <= 8) to try to increase saturation by this many bits
  /// @pre xsat <= 8
  static std::array<std::uint32_t,3> uint32_to_rgb(std::uint32_t i, unsigned int xsat = 0) {
    assert(xsat <= 8);
    // RGB has channel strides of {256^2, 256, 1} for a total range of [0, 256^3)
    // map 32-bit int to RGB by using channel strides of {s^2, s, 1}, where s=(2^32)^(1/3) = 1626
    constexpr std::uint32_t s = 1626;
    constexpr std::uint32_t ss = s*s;
    std::array<std::uint32_t,3> rgb_int; // 0 -> R, 1 -> G, 2 -> B
    rgb_int[0] = i/ss;
    rgb_int[1] = (i-rgb_int[0]*ss)/s;
    rgb_int[2] = (i-rgb_int[0]*ss - rgb_int[1]*s);
    for(int k=0; k!=3; ++k) rgb_int[k]=(rgb_int[k]*256)/s;  // rescale to 8 bits per channel

    // increase saturation, if needed
    if (xsat > 0) {
      // for each channel saturate independently
      for (int k = 0; k != 3; ++k) {
        int nsatbits = xsat;
        while (nsatbits > 0) {
          // do we have this many unset significant bits in this channel?
          bool can_sat = true;
          if ((rgb_int[k] >> (8 - nsatbits)) > 0) can_sat = false;
          if (can_sat) {
            rgb_int[k] = rgb_int[k] == 0 ? (1<<(nsatbits-1)) : (rgb_int[k] << nsatbits);
            break;
          } else
            nsatbits--;
        }
      }
    }

    return rgb_int;
  }

  /// options for generating dot file
  template <typename Char, typename Traits>
  struct DotOptions {
    /// vertex labels
    std::vector<std::basic_string<Char, Traits>> labels = {};
    /// vertex xlabels
    std::vector<std::optional<std::basic_string<Char, Traits>>> xlabels = {};
    /// vertex texlabels
    std::vector<std::optional<std::basic_string<Char, Traits>>> texlabels = {};
    /// if true, display colored vertices using color values
    /// to RGB colors;
    ///        if false, color value is appended to the vertex label; by default
    ///        set to true/false if vertex_labels are given/not given
    std::optional<bool> display_colors = std::nullopt;
    /// specifies the saturation applied to fillcolor of colored vertices
    /// 0 means use same fillcolor as the vertex color
    /// n>0 shifts each channel of RGB color left by n bits
    std::uint16_t fillcolor_saturation_nbits = 3;
    /// if nonnull, specifies a callable that maps vertex ordinal to subgraph ordinal; returning std::nullopt means vertex is not part of a subgraph
    std::function<std::optional<std::uint64_t>(std::size_t)> vertex_to_subgraph = {};
    /// extra content for node section
    std::basic_string<Char, Traits> nodeextras = {};
  };
  template <typename Char, typename Traits>
  static DotOptions<Char, Traits> make_default_dot_options() { return {};}

  /// @brief  writes a dot file, optionally using user-supplied labels and
  /// converting color values to colors
  /// @tparam Char a character type
  /// @tparam Traits a stream traits type
  /// @param os the output stream
  /// @param options options for generating dot file
  template <typename Char, typename Traits>
  void write_dot(std::basic_ostream<Char, Traits>& os,
                 DotOptions<Char,Traits> options = make_default_dot_options<Char, Traits>()) {
    using std::size;
    const auto nvertices = size(options.labels);
    const bool have_labels = nvertices > 0;
    const bool have_xlabels = size(options.xlabels) > 0;
    assert(!have_xlabels || size(options.xlabels) == nvertices);
    const bool have_texlabels = size(options.texlabels) > 0;
    assert(!have_texlabels || size(options.texlabels) == nvertices);
    const bool rgb_colors = options.display_colors.value_or(have_labels);

    remove_duplicate_edges();

    // converts int color to RGB string; increase xsat (xsat <= 8) to try to increase saturation by this many bits
    auto int_to_rgb = [](std::uint32_t i, unsigned int xsat = 0) {
      auto rgb_int8 = uint32_to_rgb(i, xsat);
      // map back to an int
      i = rgb_int8[0] * (1<<16) + rgb_int8[1] * (1<<8) + rgb_int8[2];

      std::basic_stringstream<Char> stream;
      // Set locale of this stream to C to avoid any kind of thousands separator
      stream.imbue(std::locale::classic());
      stream << std::setfill(Char('0')) << std::setw(6) << std::hex
             << i;
      return stream.str();
    };

    os << "graph g {\nnode [ style=filled, penwidth=2, margin=0" << (options.nodeextras.empty() ? "" : ",") << options.nodeextras << "];\n";

    using clusters_t = std::multimap<std::int64_t,std::size_t>;  // cluster id -> vertex ordinals
    clusters_t clusters;
    ranges::for_each(ranges::views::iota(0ul, vertices.size()),
                       [&](auto vertex_ordinal) {
                       // -1 is reserved
                       const std::int64_t cluster_ordinal = options.vertex_to_subgraph ? options.vertex_to_subgraph(vertex_ordinal).value_or(-1) : -1;
                       clusters.emplace(cluster_ordinal, vertex_ordinal);
                     });

    std::optional<std::int64_t> current_cluster_ordinal;
    for (auto& [cluster_ordinal, vertex_ordinal] : clusters) {
      auto vord = vertex_ordinal;
      Vertex& v = vertices.at(vertex_ordinal);

      // start of a new cluster? start new subgraph
      if (!current_cluster_ordinal.has_value() || (current_cluster_ordinal.has_value() && cluster_ordinal != current_cluster_ordinal.value())) {
        if (current_cluster_ordinal.value_or(-1) != -1) os << "}\n";
        current_cluster_ordinal = cluster_ordinal;
        if (current_cluster_ordinal.value_or(-1) != -1) {
          os << "subgraph cluster" << cluster_ordinal << " {\n";
        }
      }

      os << "v" << vord << " [ label=\"";
      if (have_labels) {
        assert(vord < nvertices);
        os << options.labels[vord];
      } else
        os << vord;
      if (have_xlabels) {
        assert(vord < nvertices);
        if (options.xlabels[vord].has_value()) {
          os << "\", xlabel=\"" << options.xlabels[vord].value();
        }
      }
      if (have_texlabels) {
        assert(vord < nvertices);
        if (options.texlabels[vord].has_value()) {
          os << "\", texlbl=\"" << options.texlabels[vord].value();
        }
      }
      if (rgb_colors) {
        auto color = int_to_rgb(v.color);
        auto satcolor = int_to_rgb(v.color, options.fillcolor_saturation_nbits);  // use saturated color for fill since transparency is not supported by all renderers (e.g. TikZ)
        os << "\", color=\"#" << color << "\", fillcolor=\"#" << satcolor << "\" ];\n";
      } else {
        os << ":" << v.color << "\"];\n";
      }
      for (auto ei = v.edges.begin();
           ei != v.edges.end(); ei++) {
        const unsigned int vord2 = *ei;

        bool declare_edge = false;
        // strangely vertices get assigned to clusters where their first edge is appears
        // so introduce all edges involving this vertex to vertices in this or later cluster
        if (options.vertex_to_subgraph) {
          const std::int64_t v2_cluster_ordinal = options.vertex_to_subgraph(vord2).value_or(-1);
          declare_edge = cluster_ordinal < v2_cluster_ordinal || (cluster_ordinal == v2_cluster_ordinal && vord < vord2);
        }
        else {
          // if not using clusters print edges from lesser to greater ordinal
          declare_edge = vord < vord2;
        }
        if (declare_edge)
          os << "v" << vord << " -- v" << vord2 << "\n";
      }
    }

    // close the last subgraph
    if (current_cluster_ordinal.has_value() && current_cluster_ordinal.value() != -1)
      os << "}\n";

    // end of graph
    os << "}" << std::endl;
  }

  /**
   * \copydoc AbstractGraph::is_automorphism(const std::vector<unsigned int>&
   * perm) const
   */
  bool is_automorphism(const std::vector<unsigned int>& perm) const;

  /**
   * \copydoc AbstractGraph::get_hash()
   */
  virtual unsigned int get_hash();

  /**
   * \copydoc AbstractGraph::get_hash64()
   */
  virtual std::uint64_t get_hash64();

  /**
   * Return the number of vertices in the graph.
   */
  unsigned int get_nof_vertices() const { return vertices.size(); }

  /**
   * \copydoc AbstractGraph::permute(const unsigned int* const perm) const
   */
  Graph* permute(const unsigned int* const perm) const;
  Graph* permute(const std::vector<unsigned int>& perm) const;

  /**
   * Add a new vertex with color \a color in the graph and return its index.
   */
  unsigned int add_vertex(const unsigned int color = 0);

  /**
   * Add an edge between vertices \a v1 and \a v2.
   * Duplicate edges between vertices are ignored but try to avoid introducing
   * them in the first place as they are not ignored immediately but will
   * consume memory and computation resources for a while.
   */
  void add_edge(const unsigned int v1, const unsigned int v2);

  /**
   * Change the color of the vertex \a vertex to \a color.
   */
  void change_color(const unsigned int vertex, const unsigned int color);

  /**
   * Compare this graph with the graph \a other.
   * Returns 0 if the graphs are equal, and a negative (positive) integer
   * if this graph is "smaller than" ("greater than", resp.) than \a other.
   */
  int cmp(Graph& other);

  /**
   * Set the splitting heuristic used by the automorphism and canonical
   * labeling algorithm.
   * The selected splitting heuristics affects the computed canonical
   * labelings; therefore, if you want to compare whether two graphs
   * are isomorphic by computing and comparing (for equality) their
   * canonical versions, be sure to use the same splitting heuristics
   * for both graphs.
   */
  void set_splitting_heuristic(const SplittingHeuristic shs) { sh = shs; }

protected:
  /**
   * Compare this graph with the graph \a other.
   * Returns 0 if the graphs are equal, and a negative (positive) integer
   * if this graph is "smaller than" ("greater than", resp.) than \a other.
   * Note that the const version will not perform removal of duplicate edges.
   */
  int cmp(const Graph& other) const;
};

/**
* WARNING: Use only if you are aware of the caveats regarding
* the use of the const-qualified Graph::cmp function.
*
* Accessor for the const Graph::cmp implementation
*/
struct ConstGraphCmp {
  static auto cmp(const Graph &lhs, const Graph &rhs) {
    return lhs.cmp(rhs);
  }
};

/**
 * \brief The class for directed, vertex colored graphs.
 *
 * Multiple edges between vertices are not allowed (i.e., are ignored).
 */
class Digraph : public AbstractGraph {
 public:
  /**
   * The possible splitting heuristics.
   * The selected splitting heuristics affects the computed canonical
   * labelings; therefore, if you want to compare whether two graphs
   * are isomorphic by computing and comparing (for equality) their
   * canonical versions, be sure to use the same splitting heuristics
   * for both graphs.
   */
  typedef enum {
    /** First non-unit cell.
     * Very fast but may result in large search spaces on difficult graphs.
     * Use for large but easy graphs. */
    shs_f = 0,
    /** First smallest non-unit cell.
     * Fast, should usually produce smaller search spaces than shs_f. */
    shs_fs,
    /** First largest non-unit cell.
     * Fast, should usually produce smaller search spaces than shs_f. */
    shs_fl,
    /** First maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_fm,
    /** First smallest maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_fsm,
    /** First largest maximally non-trivially connected non-unit cell.
     * Not so fast, should usually produce smaller search spaces than shs_f,
     * shs_fs, and shs_fl. */
    shs_flm
  } SplittingHeuristic;

 protected:
  class Vertex {
   public:
    Vertex();
    ~Vertex();
    void add_edge_to(const unsigned int dest_vertex);
    void add_edge_from(const unsigned int source_vertex);
    void remove_duplicate_edges(std::vector<bool>& tmp);
    void sort_edges();
    unsigned int color;
    boost::container::small_vector<unsigned int, 8> edges_out;
    boost::container::small_vector<unsigned int, 8> edges_in;
    unsigned int nof_edges_in() const { return edges_in.size(); }
    unsigned int nof_edges_out() const { return edges_out.size(); }
  };
  std::vector<Vertex> vertices;
  void remove_duplicate_edges();

  /** \internal
   * Partition independent invariant.
   * Returns the color of the vertex.
   * Time complexity: O(1).
   */
  static unsigned int vertex_color_invariant(const Digraph* const g,
                                             const unsigned int v);
  /** \internal
   * Partition independent invariant.
   * Returns the indegree of the vertex.
   * DUPLICATE EDGES MUST HAVE BEEN REMOVED BEFORE.
   * Time complexity: O(1).
   */
  static unsigned int indegree_invariant(const Digraph* const g,
                                         const unsigned int v);
  /** \internal
   * Partition independent invariant.
   * Returns the outdegree of the vertex.
   * DUPLICATE EDGES MUST HAVE BEEN REMOVED BEFORE.
   * Time complexity: O(1).
   */
  static unsigned int outdegree_invariant(const Digraph* const g,
                                          const unsigned int v);
  /** \internal
   * Partition independent invariant.
   * Returns 1 if there is an edge from the vertex to itself, 0 if not.
   * Time complexity: O(k), where k is the number of edges leaving the vertex.
   */
  static unsigned int selfloop_invariant(const Digraph* const g,
                                         const unsigned int v);

  /** \internal
   * Refine the partition \a p according to
   * the partition independent invariant \a inv.
   */
  bool refine_according_to_invariant(unsigned int (*inv)(const Digraph* const g,
                                                         const unsigned int v));

  /*
   * Routines needed when refining the partition p into equitable
   */
  bool split_neighbourhood_of_unit_cell(Partition::Cell* const);
  bool split_neighbourhood_of_cell(Partition::Cell* const);

  /** \internal
   * \copydoc AbstractGraph::is_equitable() const
   */
  bool is_equitable() const;

  /* Splitting heuristics, documented in more detail in the cc-file. */
  SplittingHeuristic sh;
  Partition::Cell* find_next_cell_to_be_splitted(Partition::Cell* cell);
  Partition::Cell* sh_first();
  Partition::Cell* sh_first_smallest();
  Partition::Cell* sh_first_largest();
  Partition::Cell* sh_first_max_neighbours();
  Partition::Cell* sh_first_smallest_max_neighbours();
  Partition::Cell* sh_first_largest_max_neighbours();

  void make_initial_equitable_partition();

  void initialize_certificate();

  bool is_automorphism(unsigned int* const perm);

  void sort_edges();

  bool nucr_find_first_component(const unsigned int level);
  bool nucr_find_first_component(const unsigned int level,
                                 std::vector<unsigned int>& component,
                                 unsigned int& component_elements,
                                 Partition::Cell*& sh_return);

 public:
  /**
   * Create a new directed graph with \a N vertices and no edges.
   */
  Digraph(const unsigned int N = 0);

  /**
   * Destroy the graph.
   */
  ~Digraph();

  /**
   * Read the graph from the file \a fp in a variant of the DIMACS format.
   * See the <A href="http://www.tcs.hut.fi/Software/bliss/">bliss website</A>
   * for the definition of the file format.
   * Note that in the DIMACS file the vertices are numbered from 1 to N while
   * in this C++ API they are from 0 to N-1.
   * Thus the vertex n in the file corresponds to the vertex n-1 in the API.
   * \param fp      the file stream for the graph file
   * \param errstr  if non-null, the possible error messages are printed
   *                in this file stream
   * \return        a new Digraph object or 0 if reading failed for some
   *                reason
   */
  static Digraph* read_dimacs(FILE* const fp, FILE* const errstr = stderr);

  /**
   * \copydoc AbstractGraph::write_dimacs(FILE * const fp)
   */
  void write_dimacs(FILE* const fp);

  /**
   * \copydoc AbstractGraph::write_dot(FILE *fp)
   */
  void write_dot(FILE* const fp);

  /**
   * \copydoc AbstractGraph::write_dot(const char * const file_name)
   */
  void write_dot(const char* const file_name);

  /**
   * \copydoc AbstractGraph::is_automorphism(const std::vector<unsigned int>&
   * perm) const
   */
  bool is_automorphism(const std::vector<unsigned int>& perm) const;

  /**
   * \copydoc AbstractGraph::get_hash()
   */
  virtual unsigned int get_hash();

  /**
   * \copydoc AbstractGraph::get_hash64()
   */
  virtual std::uint64_t get_hash64();

  /**
   * Return the number of vertices in the graph.
   */
  unsigned int get_nof_vertices() const { return vertices.size(); }

  /**
   * Add a new vertex with color 'color' in the graph and return its index.
   */
  unsigned int add_vertex(const unsigned int color = 0);

  /**
   * Add an edge from the vertex \a source to the vertex \a target.
   * Duplicate edges are ignored but try to avoid introducing
   * them in the first place as they are not ignored immediately but will
   * consume memory and computation resources for a while.
   */
  void add_edge(const unsigned int source, const unsigned int target);

  /**
   * Change the color of the vertex 'vertex' to 'color'.
   */
  void change_color(const unsigned int vertex, const unsigned int color);

  /**
   * Compare this graph with the graph \a other.
   * Returns 0 if the graphs are equal, and a negative (positive) integer
   * if this graph is "smaller than" ("greater than", resp.) than \a other.
   */
  int cmp(Digraph& other);

  /**
   * Set the splitting heuristic used by the automorphism and canonical
   * labeling algorithm.
   * The selected splitting heuristics affects the computed canonical
   * labelings; therefore, if you want to compare whether two graphs
   * are isomorphic by computing and comparing (for equality) their
   * canonical versions, be sure to use the same splitting heuristics
   * for both graphs.
   */
  void set_splitting_heuristic(SplittingHeuristic shs) { sh = shs; }

  /**
   * \copydoc AbstractGraph::permute(const unsigned int* const perm) const
   */
  Digraph* permute(const unsigned int* const perm) const;
  Digraph* permute(const std::vector<unsigned int>& perm) const;
};

}  // namespace bliss

#endif
