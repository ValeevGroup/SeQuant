//
// Created by Eduard Valeyev on 2019-03-27.
//

#ifndef SEQUANT_BLISS_HPP
#define SEQUANT_BLISS_HPP

#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/algorithm/for_each.hpp>

namespace bliss {

/// Generic wrapper for a Callable to be used as a hook given to
/// bliss::AbstractGraph::find_automorphisms
/// @tparam Callable a callable type for which Callable(n,aut) is valid and
/// returns void
/// @param callable_ptr_void pointer to the Callable object
/// @param n the number of vertices
/// @param aut an automorphism generator (permutation, specified as a permuted
/// list of vertex indices)
template <typename Callable>
void aut_hook(void* callable_ptr_void, const unsigned int n,
              const unsigned int* aut) {
  Callable* callable_ptr = reinterpret_cast<Callable*>(callable_ptr_void);
  (*callable_ptr)(n, aut);
}

/// prints automorphism generators to a stream
/// \tparam PermutationSequence a sequence of permutations, each represented as
/// a directly-addressable sequence \tparam Stream a standard stream type
/// \tparam VectorOfStrings a vector of strings
/// \param aut_generators the sequence of automorphism generators
/// \param stream the stream to print to
/// \param vlabels the labels to use to denote vertices; if empty, will use
/// vertex ordinals
template <typename PermutationSequence, typename Stream,
          typename VectorOfStrings>
void print_auts(const PermutationSequence& aut_generators, Stream& stream,
                const VectorOfStrings& vlabels) {
  ranges::for_each(aut_generators, [&stream, &vlabels](auto&& gen) {
    // see bliss::print_permutation
    auto print = [&stream, &vlabels,
                  use_labels = !vlabels.empty()](const auto& perm) {
      const unsigned int offset = 0;
      const unsigned int N = perm.size();
      if (!vlabels.empty()) assert(vlabels.size() == N);
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
               << (use_labels ? vlabels.at(i) : std::to_wstring(i + offset))
               << ",";
        j = perm[i];
        while (j != i) {
          stream << (use_labels ? vlabels.at(j) : std::to_wstring(j + offset));
          j = perm[j];
          if (j != i) stream << ",";
        }
        stream << ")";
      }
    };

    print(gen);
    stream << std::endl;
  });
}

}  // namespace bliss

#endif  // SEQUANT_BLISS_HPP
