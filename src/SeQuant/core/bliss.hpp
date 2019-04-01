//
// Created by Eduard Valeyev on 2019-03-27.
//

#ifndef SEQUANT_BLISS_HPP
#define SEQUANT_BLISS_HPP

#include "../../../external/bliss/graph.hh"

namespace bliss {

/// Generic wrapper for a Callable to be used as a hook given to bliss::AbstractGraph::find_automorphisms
/// @tparam Callable a callable type for which Callable(n,aut) is valid and returns void
/// @param callable_ptr_void pointer to the Callable object
/// @param n the number of vertices
/// @param aut an automorphism generator (permutation, specified as a permuted list of vertex indices)
template <typename Callable>
void aut_hook(void* callable_ptr_void, const unsigned int n,
              const unsigned int *aut) {
  Callable* callable_ptr = reinterpret_cast<Callable*>(callable_ptr_void);
  (*callable_ptr)(n, aut);
}

}  // namespace bliss

#endif //SEQUANT_BLISS_HPP
