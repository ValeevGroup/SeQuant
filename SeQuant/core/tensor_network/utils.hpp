//
// Created to resolve duplicate symbol issues between tensor_network_v2 and
// tensor_network_v3 Contains shared utility functions and structs

#ifndef SEQUANT_TENSOR_NETWORK_UTILITIES_HPP
#define SEQUANT_TENSOR_NETWORK_UTILITIES_HPP

#include <SeQuant/core/abstract_tensor.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>
#include <type_traits>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/view/any_view.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {

inline bool tensors_commute(const AbstractTensor &lhs,
                            const AbstractTensor &rhs) {
  // tensors commute if their colors are different or either one of them
  // is a c-number
  return !(color(lhs) == color(rhs) && !is_cnumber(lhs) && !is_cnumber(rhs));
}

struct TensorBlockCompare {
  bool operator()(const AbstractTensor &lhs, const AbstractTensor &rhs) const {
    if (label(lhs) != label(rhs)) {
      return label(lhs) < label(rhs);
    }

    if (bra_rank(lhs) != bra_rank(rhs)) {
      return bra_rank(lhs) < bra_rank(rhs);
    }
    if (ket_rank(lhs) != ket_rank(rhs)) {
      return ket_rank(lhs) < ket_rank(rhs);
    }
    if (aux_rank(lhs) != aux_rank(rhs)) {
      return aux_rank(lhs) < aux_rank(rhs);
    }

    // Note: Accessing bra, ket and aux individually is a lot faster
    // than accessing the combined index() object
#define SEQUANT_CHECK_IDX_GROUP(group)                                  \
  auto lhs_##group = lhs._##group();                                    \
  auto rhs_##group = rhs._##group();                                    \
  auto lhs_##group##_end = lhs_##group.end();                           \
  auto rhs_##group##_end = rhs_##group.end();                           \
  for (auto lhs_it = lhs_##group.begin(), rhs_it = rhs_##group.begin(); \
       lhs_it != lhs_##group##_end && rhs_it != rhs_##group##_end;      \
       ++lhs_it, ++rhs_it) {                                            \
    if (lhs_it->space() != rhs_it->space()) {                           \
      return lhs_it->space() < rhs_it->space();                         \
    }                                                                   \
  }

    SEQUANT_CHECK_IDX_GROUP(bra);
    SEQUANT_CHECK_IDX_GROUP(ket);
    SEQUANT_CHECK_IDX_GROUP(aux);

    // Tensors are identical
    return false;
  }
};

/// Compares tensors based on their label and orders them according to the order
/// of the given cardinal tensor labels. If two tensors can't be discriminated
/// via their label, they are compared based on regular
/// AbstractTensor::operator< or based on their tensor block (the spaces of
/// their indices) - depending on the configuration. If this doesn't
/// discriminate the tensors, they are considered equal
template <typename CardinalLabels>
struct CanonicalTensorCompare {
  const CardinalLabels &labels;
  bool blocks_only;

  CanonicalTensorCompare(const CardinalLabels &labels, bool blocks_only)
      : labels(labels), blocks_only(blocks_only) {}

  void set_blocks_only(bool blocks_only) { this->blocks_only = blocks_only; }

  bool operator()(const AbstractTensorPtr &lhs_ptr,
                  const AbstractTensorPtr &rhs_ptr) const {
    assert(lhs_ptr);
    assert(rhs_ptr);
    const AbstractTensor &lhs = *lhs_ptr;
    const AbstractTensor &rhs = *rhs_ptr;

    if (!tensors_commute(lhs, rhs)) {
      return false;
    }

    const auto get_label = [](const auto &t) {
      if (label(t).back() == adjoint_label) {
        // grab base label if adjoint label is present
        return label(t).substr(0, label(t).size() - 1);
      }
      return label(t);
    };

    const auto lhs_it = std::find(labels.begin(), labels.end(), get_label(lhs));
    const auto rhs_it = std::find(labels.begin(), labels.end(), get_label(rhs));

    if (lhs_it != rhs_it) {
      // At least one of the tensors is a cardinal one
      // -> Order by the occurrence in the cardinal label list
      return std::distance(labels.begin(), lhs_it) <
             std::distance(labels.begin(), rhs_it);
    }

    // Either both are the same cardinal tensor or none is a cardinal tensor
    if (blocks_only) {
      TensorBlockCompare cmp;
      return cmp(lhs, rhs);
    } else {
      return lhs < rhs;
    }
  }

  template <typename AbstractTensorPtr_T>
  bool operator()(const AbstractTensorPtr_T &lhs_ptr,
                  const AbstractTensorPtr_T &rhs_ptr) const {
    return (*this)(lhs_ptr.first, rhs_ptr.first);
  }
};

template <typename ArrayLike, typename Permutation>
auto permute(const ArrayLike &vector, const Permutation &perm) {
  using std::size;
  auto sz = size(vector);
  std::decay_t<decltype(vector)> pvector(sz);
  for (size_t i = 0; i != sz; ++i) pvector[perm[i]] = vector[i];
  return pvector;
}

template <typename ReplacementMap>
void apply_index_replacements(AbstractTensor &tensor,
                              const ReplacementMap &replacements,
                              const bool self_consistent) {
#ifndef NDEBUG
  // assert that tensors' indices are not tagged since going to tag indices
  assert(ranges::none_of(
      indices(tensor), [](const Index &idx) { return idx.tag().has_value(); }));
#endif

  bool pass_mutated;
  do {
    pass_mutated = transform_indices(tensor, replacements);
  } while (self_consistent && pass_mutated);  // transform till stops changing

  reset_tags(tensor);
}

template <typename ArrayLike, typename ReplacementMap>
void apply_index_replacements(ArrayLike &tensors,
                              const ReplacementMap &replacements,
                              const bool self_consistent) {
  for (auto &tensor : tensors) {
    apply_index_replacements(*tensor, replacements, self_consistent);
  }
}

template <typename Container>
void order_to_indices(Container &container) {
  std::vector<std::size_t> indices;
  indices.resize(container.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::sort(indices.begin(), indices.end(),
            [&container](std::size_t lhs, std::size_t rhs) {
              return container[lhs] < container[rhs];
            });
  // Overwrite container contents with indices
  std::copy(indices.begin(), indices.end(), container.begin());
}

template <bool stable, typename Container, typename Comparator>
void sort_via_indices(Container &container, const Comparator &cmp) {
  std::vector<std::size_t> indices;
  indices.resize(container.size());
  std::iota(indices.begin(), indices.end(), 0);

  if constexpr (stable) {
    std::stable_sort(indices.begin(), indices.end(), cmp);
  } else {
    std::sort(indices.begin(), indices.end(), cmp);
  }

  // Bring elements in container into the order given by indices
  // (the association is container[k] = container[indices[k]])
  // -> implementation from https://stackoverflow.com/a/838789

  for (std::size_t i = 0; i < container.size(); ++i) {
    if (indices[i] == i) {
      // This element is already where it is supposed to be
      continue;
    }

    // Find the offset of the index pointing to i
    // -> since we are going to change the content of the vector at position i,
    // we have to update the index-mapping referencing i to point to the new
    // location of the element that used to be at position i
    std::size_t k;
    for (k = i + 1; k < container.size(); ++k) {
      if (indices[k] == i) {
        break;
      }
    }
    std::swap(container[i], container[indices[i]]);
    std::swap(indices[i], indices[k]);
  }
}

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_UTILITIES_HPP
