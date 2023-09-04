#ifndef SEQUANT_CORE_UTILITY_INDICES_HPP
#define SEQUANT_CORE_UTILITY_INDICES_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>

#include <set>
#include <vector>

namespace sequant {

namespace detail {

/// This function is equal to std::remove in case the given container
/// contains none or only a single occurrence of the given element.
/// If the given element is contained multiple times, only the first
/// occurrence is removed.
/// If it is known that there is only a single occurrence of the given
/// element, this function can be more efficient than std::remove as it
/// can terminate as soon as the element has been found instead of having
/// to traverse the entire container.
template <typename Container, typename Element>
void remove_one(Container& container, const Element& e) {
  auto iter = std::find(container.begin(), container.end(), e);

  if (iter != container.end()) {
    container.erase(iter);
  }
}

}  // namespace detail

/// A composite type for holding different named groups of indices
template <typename Container = std::vector<Index>>
struct IndexGroups {
  Container bra;
  Container ket;

  bool operator==(const IndexGroups<Container>& other) const {
    return bra == other.bra && ket == other.ket;
  }

  bool operator!=(const IndexGroups<Container>& other) const {
    return !(*this == other);
  }
};

template <typename Container = std::vector<Index>>
IndexGroups<Container> get_unique_indices(const ExprPtr& expr);

template <typename Container = std::vector<Index>>
IndexGroups<Container> get_unique_indices(const Constant&) {
  return {};
}

template <typename Container = std::vector<Index>>
IndexGroups<Container> get_unique_indices(const Variable&) {
  return {};
}

/// Obtains the set of unique (non-repeated) indices used in the given tensor
template <typename Container = std::vector<Index>>
IndexGroups<Container> get_unique_indices(const Tensor& tensor) {
  IndexGroups<Container> groups;
  std::set<Index> encounteredIndices;

  for (const Index& current : tensor.bra()) {
    if (encounteredIndices.find(current) == encounteredIndices.end()) {
      groups.bra.push_back(current);
      encounteredIndices.insert(current);
    } else {
      // There can't be indices in bra at this point, so we don't have to remove
      // from that
      detail::remove_one(groups.bra, current);
    }
  }

  for (const Index& current : tensor.ket()) {
    if (encounteredIndices.find(current) == encounteredIndices.end()) {
      groups.ket.push_back(current);
      encounteredIndices.insert(current);
    } else {
      detail::remove_one(groups.bra, current);
      detail::remove_one(groups.ket, current);
    }
  }

  return groups;
}

/// Obtains the set of unique (non-repeated) indices used in the given sum
/// The assumption here is that the set of unique indices of different addends
/// is the same (which it should be anyway as otherwise it wouldn't make sense
/// to add these terms in the first place)
template <typename Container = std::vector<Index>>
IndexGroups<Container> get_unique_indices(const Sum& sum) {
  // In order for the sum to be valid, all summands must have the same
  // external indices, so it suffices to look only at the first one
  return sum.summands().empty() ? IndexGroups<Container>{}
                                : get_unique_indices<Container>(sum.summand(0));
}

/// Obtains the set of unique (non-repeated) indices used in the given product
template <typename Container = std::vector<Index>>
IndexGroups<Container> get_unique_indices(const Product& product) {
  std::set<Index> encounteredIndices;
  IndexGroups<Container> groups;

  for (const ExprPtr& current : product) {
    IndexGroups<Container> currentGroups =
        get_unique_indices<Container>(current);

    for (Index& current : currentGroups.bra) {
      if (encounteredIndices.find(current) == encounteredIndices.end()) {
        encounteredIndices.insert(current);
        groups.bra.push_back(std::move(current));
      } else {
        detail::remove_one(groups.bra, current);
        detail::remove_one(groups.ket, current);
      }
    }

    // Same for ket indices
    for (Index& current : currentGroups.ket) {
      if (encounteredIndices.find(current) == encounteredIndices.end()) {
        encounteredIndices.insert(current);
        groups.ket.push_back(std::move(current));
      } else {
        detail::remove_one(groups.bra, current);
        detail::remove_one(groups.ket, current);
      }
    }
  }

  return groups;
}

/// Obtains the set of unique (non-repeated) indices used in the given
/// expression
template <typename Container>
IndexGroups<Container> get_unique_indices(const ExprPtr& expr) {
  if (expr.is<Constant>()) {
    return get_unique_indices<Container>(expr.as<Constant>());
  } else if (expr.is<Variable>()) {
    return get_unique_indices<Container>(expr.as<Variable>());
  } else if (expr.is<Tensor>()) {
    return get_unique_indices<Container>(expr.as<Tensor>());
  } else if (expr.is<Sum>()) {
    return get_unique_indices<Container>(expr.as<Sum>());
  } else if (expr.is<Product>()) {
    return get_unique_indices<Container>(expr.as<Product>());
  } else {
    throw std::runtime_error(
        "Encountered unsupported expression type in get_unique_indices");
  }
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_INDICES_HPP
