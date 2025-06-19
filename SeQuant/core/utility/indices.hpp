#ifndef SEQUANT_CORE_UTILITY_INDICES_HPP
#define SEQUANT_CORE_UTILITY_INDICES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <range/v3/view.hpp>

#include <algorithm>
#include <optional>
#include <set>
#include <type_traits>
#include <vector>

namespace sequant {

namespace detail {
template <typename Range>
struct not_in {
  const Range& range;

  not_in(const Range& range) : range(range) {}

  template <typename T>
  bool operator()(const T& element) const {
    return std::find(range.begin(), range.end(), element) == range.end();
  }
};

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
  Container aux;

  bool operator==(const IndexGroups<Container>& other) const {
    return bra == other.bra && ket == other.ket && aux == other.aux;
  }

  bool operator!=(const IndexGroups<Container>& other) const {
    return !(*this == other);
  }
};

/// A composite type for holding tensor-of-tensor indices
template <typename Container = std::vector<Index>>
struct TensorOfTensorIndices {
  Container outer;
  Container inner;
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

/// @returns Lists of non-contracted indices arising when contracting the two
/// given tensors in the order bra, ket, auxiliary
template <typename Container = container::svector<Index>>
IndexGroups<Container> get_uncontracted_indices(const Tensor& t1,
                                                const Tensor& t2) {
  static_assert(std::is_same_v<typename Container::value_type, Index>);

  IndexGroups<Container> groups;

  // Bra indices
  std::copy_if(t1.bra().begin(), t1.bra().end(), std::back_inserter(groups.bra),
               detail::not_in{t2.ket()});
  std::copy_if(t2.bra().begin(), t2.bra().end(), std::back_inserter(groups.bra),
               detail::not_in{t1.ket()});

  // Ket indices
  std::copy_if(t1.ket().begin(), t1.ket().end(), std::back_inserter(groups.ket),
               detail::not_in{t2.bra()});
  std::copy_if(t2.ket().begin(), t2.ket().end(), std::back_inserter(groups.ket),
               detail::not_in{t1.bra()});

  // Auxiliary indices
  std::copy_if(t1.aux().begin(), t1.aux().end(), std::back_inserter(groups.aux),
               detail::not_in{t2.aux()});
  std::copy_if(t2.aux().begin(), t2.aux().end(), std::back_inserter(groups.aux),
               detail::not_in{t1.aux()});

  return groups;
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

  for (const Index& current : tensor.aux()) {
    if (encounteredIndices.find(current) == encounteredIndices.end()) {
      groups.aux.push_back(current);
      encounteredIndices.insert(current);
    } else {
      detail::remove_one(groups.bra, current);
      detail::remove_one(groups.ket, current);
      detail::remove_one(groups.aux, current);
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
        detail::remove_one(groups.aux, current);
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
        detail::remove_one(groups.aux, current);
      }
    }

    // Same for aux indices
    for (Index& current : currentGroups.aux) {
      if (encounteredIndices.find(current) == encounteredIndices.end()) {
        encounteredIndices.insert(current);
        groups.aux.push_back(std::move(current));
      } else {
        detail::remove_one(groups.bra, current);
        detail::remove_one(groups.ket, current);
        detail::remove_one(groups.aux, current);
      }
    }
  }

  return groups;
}

/// Obtains the set of unique (non-repeated) indices used in the given
/// expression
template <typename Container>
IndexGroups<Container> get_unique_indices(const Expr& expr) {
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
template <typename Container>
IndexGroups<Container> get_unique_indices(const ExprPtr& expr) {
  return get_unique_indices<Container>(*expr);
}

template <typename Container = std::vector<Index>, typename Rng>
TensorOfTensorIndices<Container> tot_indices(Rng const& idxs) {
  using ranges::not_fn;
  using ranges::views::concat;
  using ranges::views::filter;
  using ranges::views::join;
  using ranges::views::transform;

  // Container indep_idxs;

  TensorOfTensorIndices<Container> result;
  auto& outer = result.outer;

  for (auto&& i : idxs | transform(&Index::proto_indices) | join)
    if (!ranges::contains(outer, i)) outer.emplace_back(i);

  for (auto&& i : idxs | filter(not_fn(&Index::has_proto_indices)))
    if (!ranges::contains(outer, i)) outer.emplace_back(i);

  auto& inner = result.inner;
  for (auto&& i : idxs | filter(&Index::has_proto_indices))
    inner.emplace_back(i);

  return result;
}

///
/// Does the numeric comparison of the index suffixes using less-than operator.
///
/// \param idx1
/// \param idx2
/// \return True if the numeric suffix of \c idx1 is less than that of \c idx2.
///
inline bool suffix_compare(Index const& idx1, Index const& idx2) {
  auto&& s1 = idx1.suffix();
  auto&& s2 = idx2.suffix();
  return (s1 && s2) && s1.value() < s2.value();
}

///
/// Given a range of Index returns a comma-separated string of
/// their full labels.
///   eg. [a_1^{i_1,i_2},a_2^{i_2,i_3}] -> "a_1i_1i_2,a_2i_2i_3"
///   eg. [i_1, i_2] -> "i_1,i_2"
///
template <typename Rng, typename Idx = ranges::range_value_t<Rng>,
          typename = std::enable_if_t<std::is_same_v<Idx, Index>>>
std::string csv_labels(Rng&& idxs) {
  using ranges::views::concat;
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::single;
  using ranges::views::transform;

  auto str = [](Index const& i) {
    auto v = concat(single(i.label()),                             //
                    i.proto_indices() | transform(&Index::label))  //
             | join;
    return sequant::to_string(v | ranges::to<std::wstring>);
  };

  return std::forward<Rng>(idxs)  //
         | transform(str)         //
         | intersperse(",")       //
         | join                   //
         | ranges::to<std::string>;
}

template <typename Container = container::svector<container::svector<Index>>>
Container external_indices(const Expr& expr) {
  if (!expr.is<Sum>() && !expr.is<Product>() && !expr.is<Tensor>()) {
    return {};
  }

  if (expr.is<Tensor>()) {
    const Tensor& tensor = expr.as<Tensor>();

    Container cont(std::max(tensor.bra_rank(),
                            std::max(tensor.ket_rank(), tensor.aux_rank())));
    for (std::size_t i = 0; i < tensor.ket_rank(); ++i) {
      cont.at(i).push_back(tensor.ket()[i]);
    }
    for (std::size_t i = 0; i < tensor.bra_rank(); ++i) {
      cont.at(i).push_back(tensor.bra()[i]);
    }
    for (std::size_t i = 0; i < tensor.aux_rank(); ++i) {
      cont.at(i).push_back(tensor.aux()[i]);
    }

    return cont;
  }

  std::optional<Tensor> symmetrizer;
  expr.visit(
      [&](const ExprPtr& expr) {
        if (expr.is<Tensor>() && (expr.as<Tensor>().label() == L"S" ||
                                  expr.as<Tensor>().label() == L"A")) {
          assert(!symmetrizer.has_value() ||
                 symmetrizer.value() == expr.as<Tensor>());
          symmetrizer = expr.as<Tensor>();
        }
      },
      true);

  if (symmetrizer.has_value()) {
    // Generate external index list from symmetrization operator
    return external_indices(symmetrizer.value());
  }

  IndexGroups groups = get_unique_indices<container::svector<Index>>(expr);

  Container cont(std::max(groups.bra.size(),
                          std::max(groups.ket.size(), groups.aux.size())));

  for (std::size_t i = 0; i < groups.ket.size(); ++i) {
    cont.at(i).push_back(groups.ket[i]);
  }
  for (std::size_t i = 0; i < groups.bra.size(); ++i) {
    cont.at(i).push_back(groups.bra[i]);
  }
  for (std::size_t i = 0; i < groups.aux.size(); ++i) {
    cont.at(i).push_back(groups.aux[i]);
  }

  return cont;
}

template <typename Container = container::svector<container::svector<Index>>>
Container external_indices(const ExprPtr& expr) {
  return external_indices<Container>(*expr);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_INDICES_HPP
