//
// Created by Eduard Valeyev on 2019-03-24.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include <regex>
#include <type_traits>

#include <range/v3/functional/identity.hpp>

namespace sequant {

template <typename T>
using get_support = decltype(std::get<0>(std::declval<T>()));

template <typename T>
constexpr bool is_tuple_like_v = meta::is_detected_v<get_support, T>;

struct TensorBlockIndexComparer {
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const {
    return compare<T>(lhs, rhs) < 0;
  }

  template <typename T>
  int compare(const T& lhs, const T& rhs) const {
    if constexpr (is_tuple_like_v<T>) {
      static_assert(
          std::tuple_size_v<T> == 2,
          "TensorBlockIndexComparer can only deal with tuple-like objects "
          "of size 2");
      const auto& lhs_first = std::get<0>(lhs);
      const auto& lhs_second = std::get<1>(lhs);
      const auto& rhs_first = std::get<0>(rhs);
      const auto& rhs_second = std::get<1>(rhs);

      static_assert(std::is_same_v<std::decay_t<decltype(lhs_first)>, Index>,
                    "TensorBlockIndexComparer can only work with indices");
      static_assert(std::is_same_v<std::decay_t<decltype(lhs_second)>, Index>,
                    "TensorBlockIndexComparer can only work with indices");
      static_assert(std::is_same_v<std::decay_t<decltype(rhs_first)>, Index>,
                    "TensorBlockIndexComparer can only work with indices");
      static_assert(std::is_same_v<std::decay_t<decltype(rhs_second)>, Index>,
                    "TensorBlockIndexComparer can only work with indices");

      // First compare only index spaces of equivalent pairs
      int res = compare_spaces(lhs_first, rhs_first);
      if (res != 0) {
        return res;
      }

      res = compare_spaces(lhs_second, rhs_second);
      if (res != 0) {
        return res;
      }

      // Then consider tags of equivalent pairs
      res = compare_tags(lhs_first, rhs_first);
      if (res != 0) {
        return res;
      }

      res = compare_tags(lhs_second, rhs_second);
      return res;
    } else {
      static_assert(std::is_same_v<std::decay_t<T>, Index>,
                    "TensorBlockIndexComparer can only work with indices");

      int res = compare_spaces(lhs, rhs);
      if (res != 0) {
        return res;
      }

      res = compare_tags(lhs, rhs);
      return res;
    }
  }

  int compare_spaces(const Index& lhs, const Index& rhs) const {
    if (lhs.space() != rhs.space()) {
      return lhs.space() < rhs.space() ? -1 : 1;
    }

    if (lhs.has_proto_indices() != rhs.has_proto_indices()) {
      return lhs.has_proto_indices() ? -1 : 1;
    }

    if (lhs.proto_indices().size() != rhs.proto_indices().size()) {
      return lhs.proto_indices().size() < rhs.proto_indices().size() ? -1 : 1;
    }

    for (std::size_t i = 0; i < lhs.proto_indices().size(); ++i) {
      const auto& lhs_proto = lhs.proto_indices()[i];
      const auto& rhs_proto = rhs.proto_indices()[i];

      int res = compare_spaces(lhs_proto, rhs_proto);
      if (res != 0) {
        return res;
      }
    }

    // Index spaces are equal
    return 0;
  }

  int compare_tags(const Index& lhs, const Index& rhs) const {
    if (!lhs.tag().has_value() || !rhs.tag().has_value()) {
      // We only compare tags if both indices have a tag
      return 0;
    }

    const int lhs_tag = lhs.tag().value<int>();
    const int rhs_tag = rhs.tag().value<int>();

    if (lhs_tag != rhs_tag) {
      return lhs_tag < rhs_tag ? -1 : 1;
    }

    return 0;
  }
};

struct TensorIndexComparer {
  template <typename T>
  bool operator()(const T& lhs, const T& rhs) const {
    TensorBlockIndexComparer block_comp;

    int res = block_comp.compare<T>(lhs, rhs);

    if (res != 0) {
      return res < 0;
    }

    // Fall back to regular index compare to break the tie
    if constexpr (is_tuple_like_v<T>) {
      static_assert(std::tuple_size_v<T> == 2,
                    "TensorIndexComparer can only deal with tuple-like objects "
                    "of size 2");

      const Index& lhs_first = std::get<0>(lhs);
      const Index& lhs_second = std::get<1>(lhs);
      const Index& rhs_first = std::get<0>(rhs);
      const Index& rhs_second = std::get<1>(rhs);

      if (lhs_first != rhs_first) {
        return lhs_first < rhs_first;
      }

      return lhs_second < rhs_second;
    } else {
      return lhs < rhs;
    }
  }
};

TensorCanonicalizer::~TensorCanonicalizer() = default;

std::pair<container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>*,
          std::unique_lock<std::recursive_mutex>>
TensorCanonicalizer::instance_map_accessor() {
  static container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>
      map_;
  static std::recursive_mutex mtx_;
  return std::make_pair(&map_, std::unique_lock<std::recursive_mutex>{mtx_});
}

container::vector<std::wstring>&
TensorCanonicalizer::default_cardinal_tensor_labels_accessor() {
  // {antisymm_label, symm_label} is the default
  static container::vector<std::wstring> default_ctlabels_{
      reserved::antisymm_label(), reserved::symm_label()};
  return default_ctlabels_;
}

container::vector<std::wstring>&
TensorCanonicalizer::cardinal_tensor_labels_accessor() {
  static container::vector<std::wstring> ctlabels_ =
      default_cardinal_tensor_labels_accessor();
  return ctlabels_;
}

void TensorCanonicalizer::set_cardinal_tensor_labels(
    const container::vector<std::wstring>& labels) {
  // check for duplicates
#ifdef SEQUANT_ASSERT_ENABLED
  // check for duplicates within user provided labels
  auto sorted_labels = labels;
  ranges::sort(sorted_labels);
  auto duplicate = ranges::adjacent_find(sorted_labels);
  SEQUANT_ASSERT(duplicate == sorted_labels.end() &&
                 "cardinal tensor labels must not contain duplicates");

  // check if any label conflicts with existing ones
  const auto& existing = cardinal_tensor_labels_accessor();
  for (const auto& label : labels) {
    auto conflict = ranges::find(existing, label);
    SEQUANT_ASSERT(conflict == existing.end() &&
                   "cardinal tensor labels must not contain duplicates");
  }
#endif
  auto& ctlabels = cardinal_tensor_labels_accessor();
  // get defaults
  ctlabels = default_cardinal_tensor_labels_accessor();
  // append
  ctlabels.insert(ctlabels.end(), labels.begin(), labels.end());
}

void TensorCanonicalizer::reset_cardinal_tensor_labels() {
  cardinal_tensor_labels_accessor() = default_cardinal_tensor_labels_accessor();
}

void TensorCanonicalizer::clear_all_cardinal_tensor_labels() {
  cardinal_tensor_labels_accessor().clear();
}

std::shared_ptr<TensorCanonicalizer>
TensorCanonicalizer::nondefault_instance_ptr(std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  // look for label-specific canonicalizer
  auto it = map_ptr->find(std::wstring{label});
  if (it != map_ptr->end()) {
    return it->second;
  } else
    return {};
}

std::shared_ptr<TensorCanonicalizer> TensorCanonicalizer::instance_ptr(
    std::wstring_view label) {
  auto result = nondefault_instance_ptr(label);
  if (!result)  // not found? look for default
    result = nondefault_instance_ptr(L"");
  return result;
}

std::shared_ptr<TensorCanonicalizer> TensorCanonicalizer::instance(
    std::wstring_view label) {
  auto inst_ptr = instance_ptr(label);
  if (!inst_ptr)
    throw std::runtime_error(
        "must first register canonicalizer via "
        "TensorCanonicalizer::register_instance(...)");
  return inst_ptr;
}

void TensorCanonicalizer::register_instance(
    std::shared_ptr<TensorCanonicalizer> can, std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  (*map_ptr)[std::wstring{label}] = can;
}

bool TensorCanonicalizer::try_register_instance(
    std::shared_ptr<TensorCanonicalizer> can, std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  if (!map_ptr->contains(std::wstring{label})) {
    (*map_ptr)[std::wstring{label}] = can;
    return true;
  } else
    return false;
}

void TensorCanonicalizer::deregister_instance(std::wstring_view label) {
  auto&& [map_ptr, lock] = instance_map_accessor();
  auto it = map_ptr->find(std::wstring{label});
  if (it != map_ptr->end()) {
    map_ptr->erase(it);
  }
}

TensorCanonicalizer::index_comparer_t TensorCanonicalizer::index_comparer_ =
    TensorIndexComparer{};

TensorCanonicalizer::index_pair_comparer_t
    TensorCanonicalizer::index_pair_comparer_ = TensorIndexComparer{};

const TensorCanonicalizer::index_comparer_t&
TensorCanonicalizer::index_comparer() {
  return index_comparer_;
}

void TensorCanonicalizer::index_comparer(index_comparer_t comparer) {
  index_comparer_ = std::move(comparer);
}

const TensorCanonicalizer::index_pair_comparer_t&
TensorCanonicalizer::index_pair_comparer() {
  return index_pair_comparer_;
}

void TensorCanonicalizer::index_pair_comparer(index_pair_comparer_t comparer) {
  index_pair_comparer_ = std::move(comparer);
}

ExprPtr NullTensorCanonicalizer::apply(AbstractTensor&) const { return {}; }

void DefaultTensorCanonicalizer::tag_indices(AbstractTensor& t) const {
  // tag all indices as ext->true/ind->false
  ranges::for_each(slots(t), [this](auto& idx) {
    auto it = external_indices_.find(idx);
    auto is_ext = it != external_indices_.end();
    idx.tag().assign(
        is_ext ? 0 : 1);  // ext -> 0, int -> 1, so ext will come before
  });
}

ExprPtr DefaultTensorCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  auto result =
      this->apply(t, this->index_comparer_, this->index_pair_comparer_);

  reset_tags(t);

  return result;
}

template <typename Callable, typename... Args>
using suitable_call_operator =
    decltype(std::declval<Callable>()(std::declval<Args>()...));

ExprPtr TensorBlockCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  auto result = DefaultTensorCanonicalizer::apply(t, TensorBlockIndexComparer{},
                                                  TensorBlockIndexComparer{});

  reset_tags(t);

  return result;
}

}  // namespace sequant
