//
// Created by Eduard Valeyev on 2019-03-24.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include <algorithm>
#include <regex>
#include <type_traits>
#include <vector>

#include <range/v3/algorithm/adjacent_find.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/lexicographical_compare.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/range/access.hpp>

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
  // The map is seeded with DefaultTensorCanonicalizer as the default default
  // (label L""), so a bare Tensor canonicalizes (including the
  // braket-orientation fold, now part of DefaultTensorCanonicalizer::apply)
  // even when no canonicalizer was registered explicitly. Explicit
  // register_instance calls override the seed as before.
  static container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>
      map_ = [] {
        container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>> m;
        m.emplace(L"", std::make_shared<DefaultTensorCanonicalizer>());
        return m;
      }();
  static std::recursive_mutex mtx_;
  return std::make_pair(&map_, std::unique_lock<std::recursive_mutex>{mtx_});
}

container::vector<std::wstring>&
TensorCanonicalizer::default_cardinal_tensor_labels_accessor() {
  // {antisymm_label, symm_label, transposition_label} is the default
  static container::vector<std::wstring> default_ctlabels_{
      reserved::antisymm_label(), reserved::symm_label(),
      reserved::transposition_label()};
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
  if constexpr (assert_enabled()) {
    // check for duplicates within user provided labels
    auto sorted_labels = labels;
    ranges::sort(sorted_labels);
    [[maybe_unused]] auto duplicate = ranges::adjacent_find(sorted_labels);
    SEQUANT_ASSERT(duplicate == sorted_labels.end() &&
                   "cardinal tensor labels must not contain duplicates");

    // check if any label conflicts with existing ones
    const auto& existing = cardinal_tensor_labels_accessor();
    for (const auto& label : labels) {
      [[maybe_unused]] auto conflict = ranges::find(existing, label);
      SEQUANT_ASSERT(conflict == existing.end() &&
                     "cardinal tensor labels must not contain duplicates");
    }
  }
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
    throw Exception(
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

bool braket_orientation_pinned(const AbstractTensor& t) {
  const auto lbl = t._label();
  return lbl == reserved::antisymm_label() || lbl == reserved::symm_label() ||
         lbl == reserved::transposition_label();
}

bool prefer_swapped_braket(const AbstractTensor& t) {
  const TensorBlockIndexComparer space_cmp;
  auto space_less = [&space_cmp](const Index& a, const Index& b) {
    return space_cmp.compare_spaces(a, b) < 0;
  };
  auto sorted = [](auto&& rng, auto&& less) {
    std::vector<Index> v;
    for (const auto& idx : rng) v.push_back(idx);
    ranges::sort(v, less);
    return v;
  };

  // Space level: the space-lexicographically larger bundle belongs in the
  // bra (the historical convention: e.g. the half-tensor X{;a;x} folds into
  // X{a;;x}).
  const auto bra_by_space = sorted(t._bra(), space_less);
  const auto ket_by_space = sorted(t._ket(), space_less);
  if (ranges::lexicographical_compare(bra_by_space, ket_by_space, space_less))
    return true;
  if (ranges::lexicographical_compare(ket_by_space, bra_by_space, space_less))
    return false;

  // Full space tie: break it on the index labels, keeping the
  // label-lexicographically SMALLER bundle in the bra, so label-ascending
  // spellings (e.g. g{p1,p2;p3,p4}) remain canonical as written. Identical
  // bundles (diagonal trace T{p,q;p,q}) compare equal and never swap.
  const auto bra_full = sorted(t._bra(), std::less<Index>{});
  const auto ket_full = sorted(t._ket(), std::less<Index>{});
  return ranges::lexicographical_compare(ket_full, bra_full);
}

namespace {

/// applies the canonical braket orientation (prefer_swapped_braket) to a
/// braket-foldable tensor: Symm braket swaps freely, Conjugate braket swaps
/// with the elementwise-conjugation marker toggled (T{q;p} = conj(T{p;q})).
/// Operator-valued tensors (swap exchanges creators/annihilators) and the
/// reserved bookkeeping operators (orientation defines/extracts external
/// indices) are left untouched.
/// @return true if bra and ket were swapped
bool apply_canonical_braket_orientation(AbstractTensor& t) {
  const auto bks = t._braket_symmetry();
  const bool foldable =
      (bks == BraKetSymmetry::Symm || bks == BraKetSymmetry::Conjugate) &&
      t._is_cnumber() && !braket_orientation_pinned(t);
  if (!foldable || !prefer_swapped_braket(t)) return false;
  t._swap_bra_ket();
  if (bks == BraKetSymmetry::Conjugate) t._conjugate();
  return true;
}

}  // namespace

ExprPtr DefaultTensorCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  // pick the canonical braket orientation of braket-foldable tensors (same
  // fold as TensorBlockCanonicalizer::apply and
  // TensorNetworkV3::canonicalize_graph): a bare tensor's canonicalization
  // must spell one value one way regardless of the route it took
  apply_canonical_braket_orientation(t);

  auto result =
      this->apply(t, this->index_comparer_, this->index_pair_comparer_);

  reset_tags(t);

  return result;
}

template <typename Callable, typename... Args>
using suitable_call_operator =
    decltype(std::declval<Callable>()(std::declval<Args>()...));

bool TensorBlockCanonicalizer::orient_braket_by_color(AbstractTensor& t) const {
  // The canonical orientation is a pure function of the tensor's content
  // (prefer_swapped_braket: sorted-bundle lexicographic comparison over full
  // Index values -- space, then label), shared with the bundle swap in
  // TensorNetworkV3::canonicalize_graph so the per-tensor and network
  // canonicalization routes spell one value one way. E.g. a half-tensor
  // X{;a;x} folds into X{a;;x}, and same-space bundles tie-break on the
  // (final, for a lone tensor: all-named) index labels.
  if (prefer_swapped_braket(t)) {
    t._swap_bra_ket();
    return true;
  }
  return false;
}

bool TensorBlockCanonicalizer::fold_conjugate_braket(AbstractTensor& t) const {
  // For a Conjugate tensor a bra<->ket swap is not free -- by the symmetry
  // relation T{q;p} = conj(T{p;q}) the swapped spelling represents the
  // conjugate value -- so the reorientation toggles the tensor's own
  // elementwise-conjugation marker (_conjugate()), keeping the represented
  // value invariant. Same color rule as the Symm fold; also reports whether
  // it swapped (legacy byproduct channel, to be retired once every consumer
  // reads the marker off the tensor).
  if (t._braket_symmetry() != BraKetSymmetry::Conjugate) return false;
  // value identity only holds for c-number tensors (an operator-valued
  // tensor's bra<->ket swap exchanges creators and annihilators)
  if (!t._is_cnumber()) return false;
  // reserved bookkeeping operators must keep their orientation
  if (braket_orientation_pinned(t)) return false;
  const bool swapped = orient_braket_by_color(t);
  if (swapped) t._conjugate();
  return swapped;
}

ExprPtr TensorBlockCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  // pick the canonical braket orientation (shared with
  // DefaultTensorCanonicalizer::apply and
  // TensorNetworkV3::canonicalize_graph)
  apply_canonical_braket_orientation(t);

  auto result = DefaultTensorCanonicalizer::apply(t, TensorBlockIndexComparer{},
                                                  TensorBlockIndexComparer{});

  reset_tags(t);

  return result;
}

}  // namespace sequant
