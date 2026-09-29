//
// Created by Eduard Valeyev on 2019-03-24.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include <compare>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

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

TensorCanonicalizer::index_comparer_t
TensorCanonicalizer::default_index_comparer() {
  return TensorIndexComparer{};
}

TensorCanonicalizer::index_pair_comparer_t
TensorCanonicalizer::default_index_pair_comparer() {
  return TensorIndexComparer{};
}

const std::shared_ptr<NullTensorCanonicalizer>&
NullTensorCanonicalizer::instance() {
  static const auto result = std::make_shared<NullTensorCanonicalizer>();
  return result;
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

namespace {
/// @return the canonicalization phase byproduct @p phase -- the convention
///         every TensorCanonicalizer::apply() returns, nullptr for +1 and
///         Constant(-1) for -1 -- multiplied by @p sign, spelled in the same
///         convention so that a +1 product stays nullptr
ExprPtr multiply_phase(ExprPtr phase, std::int8_t sign) {
  if (sign == 1) return phase;
  if (!phase) return ex<Constant>(-1);
  SEQUANT_ASSERT(phase->is<Constant>() && phase->as<Constant>().value() == -1);
  return {};
}
}  // namespace

bool braket_orientation_pinned(const AbstractTensor& t) {
  const auto lbl = t._label();
  return lbl == reserved::antisymm_label() || lbl == reserved::symm_label() ||
         lbl == reserved::transposition_label();
}

bool braket_foldable(const AbstractTensor& t) {
  return braket_swap_sign(t._braket_symmetry()).has_value() &&
         !braket_orientation_pinned(t);
}

std::int8_t DefaultTensorCanonicalizer::canonicalize_braket(AbstractTensor& t,
                                                            bool fold_signed) {
  if (!braket_foldable(t)) {
    return 1;
  }
  const auto bks = t._braket_symmetry();
  if (!fold_signed && bks == BraKetSymmetry::Antisymm) {
    return 1;
  }

  // the sign every respelling below contributes, for the caller to record
  std::int8_t sign = 1;

  // bra<->ket exchange is a symmetry for Symm/Antisymm tensors, so pick a
  // canonical orientation, a plain respelling. The choice is governed by the
  // canonical "colors" of the bra and ket bundles -- i.e. their index spaces,
  // not the index labels -- so the result is label-independent. Bundles with
  // identical spaces (e.g. g{p,q;r,s}) compare equal and are left untouched.
  const TensorBlockIndexComparer cmp;
  auto space_less = [&cmp](const Index& a, const Index& b) {
    return cmp.compare_spaces(a, b) < 0;
  };

  auto bra = mutable_bra_range(t);
  auto ket = mutable_ket_range(t);

  // Compare the bundles by their space sequences *sorted by color*, so the
  // decision is independent of the within-bundle index order. Column/perm
  // symmetry can permute the bra (and ket) order without changing the tensor,
  // and a comparison over the as-given order could otherwise pick different
  // orientations for equivalent inputs.
  std::vector<Index> bra_spaces, ket_spaces;
  for (auto&& idx : bra) bra_spaces.push_back(idx);
  for (auto&& idx : ket) ket_spaces.push_back(idx);

  ranges::sort(bra_spaces, space_less);
  ranges::sort(ket_spaces, space_less);

  // canonical orientation: the bundle whose spaces are lexicographically
  // larger goes to bra (three-way compare, so a full tie is detected without
  // re-comparing in reverse)
  const auto space_order = std::lexicographical_compare_three_way(
      bra_spaces.begin(), bra_spaces.end(), ket_spaces.begin(),
      ket_spaces.end(), [&cmp](const Index& a, const Index& b) {
        return cmp.compare_spaces(a, b) <=> 0;
      });
  const bool swap = space_order < 0;

  if (swap) {
    t._swap_bra_ket();
    // Symm: +1, Antisymm: -1
    sign = static_cast<std::int8_t>(sign * *braket_swap_sign(bks));
  }
  return sign;
}

ExprPtr DefaultTensorCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  const auto braket_sign = canonicalize_braket(t);

  const auto ctx = get_default_context_snapshot();
  auto result = this->apply(t, ctx.index_comparer(), ctx.index_pair_comparer());

  reset_tags(t);

  return multiply_phase(std::move(result), braket_sign);
}

ExprPtr TensorBlockCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  const auto braket_sign = canonicalize_braket(t, fold_signed_braket_);

  auto result = DefaultTensorCanonicalizer::apply(t, TensorBlockIndexComparer{},
                                                  TensorBlockIndexComparer{});

  reset_tags(t);

  return multiply_phase(std::move(result), braket_sign);
}

}  // namespace sequant
