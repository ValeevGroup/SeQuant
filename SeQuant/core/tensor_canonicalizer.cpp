//
// Created by Eduard Valeyev on 2019-03-24.
//

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <compare>
#include <cstdint>
#include <memory>
#include <mutex>
#include <type_traits>
#include <vector>

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
  static bool initialized_ = false;

  std::unique_lock lock(mtx_);
  if (!initialized_) {
    // Ensure DefaultTensorCanonicalizer is installed as the default
    // canonicalizer by default
    map_.emplace(L"", std::make_shared<DefaultTensorCanonicalizer>());
    initialized_ = true;
  }

  return std::make_pair(&map_, std::move(lock));
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
    SEQUANT_ASSERT(!has_duplicates(labels) &&
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

std::int8_t kramers_conjugate_mark(AbstractTensor& t) {
  auto* e = dynamic_cast<Expr*>(&t);
  if (!e || !dynamic_cast<Tensor*>(&t))
    throw std::logic_error(
        "kramers_conjugate_mark: only a Tensor carries the core conjugation "
        "states");
  // T{a;b} -> T⁺{b;a} = conj T{a;b}: the adjoint state with the bundles
  // exchanged is the core model's spelling of the elementwise conjugate. A
  // definite hermiticity consumes the state into the exchange itself
  // (T{b;a} = s conj T{a;b}) at the returned sign.
  return e->adjoint();
}

bool kramers_orientation_free(const AbstractTensor& t) {
  if (!t._is_cnumber() || braket_orientation_pinned(t)) return false;
  const auto bks = t._braket_symmetry();
  return braket_swap_sign(bks).has_value() ||
         bks == BraKetSymmetry::Conjugate ||
         bks == BraKetSymmetry::AntiConjugate;
}

std::pair<bool, std::int8_t> kramers_uprow_exchange(AbstractTensor& t) {
  if (!kramers_foldable(t)) return {false, 1};
  const auto bks = t._braket_symmetry();
  if (bks != BraKetSymmetry::Conjugate && bks != BraKetSymmetry::AntiConjugate)
    return {false, 1};
  const auto isr = get_default_context().index_space_registry();
  if (!isr) return {false, 1};
  auto n_down = [&isr](auto slots) {
    std::size_t n = 0;
    for (const Index& idx : slots)
      if (isr->kramers_partner(idx.space()) &&
          !isr->kramers_canonical(idx.space()))
        ++n;
    return n;
  };
  if (n_down(t._bra()) <= n_down(t._ket())) return {false, 1};
  auto* e = dynamic_cast<Expr*>(&t);
  if (!e) return {false, 1};
  // the definite hermiticity consumes the state: a bundle exchange at this
  // sign, the spelling of s conj T{a;b}
  const auto sign = e->adjoint();
  SEQUANT_ASSERT(!dynamic_cast<Tensor&>(t).adjointed());
  return {true, sign};
}

namespace {
/// the bra and ket of @p t in its VALUE orientation: a kept adjoint state
/// spells conj T{a;b} as T⁺{b;a} (see kramers_conjugate_mark), so the
/// value's bundles are the spelled ones exchanged back
std::pair<container::svector<Index>, container::svector<Index>>
kramers_value_bundles(const AbstractTensor& t) {
  container::svector<Index> bra, ket;
  for (const Index& idx : t._bra()) bra.push_back(idx);
  for (const Index& idx : t._ket()) ket.push_back(idx);
  if (const auto* tensor = dynamic_cast<const Tensor*>(&t);
      tensor && tensor->adjointed())
    std::swap(bra, ket);
  return {std::move(bra), std::move(ket)};
}
}  // namespace

bool kramers_foldable(const AbstractTensor& t) {
  return t._is_cnumber() &&
         t._kramers_symmetry() == KramersSymmetry::TimeReversal &&
         !braket_orientation_pinned(t);
}

bool kramers_flip_slots_deep(AbstractTensor& t) {
  const auto isr = get_default_context().index_space_registry();
  if (!isr) return false;
  bool flipped = false;
  auto flip = [&](auto&& slots) {
    for (auto& idx : slots) {
      auto f = kramers_flipped_deep(idx, *isr);
      if (f == idx) continue;
      const bool tagged = idx.tag().has_value();
      idx = std::move(f);
      if (tagged) idx.tag().assign(0);
      flipped = true;
    }
  };
  flip(t._bra_mutable());
  flip(t._ket_mutable());
  flip(t._aux_mutable());
  return flipped;
}

bool kramers_flip_slots(AbstractTensor& t) {
  const auto isr = get_default_context().index_space_registry();
  if (!isr) return false;
  bool flipped = false;
  // in place (not via _transform_indices: the block canonicalizer tags every
  // slot and Index::transform skips tagged indices); a tag present on the
  // slot is carried over
  auto flip = [&](auto&& slots) {
    for (auto& idx : slots) {
      auto f = kramers_flipped(idx, *isr);
      if (!f) continue;
      const bool tagged = idx.tag().has_value();
      idx = std::move(*f);
      if (tagged) idx.tag().assign(0);
      flipped = true;
    }
  };
  flip(t._bra_mutable());
  flip(t._ket_mutable());
  flip(t._aux_mutable());
  return flipped;
}

bool kramers_union_index(const Index& idx, const IndexSpaceRegistry& isr) {
  const auto& sp = idx.space();
  if (isr.kramers_partner(sp)) return false;  // a flavoured index
  // sp is the union of a Kramers-partnered pair of the same type: its quantum
  // numbers are EXACTLY the union of the pair's (the pair differs from sp in
  // the spin sector alone). A spin-free space of the same type that carries a
  // trait bit no partnered pair has (an AO/PAO-like space without flavoured
  // clones of its own) is not a union, even though a flavoured space's
  // quantum numbers are a subset of its own.
  for (const auto& s : isr) {
    if (s.type() != sp.type()) continue;
    const auto partner = isr.kramers_partner(s);
    if (!partner) continue;
    if ((s.qns() | partner->qns()) == sp.qns()) return true;
  }
  return false;
}

bool has_kramers_union_slot(const AbstractTensor& t) {
  const auto isr = get_default_context().index_space_registry();
  if (!isr) return false;
  auto any_union = [&](auto slots) {
    for (const Index& idx : slots)
      if (kramers_union_index(idx, *isr)) return true;
    return false;
  };
  return any_union(t._bra()) || any_union(t._ket()) || any_union(t._aux());
}

std::wstring kramers_flavor_key(const AbstractTensor& t, bool flipped) {
  const auto isr = get_default_context().index_space_registry();
  auto bundle = [&](auto slots) {
    std::wstring b;
    for (const Index& idx : slots) {
      if (!isr || !isr->kramers_partner(idx.space())) {
        b += L'-';
        continue;
      }
      const bool down = !isr->kramers_canonical(idx.space());
      b += (down != flipped) ? L'b' : L'a';  // up 'a' orders before down 'b'
    }
    std::sort(b.begin(), b.end());
    return b;
  };
  // the value orientation (a kept adjoint state exchanges the bundles), and
  // for a tensor whose bra<->ket exchange the fold may spell (a bare swap or
  // the adjoint itself) the two bundles in a canonical order, so the key is
  // orientation-invariant
  const auto [bra_v, ket_v] = kramers_value_bundles(t);
  std::wstring bra = bundle(bra_v), ket = bundle(ket_v), aux = bundle(t._aux());
  if (kramers_orientation_free(t) && ket < bra) std::swap(bra, ket);
  std::wstring key(t._label());
  key += L'|';
  key += bra;
  key += L'|';
  key += ket;
  key += L'|';
  key += aux;
  return key;
}

bool kramers_noncanonical(const AbstractTensor& t) {
  const auto isr = get_default_context().index_space_registry();
  if (!isr) return false;
  std::size_t n_up = 0, n_down = 0;
  auto count = [&](auto slots) {
    for (const Index& idx : slots) {
      if (!isr->kramers_partner(idx.space())) continue;
      if (isr->kramers_canonical(idx.space()))
        ++n_up;
      else
        ++n_down;
    }
  };
  count(t._bra());
  count(t._ket());
  count(t._aux());
  if (n_up + n_down == 0) return false;
  if (n_down != n_up) return n_down > n_up;
  return kramers_flavor_key(t, true) < kramers_flavor_key(t, false);
}

int canonicalize_kramers(AbstractTensor& t, bool mark) {
  if (!kramers_foldable(t)) return 1;
  const auto isr = get_default_context().index_space_registry();
  if (!isr) return 1;
  // orientation/permutation-invariant decision (see kramers_noncanonical);
  // count the down slots for the phase
  if (!kramers_noncanonical(t)) return 1;  // nothing to fold
  int n_down = 0;
  auto visit = [&](const Index& idx) {
    if (isr->kramers_partner(idx.space()) &&
        !isr->kramers_canonical(idx.space()))
      ++n_down;
  };
  for (const auto& idx : t._bra()) visit(idx);
  for (const auto& idx : t._ket()) visit(idx);
  for (const auto& idx : t._aux()) visit(idx);
  kramers_flip_slots(t);
  int phase = (n_down % 2) ? -1 : 1;
  if (mark) phase *= kramers_conjugate_mark(t);
  return phase;
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
  bool swap = space_order < 0;

  // Kramers (time-reversal) tensors: prefer the orientation whose bra
  // carries fewer down-flavored indices, so the up-row spelling is reached
  // by the braket move and "first flavored slot up" is a braket-invariant
  // notion for the Kramers fold. Ties fall through to the space criterion.
  if (t._kramers_symmetry() == KramersSymmetry::TimeReversal) {
    if (const auto isr = get_default_context().index_space_registry()) {
      auto n_down = [&isr](const std::vector<Index>& v) {
        std::size_t n = 0;
        for (const auto& idx : v)
          if (isr->kramers_partner(idx.space()) &&
              !isr->kramers_canonical(idx.space()))
            ++n;
        return n;
      };
      const auto nb = n_down(bra_spaces), nk = n_down(ket_spaces);
      if (nb != nk) swap = nb > nk;
    }
  }

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

  auto result =
      this->apply(t, this->index_comparer_, this->index_pair_comparer_);

  reset_tags(t);

  return multiply_phase(std::move(result), braket_sign);
}

ExprPtr TensorBlockCanonicalizer::apply(AbstractTensor& t) const {
  tag_indices(t);

  std::int8_t braket_sign = canonicalize_braket(t, fold_signed_braket_);
  const int kramers_phase = fold_kramers_ ? canonicalize_kramers(t) : 1;
  // the flipped spelling may prefer the other braket orientation
  if (fold_kramers_)
    braket_sign = static_cast<std::int8_t>(
        braket_sign * canonicalize_braket(t, fold_signed_braket_));

  auto result = DefaultTensorCanonicalizer::apply(t, TensorBlockIndexComparer{},
                                                  TensorBlockIndexComparer{});

  reset_tags(t);

  // combine the braket respelling sign with the Kramers fold phase
  return multiply_phase(std::move(result),
                        static_cast<std::int8_t>(braket_sign * kramers_phase));
}

}  // namespace sequant
