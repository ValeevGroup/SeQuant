//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT_EXPRESSIONS_TENSOR_HPP
#define SEQUANT_EXPRESSIONS_TENSOR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/labeled.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/utility/strong.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

namespace sequant {

// strong type wrapper for objects associated with bra
DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(bra);
// strong type wrapper for objects associated with ket
DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(ket);
// strong type wrapper for objects associated with aux
DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(aux);

/// @brief a Tensor is an instance of AbstractTensor over a scalar field, i.e.
/// Tensors have commutative addition and product operations
class Tensor : public Expr, public AbstractTensor, public MutatableLabeled {
 private:
  using index_container_type = container::svector<Index>;
  static auto make_indices(IndexList indices) { return indices; }
  static auto make_indices(WstrList index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label : index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
  static auto make_indices(
      std::initializer_list<const wchar_t *> index_labels) {
    index_container_type result;
    result.reserve(index_labels.size());
    for (const auto &label : index_labels) {
      result.push_back(Index{label});
    }
    return result;
  }
  template <typename IndexRange>
  static index_container_type make_indices(IndexRange &&indices) {
    if constexpr (std::is_same_v<index_container_type,
                                 std::decay_t<IndexRange>>) {
      return std::forward<IndexRange>(indices);
    } else {
      using ranges::begin;
      using ranges::end;
      return index_container_type(begin(indices), end(indices));
    }
  }

  /// @return concatenated view of the bra and ket index ranges
  auto braket() { return ranges::views::concat(bra_, ket_); }

  /// @return concatenated view of bra, ket, and aux index ranges
  auto indices() { return ranges::views::concat(bra_, ket_, aux_); }

  /// asserts that @p label is not reserved
  /// @note Tensor with reserved labels are constructed using friends of Tensor
  ///       or are to be avoided to avoid confusion with other tensorials like
  ///       NormalOperator
  /// @param label a Tensor label candidate
  void assert_nonreserved_label(std::wstring_view label) const;

  // utility for dispatching to private ctor
  struct reserved_tag {};

  // list of friends who can make Tensor objects with reserved labels
  friend ExprPtr make_overlap(const Index &bra_index, const Index &ket_index);

  /// validates bra, ket, and aux indices

  // clang-format off
  /// @throw std::invalid_argument if `NDEBUG` not `#define`d and:
  ///        - `symmetry()==Symmetry::Antisymm` and `bra()` or `ket()` contains null indices, or
  ///        - `aux()` contains null indices, or
  ///        - there are duplicate indices within bra, within ket, or within aux
  // clang-format on
  void validate_indices() const {
#ifndef NDEBUG
    // antisymmetric bra or ket cannot support null indices (because it is not
    // clear what antisymmetry means if some slots can be empty; permutation of
    // 2 empty slots is supposed to do what?)
    // by analogy symmetric bra or ket should not have null indices, but limited
    // circumstances do allow null indices in such context ... but no use cases
    if (symmetry() != Symmetry::Nonsymm) {
      if (!bra_.empty() && ranges::contains(bra_, Index::null))
        throw std::invalid_argument(
            "Tensor ctor: found null indices in symmetric/antisymmetric bra");
      if (!ket_.empty() && ranges::contains(ket_, Index::null))
        throw std::invalid_argument(
            "Tensor ctor: found null indices in symmetric/antisymmetric ket");
    } else {  // asymmetric tensor

      // matching bra and ket slots cannot be both empty
      const auto braket_rank = std::min(bra_rank(), ket_rank());
      for (std::size_t r = 0; r != braket_rank; ++r) {
        if (bra_[r] == Index::null && ket_[r] == Index::null)
          throw std::invalid_argument(
              "Tensor ctor: found null indices in both matching slots of "
              "asymmetric bra and ket");
      }

      // if bra and ket differ in size, make sure unpaired slots are not empty
      if (bra_rank() != ket_rank()) {
        const auto longer_bundle_type =
            bra_rank() > ket_rank() ? SlotType::Bra : SlotType::Ket;
        auto *longer_bundle = longer_bundle_type == SlotType::Bra
                                  ? &bra_[0]
                                  : &ket_[0];  // n.b. these are contiguous
        const auto rank = std::max(bra_rank(), ket_rank());

        for (std::size_t r = braket_rank; r != rank; ++r) {
          if (longer_bundle[r] == Index::null)
            throw std::invalid_argument(
                (std::string("Tensor ctor: found null index in a slot of "
                             "asymmetric ") +
                 (longer_bundle_type == SlotType::Bra ? "bra" : "ket") +
                 " that does have a matching " +
                 (longer_bundle_type == SlotType::Bra ? "ket" : "bra") +
                 " slot")
                    .c_str());
        }
      }
    }
    if (!aux_.empty() && ranges::contains(aux_, Index::null)) {
      throw std::invalid_argument("Tensor ctor: found null aux indices");
    }
    // check for duplicates
    {
      auto throw_if_contains_duplicates = [](const auto &indices,
                                             const char *id) -> void {
        // sort via ptrs
        auto index_ptrs = indices |
                          ranges::views::transform(
                              [&](const Index &index) { return &index; }) |
                          ranges::to<container::svector<Index const *>>;
        ranges::sort(index_ptrs,
                     [](Index const *l, Index const *r) { return *l < *r; });
        if (ranges::adjacent_find(index_ptrs,
                                  [](Index const *l, Index const *r) {
                                    // N.B. multiple null indices OK
                                    return *l == *r && *l != Index::null;
                                  }) != index_ptrs.end())
          throw std::invalid_argument(
              (std::string("Tensor ctor: duplicate ") + id + " indices")
                  .c_str());
      };
      throw_if_contains_duplicates(bra_, "bra");
      throw_if_contains_duplicates(ket_, "ket");
      throw_if_contains_duplicates(aux_, "aux");
    }
#endif
  }

  /// put slots/slot bundles in canonical order:
  /// - if tensor is asymmetric but particle symmetric make sure
  ///   braket bundles are first, then bra-only braket bundles (paired with
  ///   empty ket slots or unpaired if there are no unpaired ket slots) then, if
  ///   any, ket-only braket bundles (unpaired).
  void canonicalize_slots() {
    // if tensor is particle symmetric make sure
    // braket bundles are first, then bra-only, then ket-only
    if (symmetry() == Symmetry::Nonsymm &&
        column_symmetry() == ColumnSymmetry::Symm) {
      const bool have_empty_slots =
          bra_rank() != bra_net_rank() || ket_rank() != ket_net_rank();
      if (have_empty_slots) {
        decltype(bra_)::value_type canonical_bra_;
        canonical_bra_.reserve(bra_.size());
        decltype(ket_)::value_type canonical_ket_;
        canonical_ket_.reserve(ket_.size());
        // push all braket bundles first
        for (std::size_t p = 0; p != std::min(bra_rank(), ket_rank()); ++p) {
          const auto nonempty_bra = bra_[p].nonnull();
          const auto nonempty_ket = ket_[p].nonnull();
          assert(nonempty_bra ||
                 nonempty_ket);  // validate_indices() checked this
          if (nonempty_bra && nonempty_ket) {
            // we know that no empty braket bundles exist, so move the indices
            // out
            canonical_bra_.emplace_back(std::move(bra_[p]));
            canonical_ket_.push_back(std::move(ket_[p]));
          }
        }

        std::size_t num_bra_indices = canonical_bra_.size();
        std::size_t num_ket_indices = canonical_ket_.size();
        for (auto &&idx : bra_) {
          if (idx.nonnull()) {
            canonical_bra_.emplace_back(std::move(idx));
            ++num_bra_indices;
            canonical_ket_.emplace_back(Index::null);
          }
        }
        for (auto &&idx : ket_) {
          if (idx.nonnull()) {
            canonical_ket_.emplace_back(std::move(idx));
            ++num_ket_indices;
          }
        }
        // done, assert that all nonnull indices moved out
        assert(ranges::count_if(
                   bra_, [](const Index &idx) { return idx.nonnull(); }) == 0);
        ;
        assert(ranges::count_if(
                   ket_, [](const Index &idx) { return idx.nonnull(); }) == 0);
        ;
        bra_ = sequant::bra(std::move(canonical_bra_));
        ket_ = sequant::ket(std::move(canonical_ket_));
        bra_net_rank_ = num_bra_indices;
        ket_net_rank_ = num_ket_indices;
      }
    }
  }

  template <basic_string_convertible S, range_of_castables_to_index IndexRange1,
            range_of_castables_to_index IndexRange2,
            range_of_castables_to_index IndexRange3>
  Tensor(S &&label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices,
         const aux<IndexRange3> &aux_indices, reserved_tag,
         Symmetry s = Symmetry::Nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ColumnSymmetry ps = ColumnSymmetry::Symm)
      : label_(to_wstring(std::forward<S>(label))),
        bra_(make_indices(bra_indices)),
        ket_(make_indices(ket_indices)),
        aux_(make_indices(aux_indices)),
        symmetry_(s),
        braket_symmetry_(bks),
        column_symmetry_(ps),
        bra_net_rank_(ranges::count_if(
            bra_, [](const Index &idx) { return static_cast<bool>(idx); })),
        ket_net_rank_(ranges::count_if(
            ket_, [](const Index &idx) { return static_cast<bool>(idx); })) {
    validate_indices();
    validate_symmetries();
    canonicalize_slots();
  }

  template <basic_string_convertible S>
  Tensor(S &&label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices,
         aux<index_container_type> &&aux_indices, reserved_tag,
         Symmetry s = Symmetry::Nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ColumnSymmetry ps = ColumnSymmetry::Symm)
      : label_(to_wstring(std::forward<S>(label))),
        bra_(std::move(bra_indices)),
        ket_(std::move(ket_indices)),
        aux_(std::move(aux_indices)),
        symmetry_(s),
        braket_symmetry_(bks),
        column_symmetry_(ps),
        bra_net_rank_(ranges::count_if(
            bra_, [](const Index &idx) { return static_cast<bool>(idx); })),
        ket_net_rank_(ranges::count_if(
            ket_, [](const Index &idx) { return static_cast<bool>(idx); })) {
    validate_indices();
    validate_symmetries();
    canonicalize_slots();
  }

 public:
  /// constructs an uninitialized Tensor
  /// @sa Tensor::operator bool()
  Tensor() = default;
  virtual ~Tensor();

  // clang-format off
  /// @name nontrivial constructors
  /// @note if `NDEBUG` is not `#define`d and invalid combinations of indices are found, these throw `std::invalid_argument`. Specifically, these throw if
  /// - null indices are found in any slot of antisymmetric bra or ket (it's not clear what permutation of 2 empty slots is supposed to do in general)
  /// - null indices are found in any slot of symmetric bra or ket: can assign consistent semantics if have empty slots in the shorter of the bra/ket, but not clear if there is a use case that demands this
  /// - null indices are found in both slots of bra-ket slot pairs (meaning is unclear)
  /// - null indices are found in any aux slot (why would empty slots be needed?)
  /// - asymmetric tensor has a bra/ket slot without a matching ket/bra slot; this is primarily to make operations (such as canonicalization) easier on such tensors; null indices can be used to ensure complete bra/ket bundles in asymmetric tensors
  /// - duplicate indices are found within bra, within ket, or within aux (but same index can appear in bra and ket, for example).
  /// @{
  // clang-format on

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  template <basic_string_convertible S, range_of_castables_to_index IndexRange1,
            range_of_castables_to_index IndexRange2>
  Tensor(S &&label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices, Symmetry s = Symmetry::Nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ColumnSymmetry ps = ColumnSymmetry::Symm)
      : Tensor(std::forward<S>(label), bra_indices, ket_indices, sequant::aux{},
               reserved_tag{}, s, bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param aux_indices list of aux indices (or objects that can be
  /// converted to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  template <basic_string_convertible S, range_of_castables_to_index IndexRange1,
            range_of_castables_to_index IndexRange2,
            range_of_castables_to_index IndexRange3>
  Tensor(S &&label, const bra<IndexRange1> &bra_indices,
         const ket<IndexRange2> &ket_indices,
         const aux<IndexRange3> &aux_indices, Symmetry s = Symmetry::Nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ColumnSymmetry ps = ColumnSymmetry::Symm)
      : Tensor(std::forward<S>(label), bra_indices, ket_indices, aux_indices,
               reserved_tag{}, s, bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  template <basic_string_convertible S>
  Tensor(S &&label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices,
         Symmetry s = Symmetry::Nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ColumnSymmetry ps = ColumnSymmetry::Symm)
      : Tensor(std::forward<S>(label), std::move(bra_indices),
               std::move(ket_indices), sequant::aux{}, reserved_tag{}, s, bks,
               ps) {
    assert_nonreserved_label(label_);
  }

  /// @param label the tensor label
  /// @param bra_indices list of bra indices (or objects that can be converted
  /// to indices)
  /// @param ket_indices list of ket indices (or objects that can be converted
  /// to indices)
  /// @param aux_indices list of aux indices (or objects that can be
  /// converted to indices)
  /// @param s the symmetry of bra or ket
  /// @param bks the symmetry with respect to bra-ket exchange
  /// @param ps the symmetry under exchange of particles
  template <basic_string_convertible S>
  Tensor(S &&label, bra<index_container_type> &&bra_indices,
         ket<index_container_type> &&ket_indices,
         aux<index_container_type> &&aux_indices,
         Symmetry s = Symmetry::Nonsymm,
         BraKetSymmetry bks = get_default_context().braket_symmetry(),
         ColumnSymmetry ps = ColumnSymmetry::Symm)
      : Tensor(std::forward<S>(label), std::move(bra_indices),
               std::move(ket_indices), std::move(aux_indices), reserved_tag{},
               s, bks, ps) {
    assert_nonreserved_label(label_);
  }

  /// @}

  /// @return true if the Tensor is initialized
  explicit operator bool() const { return !label_.empty(); }

  /// @return "core" label of the tensor
  std::wstring_view label() const override { return label_; }
  void set_label(std::wstring label) override {
    label_ = std::move(label);
    reset_hash_value();
  }
  /// @return the bra slot range (empty slots are occupied by null indices)
  const auto &bra() const { return bra_; }
  void set_bra(index_container_type indices) {
    bra_ = sequant::bra(std::move(indices));
    reset_hash_value();
  }
  /// @return the ket slot range (empty slots are occupied by null indices)
  const auto &ket() const { return ket_; }
  void set_ket(index_container_type indices) {
    ket_ = sequant::ket(std::move(indices));
    reset_hash_value();
  }
  /// @return the aux slot range (empty slots are occupied by null indices)
  const auto &aux() const { return aux_; }
  void set_aux(index_container_type indices) {
    aux_ = sequant::aux(std::move(indices));
    reset_hash_value();
  }
  /// @return concatenated view of the bra and ket slot ranges
  auto braket() const { return ranges::views::concat(bra_, ket_); }
  /// @return concatenated view of the bra and ket index ranges (i.e., nonempty
  /// slots)
  auto braket_indices() const {
    return ranges::views::filter(
        braket(), [](const Index &idx) { return idx.nonnull(); });
  }
  /// @return concatenated view of all indices of this tensor (bra, ket and
  /// aux)
  auto braketaux() const { return ranges::views::concat(bra_, ket_, aux_); }
  /// @return concatenated view of all slots
  auto slots() const { return ranges::views::concat(bra_, ket_, aux_); }
  /// @return concatenated view of all nonnull indices of this tensor (bra, ket
  /// and aux)
  auto braketaux_indices() const {
    return ranges::views::filter(
        braketaux(), [](const Index &idx) { return idx.nonnull(); });
  }
  /// @return concatenated view of all nonnull indices of this tensor
  auto indices() const {
    return ranges::views::filter(
        slots(), [](const Index &idx) { return idx.nonnull(); });
  }
  /// @return view of the bra+ket slots
  /// @note this is to work around broken lookup rules
  auto const_braket() const { return this->braket(); }
  /// @return view of the bra+ket indices
  /// @note this is to work around broken lookup rules
  auto const_braket_indices() const { return this->braket_indices(); }
  /// @return const view of all slots
  /// @note this is to work around broken lookup rules
  auto const_braketaux() const { return this->braketaux(); }
  /// @return const view of all slots
  /// @note this is to work around broken lookup rules
  auto const_slots() const { return this->slots(); }
  /// @return const view of all indices
  /// @note this is to work around broken lookup rules
  auto const_braketaux_indices() const { return this->braketaux_indices(); }
  /// @return const view of all indices
  /// @note this is to work around broken lookup rules
  auto const_indices() const { return this->indices(); }
  /// Returns the Symmetry object describing the symmetry of the bra and ket of
  /// the Tensor, i.e. what effect swapping indices in positions @c i  and @c j
  /// in <em>either bra or ket</em> has on the elements of the Tensor;
  /// Tensor's are <em>always assumed</em> to be particle-symmetric, i.e.
  /// swapping indices in positions @c i and @c j in <b>both bra and ket</b>;
  /// The allowed values are Symmetry::Symm, Symmetry::Antisymm, and
  /// Symmetry::Nonsymm
  /// @return the Symmetry object describing the symmetry of the bra and ket of
  /// the Tensor.
  Symmetry symmetry() const { return symmetry_; }
  /// @return the BraKetSymmetry object describing the symmetry of the Tensor
  /// under exchange of bra and ket.
  BraKetSymmetry braket_symmetry() const { return braket_symmetry_; }
  /// @return the ColumnSymmetry object describing the symmetry of the Tensor
  /// under exchange of _columns_ (i.e., pairs of matching {bra[i],ket[i]}
  /// slot bundles).
  ColumnSymmetry column_symmetry() const { return column_symmetry_; }

  /// @return number of bra slots (some may be occupied by null indices, hence
  /// this is the gross rank)
  std::size_t bra_rank() const { return bra_.size(); }
  /// @return number of nonnull bra indices (i.e. non-empty slots)
  std::size_t bra_net_rank() const { return bra_net_rank_; }
  /// @return number of ket slots (some may be occupied by null indices, hence
  /// this is the gross rank)
  std::size_t ket_rank() const { return ket_.size(); }
  /// @return number of nonnull ket indices  (i.e. non-empty slots)
  std::size_t ket_net_rank() const { return ket_net_rank_; }
  /// @return number of aux indices
  std::size_t aux_rank() const { return aux_.size(); }
  /// @return number of slots
  std::size_t num_slots() const {
    return bra_.size() + ket_.size() + aux_.size();
  }
  /// @return number of indices
  std::size_t num_indices() const {
    return std::ranges::count_if(
        slots(), [](const Index &idx) { return idx.nonnull(); });
  }
  /// @return number of indices in bra/ket
  /// @throw std::logic_error if bra and ket ranks do not match
  std::size_t rank() const {
    if (bra_rank() != ket_rank()) {
      throw std::logic_error("Tensor::rank(): bra rank != ket rank");
    }
    return bra_rank();
  }

  std::wstring to_latex() const override {
    auto &ctx = get_default_context();
    const auto bkt = ctx.braket_typesetting();
    const auto bkst = ctx.braket_slot_typesetting();

    // either rank > 1 or sum of bra and ket ranks > 1
    const bool add_bar =
        bra_rank() == ket_rank() ? rank() > 1 : bra_rank() + ket_rank() > 1;

    std::wstring core_label;
    if ((this->symmetry() == Symmetry::Antisymm) && add_bar)
      core_label += L"\\bar{";
    core_label += utf_to_latex(this->label());
    if ((this->symmetry() == Symmetry::Antisymm) && add_bar) core_label += L"}";

    switch (bkst) {
      case BraKetSlotTypesetting::Naive: {
        std::wstring result = L"{";
        result += core_label;

        // ket
        result += (bkt == BraKetTypesetting::KetSub ? L"_" : L"^");
        result += L"{";
        for (const auto &i : this->ket()) {
          result += i ? sequant::to_latex(i) : L"\\textvisiblespace";
        }
        result += L"}";

        // bra
        result += (bkt == BraKetTypesetting::BraSub ? L"_" : L"^");
        result += L"{";
        for (const auto &i : this->bra()) {
          result += i ? sequant::to_latex(i) : L"\\textvisiblespace";
        }
        result += L"}";

        // aux
        if (!this->aux_.empty()) {
          result += L"[";
          const index_container_type &__aux = this->aux();
          for (std::size_t i = 0; i < aux_rank(); ++i) {
            result += sequant::to_latex(__aux[i]);

            if (i + 1 < aux_rank()) {
              result += L",";
            }
          }
          result += L"]";
        }

        result += L"}";
        return result;
      }

      case BraKetSlotTypesetting::TensorPackage: {
        return to_latex_tensor(core_label, this->_bra(), this->_ket(),
                               this->_aux(), bkt, /* left_align = */ true);
      }
    }
    abort();  // unreachable
  }

  /// @note this performs rapid canonicalization only
  ExprPtr canonicalize(CanonicalizeOptions = {}) override;

  /// @brief adjoint of a Tensor swaps its bra and ket
  virtual void adjoint() override;

  /// Replaces indices using the index map
  /// @param index_map maps Index to Index
  /// @return true if one or more indices changed
  template <template <typename, typename, typename... Args> class Map,
            typename... Args>
  bool transform_indices(const Map<Index, Index, Args...> &index_map) {
    bool mutated = false;
    ranges::for_each(indices(), [&](auto &idx) {
      if (idx.transform(index_map)) mutated = true;
    });
    if (mutated) this->reset_hash_value();
    return mutated;
  }

  type_id_type type_id() const override { return get_type_id<Tensor>(); };

  ExprPtr clone() const override { return ex<Tensor>(*this); }

  void reset_tags() const {
    ranges::for_each(slots(), [](const auto &idx) { idx.reset_tag(); });
  }

  hash_type bra_hash_value() const {
    if (!bra_hash_value_)  // if hash not computed, or reset, recompute
      memoizing_hash();
    return *bra_hash_value_;
  }

  bool operator<(const AbstractTensor &other) const override final {
    auto *other_tensor = dynamic_cast<const Tensor *>(&other);
    if (other_tensor) {
      const Expr *other_expr = static_cast<const Expr *>(other_tensor);
      return this->static_less_than(*other_expr);
    } else
      return false;  // TODO do we compare typeid? labels? probably the latter
  }

 private:
  std::wstring label_{};
  sequant::bra<index_container_type> bra_{};
  sequant::ket<index_container_type> ket_{};
  sequant::aux<index_container_type> aux_{};
  Symmetry symmetry_ = Symmetry::Nonsymm;
  BraKetSymmetry braket_symmetry_ = BraKetSymmetry::Nonsymm;
  ColumnSymmetry column_symmetry_ = ColumnSymmetry::Nonsymm;
  mutable std::optional<hash_type>
      bra_hash_value_;  // memoized byproduct of memoizing_hash()
  bool is_adjoint_ = false;
  std::size_t bra_net_rank_;
  std::size_t ket_net_rank_;

  void validate_symmetries() {
    // (anti)symmetric bra or ket makes sense only for particle-symmetric
    // tensors
    if (symmetry_ == Symmetry::Symm || symmetry_ == Symmetry::Antisymm)
      assert(column_symmetry_ == ColumnSymmetry::Symm);
  }

  hash_type memoizing_hash() const override {
    auto compute_hash = [this]() {
      using std::begin;
      using std::end;
      auto val = hash::range(begin(bra()), end(bra()));
      bra_hash_value_ = val;
      hash::range(val, begin(ket()), end(ket()));
      hash::range(val, begin(aux()), end(aux()));
      hash::combine(val, label_);
      hash::combine(val, symmetry_);
      hash::combine(val, braket_symmetry_);
      hash::combine(val, column_symmetry_);
      // N.B. adjointness is baked into the label
      return val;
    };
    if (!hash_value_) {
      hash_value_ = compute_hash();
    } else {
      assert(*hash_value_ == compute_hash());
    }
    return *hash_value_;
  }
  void reset_hash_value() const override {
    Expr::reset_hash_value();
    bra_hash_value_.reset();
  }

  bool static_equal(const Expr &that) const override {
    const auto &that_cast = static_cast<const Tensor &>(that);
    if (this->label() == that_cast.label() &&
        this->symmetry() == that_cast.symmetry() &&
        this->braket_symmetry() == that_cast.braket_symmetry() &&
        this->column_symmetry() == that_cast.column_symmetry() &&
        this->bra_rank() == that_cast.bra_rank() &&
        this->ket_rank() == that_cast.ket_rank() &&
        this->aux_rank() == that_cast.aux_rank()) {
      // compare hash values first
      if (this->hash_value() ==
          that.hash_value())  // hash values agree -> do full comparison
        return this->bra() == that_cast.bra() &&
               this->ket() == that_cast.ket() && this->aux() == that_cast.aux();
      else
        return false;
    } else
      return false;
  }

  bool static_less_than(const Expr &that) const override {
    if (this == &that) return false;

    const auto &that_cast = static_cast<const Tensor &>(that);
    if (this->label() != that_cast.label()) {
      return this->label() < that_cast.label();
    }

    if (this->bra_rank() != that_cast.bra_rank()) {
      return this->bra_rank() < that_cast.bra_rank();
    }

    if (this->ket_rank() != that_cast.ket_rank()) {
      return this->ket_rank() < that_cast.ket_rank();
    }

    if (this->aux_rank() != that_cast.aux_rank()) {
      return this->aux_rank() < that_cast.aux_rank();
    }

    //          v1: compare hashes only
    //          return Expr::static_less_than(that);
    //          v2: compare fully
    if (this->bra_hash_value() != that_cast.bra_hash_value()) {
      return std::lexicographical_compare(
          this->bra().begin(), this->bra().end(), that_cast.bra().begin(),
          that_cast.bra().end());
    }

    if (this->ket() != that_cast.ket()) {
      return std::lexicographical_compare(
          this->ket().begin(), this->ket().end(), that_cast.ket().begin(),
          that_cast.ket().end());
    }

    return std::lexicographical_compare(this->aux().begin(), this->aux().end(),
                                        that_cast.aux().begin(),
                                        that_cast.aux().end());
  }

  Tensor *_clone() const override final { return new Tensor(*this); }
  std::shared_ptr<AbstractTensor> _clone_shared() const override final {
    return std::make_shared<Tensor>(*this);
  }

  // these implement the AbstractTensor interface
  AbstractTensor::const_any_view_randsz _bra() const override final {
    return ranges::counted_view<const Index *>(
        bra_.empty() ? nullptr : &(bra_[0]), bra_.size());
  }
  AbstractTensor::const_any_view_randsz _ket() const override final {
    return ranges::counted_view<const Index *>(
        ket_.empty() ? nullptr : &(ket_[0]), ket_.size());
  }
  AbstractTensor::const_any_view_randsz _aux() const override final {
    return ranges::counted_view<const Index *>(
        aux_.empty() ? nullptr : &(aux_[0]), aux_.size());
  }
  AbstractTensor::const_any_view_rand _braket() const override final {
    return braket();
  }
  AbstractTensor::const_any_view_rand _braketaux() const override final {
    return braketaux();
  }
  AbstractTensor::const_any_view_rand _slots() const override final {
    return slots();
  }
  std::size_t _bra_rank() const override final { return bra_rank(); }
  std::size_t _ket_rank() const override final { return ket_rank(); }
  std::size_t _bra_net_rank() const override final { return bra_net_rank(); }
  std::size_t _ket_net_rank() const override final { return ket_net_rank(); }
  std::size_t _aux_rank() const override final { return aux_rank(); }
  std::size_t _num_slots() const override final { return num_slots(); }
  std::size_t _num_indices() const override final { return num_indices(); }
  Symmetry _symmetry() const override final { return symmetry_; }
  BraKetSymmetry _braket_symmetry() const override final {
    return braket_symmetry_;
  }
  ColumnSymmetry _column_symmetry() const override final {
    return column_symmetry_;
  }
  std::size_t _color() const override final { return 0; }
  bool _is_cnumber() const override final { return true; }
  std::wstring_view _label() const override final { return label_; }
  std::wstring _to_latex() const override final { return to_latex(); }
  std::size_t _hash_value() const override final { return this->hash_value(); }
  bool _transform_indices(
      const container::map<Index, Index> &index_map) override final {
    return transform_indices(index_map);
  }
  void _reset_tags() override final { reset_tags(); }

  AbstractTensor::any_view_randsz _bra_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index *>(bra_.empty() ? nullptr : &(bra_[0]),
                                         bra_.size());
  }
  AbstractTensor::any_view_randsz _ket_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index *>(ket_.empty() ? nullptr : &(ket_[0]),
                                         ket_.size());
  }
  AbstractTensor::any_view_randsz _aux_mutable() override final {
    this->reset_hash_value();
    return ranges::counted_view<Index *>(aux_.empty() ? nullptr : &(aux_[0]),
                                         aux_.size());
  }

  /// swaps bra and ket slots
  void _swap_bra_ket() override final {
    this->reset_hash_value();
    std::swap(bra_.value(), ket_.value());
    std::swap(bra_net_rank_, ket_net_rank_);

    validate_indices();
    canonicalize_slots();
  }

};  // class Tensor

static_assert(is_tensor_v<Tensor>,
              "The Tensor class does not fulfill the requirements of the "
              "Tensor interface");

using TensorPtr = std::shared_ptr<Tensor>;

/// make_overlap tensor label is reserved since it is used by low-level SeQuant
/// machinery. Users can create make_overlap Tensor using make_overlap()
inline std::wstring overlap_label() { return L"s"; }

inline ExprPtr make_overlap(const Index &bra_index, const Index &ket_index) {
  return ex<Tensor>(Tensor(overlap_label(), bra{bra_index}, ket{ket_index},
                           aux{}, Tensor::reserved_tag{}));
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_TENSOR_HPP
