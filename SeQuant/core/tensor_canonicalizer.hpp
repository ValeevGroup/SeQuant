//
// Created by Robert Adam on 2023-09-08
//

#ifndef SEQUANT_CORE_TENSOR_CANONICALIZER_HPP
#define SEQUANT_CORE_TENSOR_CANONICALIZER_HPP

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/counted.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/zip.hpp>

#include <cstdint>
#include <functional>
#include <memory>
#include <string_view>

namespace sequant {

class AbstractTensor;

/// @return true for reserved bookkeeping operators ((anti)symmetrizer,
///         transposition) whose bra<->ket orientation defines/extracts
///         external indices: canonicalization must never reorient them.
///         Their Conjugate braket symmetry is the reserved Symm->Conjugate
///         demotion sentinel (see Tensor's constructor), not a foldable
///         value symmetry.
bool braket_orientation_pinned(const AbstractTensor& t);

/// @brief Base class for Tensor canonicalizers
/// To make custom canonicalizer make a derived class and register an instance
/// of that class with TensorCanonicalizer::register_instance
class TensorCanonicalizer {
 public:
  using index_comparer_t = std::function<bool(const Index&, const Index&)>;
  using index_pair_t = std::pair<const Index, const Index>;
  using index_pair_comparer_t =
      std::function<bool(const index_pair_t&, const index_pair_t)>;

  virtual ~TensorCanonicalizer();

  /// @return ptr to the TensorCanonicalizer object, if any, that had been
  /// previously registered via TensorCanonicalizer::register_instance()
  /// with @c label , or to the default canonicalizer, if any
  static std::shared_ptr<TensorCanonicalizer> instance_ptr(
      std::wstring_view label = L"");

  /// @return ptr to the TensorCanonicalizer object, if any, that had been
  /// previously registered via TensorCanonicalizer::register_instance()
  /// with @c label
  /// @sa instance_ptr
  static std::shared_ptr<TensorCanonicalizer> nondefault_instance_ptr(
      std::wstring_view label);

  /// @return a TensorCanonicalizer previously registered via
  /// TensorCanonicalizer::register_instance() with @c label or to the default
  /// canonicalizer
  /// @throw Exception if no canonicalizer has been registered
  static std::shared_ptr<TensorCanonicalizer> instance(
      std::wstring_view label = L"");

  /// registers @c canonicalizer to be applied to Tensor objects with label
  /// @c label ; leave the label empty if @c canonicalizer is to apply to Tensor
  /// objects with any label
  /// @note if a canonicalizer registered with label @c label exists, it is
  /// replaced
  static void register_instance(
      std::shared_ptr<TensorCanonicalizer> canonicalizer,
      std::wstring_view label = L"");

  /// tries to register @c canonicalizer to be applied to Tensor objects
  /// with label @c label ; leave the label empty if @c canonicalizer is to
  /// apply to Tensor objects with any label
  /// @return false if there is already a canonicalizer registered with @c label
  /// @sa regiter_instance
  static bool try_register_instance(
      std::shared_ptr<TensorCanonicalizer> canonicalizer,
      std::wstring_view label = L"");

  /// deregisters canonicalizer (if any) registered previously
  /// to be applied to tensors with label @c label
  static void deregister_instance(std::wstring_view label = L"");

  /// @return a list of Tensor labels with lexicographic preference (in order)
  static const auto& cardinal_tensor_labels() {
    return cardinal_tensor_labels_accessor();
  }

  /// @brief Sets cardinal tensor labels by appending to defaults
  /// @param labels a list of additional Tensor labels with lexicographic
  /// preference (in order)
  /// @note The default labels are always prepended to the provided labels
  /// @note To restore defaults only, use reset_cardinal_tensor_labels()
  static void set_cardinal_tensor_labels(
      const container::vector<std::wstring>& labels);

  /// @brief Resets cardinal tensor labels to default values only
  static void reset_cardinal_tensor_labels();

  /// @brief Clears all cardinal tensor labels including defaults (sets to empty
  /// list)
  static void clear_all_cardinal_tensor_labels();

  /// @return a side effect of canonicalization (e.g. phase), or nullptr if none
  /// @internal what should be returned if canonicalization requires
  /// complex conjugation? Special ExprPtr type (e.g. ConjOp)? Or the actual
  /// return of the canonicalization?
  /// @note canonicalization compared indices returned by index_comparer
  // TODO generalize for complex tensors
  virtual ExprPtr apply(AbstractTensor&) const = 0;

  /// @return reference to the object used to compare Index objects
  static const index_comparer_t& index_comparer();

  /// @param comparer the compare object to be used by this
  static void index_comparer(index_comparer_t comparer);

  /// @return reference to the object used to compare Index objects
  static const index_pair_comparer_t& index_pair_comparer();

  /// @param comparer the compare object to be used by this
  static void index_pair_comparer(index_pair_comparer_t comparer);

 protected:
  static inline auto mutable_bra_range(AbstractTensor& t) {
    return t._bra_mutable();
  }
  static inline auto mutable_ket_range(AbstractTensor& t) {
    return t._ket_mutable();
  }
  static inline auto mutable_aux_range(AbstractTensor& t) {
    return t._aux_mutable();
  }

  /// the object used to compare indices
  static index_comparer_t index_comparer_;
  /// the object used to compare pairs of indices
  static index_pair_comparer_t index_pair_comparer_;

 private:
  static std::pair<
      container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>*,
      std::unique_lock<std::recursive_mutex>>
  instance_map_accessor();  // map* + locked recursive mutex
  static container::vector<std::wstring>& cardinal_tensor_labels_accessor();
  static container::vector<std::wstring>&
  default_cardinal_tensor_labels_accessor();
};

/// @brief null Tensor canonicalizer does nothing
class NullTensorCanonicalizer : public TensorCanonicalizer {
 public:
  virtual ~NullTensorCanonicalizer() = default;

  ExprPtr apply(AbstractTensor&) const override;
};

/// @return whether @p t 's braket orientation is pinned: the reserved
///         (anti)symmetrization/transposition bookkeeping operators, whose
///         bra<->ket orientation defines/extracts external indices and must
///         never be reoriented
bool braket_orientation_pinned(const AbstractTensor& t);

/// @return whether the bra<->ket exchange is a respelling of @p t, i.e. a
///         braket symmetry with a braket_swap_sign() (BraKetSymmetry::Symm,
///         free, or BraKetSymmetry::Antisymm, carrying -1) on a tensor whose
///         orientation is not pinned. A Conjugate/AntiConjugate tensor's two
///         orientations are two values (`T{q;p} = s conj(T{p;q})`), so they
///         are never exchanged. Every fold site gates on this predicate and
///         records the sign DefaultTensorCanonicalizer::canonicalize_braket()
///         returns.
bool braket_foldable(const AbstractTensor& t);

/// @brief the elementwise conjugation of @p t's array, spelled on the core
/// conjugation model: the adjoint state with the bra/ket bundles exchanged,
/// i.e. `conj T{a;b}` is spelled `T⁺{b;a}` (the definition of the adjoint
/// state; the network contracts by labels, so the two are one labeled
/// array). A definite hermiticity consumes the state in the exchange itself
/// (`T{b;a} = s conj T{a;b}`), leaving the bundles exchanged; over a real
/// basis the coset rule turns it into `꙳` with the slots in place. An
/// involution. Used by the Kramers network/block fold and the Kramers
/// tracer, whose time-reversal image is the conjugate of the flavor-flipped
/// array. The Kramers orientation verdicts (kramers_flavor_key,
/// kramers_noncanonical) read a marked tensor in its VALUE orientation, so a
/// folded spelling is a fixed point of the fold.
/// @return the sign consumed by a definite hermiticity (-1 for an
///         anti-Hermitian @p t), +1 otherwise
std::int8_t kramers_conjugate_mark(AbstractTensor& t);

/// @return whether the Kramers fold may spell @p t's bra<->ket exchange: a
///         c-number, not orientation-pinned, whose exchange is a bare swap
///         (braket_foldable) or the adjoint itself (a definite hermiticity,
///         BraKetSymmetry::Conjugate/AntiConjugate, the fold's own mark)
bool kramers_orientation_free(const AbstractTensor& t);

/// @brief eval-boundary respelling of a definite-hermiticity Kramers tensor
/// whose bra carries more down-flavored slots than its ket: the bundles are
/// exchanged through adjoint(), which such a hermiticity consumes
/// (`T{b;a} = s conj T{a;b}`), so the tensor left behind is the up-row
/// spelling a leaf provider serves and the as-written value is that spelling
/// conjugated, transposed and scaled by the returned sign -- the caller
/// records `{conj, braket_swap, s}` as the retrieval transform. Not a
/// symbolic respelling: the two orientations of such a tensor are two values.
/// No-op for any other tensor (kramers_foldable() false, no definite
/// hermiticity, or a bra with no more down slots than the ket).
/// @return {whether the exchange was made, the sign it consumed}
std::pair<bool, std::int8_t> kramers_uprow_exchange(AbstractTensor& t);

/// @return whether the Kramers (time-reversal) fold applies to @p t:
///         KramersSymmetry::TimeReversal, a c-number, and not
///         orientation-pinned (reserved operators never fold)
bool kramers_foldable(const AbstractTensor& t);

/// @brief Kramers (time-reversal) fold of a single tensor: if @p t's FIRST
/// flavored slot (bra, ket, aux order) carries the non-canonical (down)
/// flavor, every flavored slot index is replaced by its Kramers partner
/// (proto indices flipped recursively, unflavored slots untouched) and, with
/// @p mark, the elementwise conjugation is spelled on the tensor
/// (kramers_conjugate_mark), preserving the value up to the returned phase:
/// T = phase * conj(T_flipped), phase = (-1)^(#slots flipped from down) times
/// the sign the mark consumes. With @p mark false only the slots are flipped
/// and the caller owns the conjugation (the eval boundary records it as the
/// leaf's retrieval transform). No-op (phase +1) if the fold does not apply
/// or the first flavored slot is already canonical. Idempotent.
/// @return the phase (+1 or -1)
int canonicalize_kramers(AbstractTensor& t, bool mark = true);

/// @brief flips every flavored slot index of @p t to its Kramers partner in
/// place (proto indices recursively, unflavored slots untouched); no marker
/// or phase bookkeeping -- an involution used by canonicalize_kramers() and
/// by consumers that must restore a folded spelling
/// @return whether any slot was flipped
bool kramers_flip_slots(AbstractTensor& t);

/// @brief the DEEP variant of kramers_flip_slots: every slot index is replaced
/// by its kramers_flipped_deep image, so the flavoured proto indices of an
/// unflavoured (e.g. Kramers-union) composite are flipped too -- the
/// whole-expression time-reversal flip the KramersFlip fold builds its
/// canonical partner with (eval_expr.cpp); no marker or phase bookkeeping
/// @return whether any slot changed
bool kramers_flip_slots_deep(AbstractTensor& t);

/// @brief whether @p idx is a Kramers-UNION index: a spin-free index (its
/// space has no Kramers partner) of a space whose flavoured subspaces ARE
/// registered Kramers partners in @p isr, i.e. the union of the ↑ and ↓
/// halves of one space (the expansion dummy of a Kramers-union CSV
/// transform, a union-contracted integral leg). The time-reversal image of a
/// tensor over such an axis permutes the axis (swaps the two halves, with a
/// sign), so it is NOT an elementwise {conj, phase} of the tensor: a leaf
/// with a union slot must not be Kramers-folded at the eval-leaf level (the
/// network-level fold, which relabels the summed union dummy consistently
/// across the term, is unaffected). Indices of spaces without flavoured
/// partners (e.g. a density-fitting auxiliary index) are not union indices.
bool kramers_union_index(const Index& idx, const IndexSpaceRegistry& isr);

/// @return true if any slot (bra, ket, aux) of @p t is a Kramers-union index
/// (see kramers_union_index) under the default context's registry
bool has_kramers_union_slot(const AbstractTensor& t);

/// @brief flavor key of @p t: label + per-bundle flavor characters
/// ('a'/'b'/'-' for up/down/unflavored, so up orders first) SORTED within
/// each bundle, the bra
/// and ket bundles ordered canonically for braket-foldable tensors -- hence
/// invariant under every symmetry the canonicalizer may exercise
/// (within-bundle permutation, bra<->ket exchange) and under index
/// relabeling
/// @param flipped if true, the key of the Kramers-flipped spelling
std::wstring kramers_flavor_key(const AbstractTensor& t, bool flipped = false);

/// @return whether @p t is spelled in its non-canonical Kramers orientation:
///         more down- than up-flavored slots, or (tie) the flipped flavor
///         key orders before its own (see kramers_flavor_key); false for
///         tensors without flavored slots
bool kramers_noncanonical(const AbstractTensor& t);

class DefaultTensorCanonicalizer : public TensorCanonicalizer {
 public:
  DefaultTensorCanonicalizer() = default;

  /// @tparam IndexContainer a Container of Index objects such that @c
  /// IndexContainer::value_type is convertible to Index (e.g. this can be
  /// std::vector or std::set , but not std::map)
  /// @param external_indices container of external Index objects
  /// @warning @c external_indices is assumed to be immutable during the
  /// lifetime of this object
  template <typename IndexContainer>
  DefaultTensorCanonicalizer(IndexContainer&& external_indices) {
    ranges::for_each(external_indices, [this](const Index& idx) {
      this->external_indices_.emplace(idx);
    });
  }
  virtual ~DefaultTensorCanonicalizer() = default;

  /// Canonicalizes the assignment of indices to bra and ket of a
  /// braket-foldable tensor (see braket_foldable()): a bare swap for
  /// BraKetSymmetry::Symm/Antisymm. Every other tensor, a
  /// BraKetSymmetry::Conjugate/AntiConjugate one among them, is left as
  /// written
  /// @param fold_signed if false, a respelling that costs a sign (an
  ///        Antisymm tensor) is left untouched
  /// @return the sign the respelling contributed (+1 or -1): the tensor as it
  ///         stood is that sign times the tensor this leaves behind. The
  ///         caller must record it, e.g. in the phase byproduct of apply()
  static std::int8_t canonicalize_braket(AbstractTensor& t,
                                         bool fold_signed = true);

  /// Implements TensorCanonicalizer::apply
  /// @note Canonicalizes @c t by sorting its bra (if @c
  /// t.symmetry()==Symmetry::Nonsymm ) or its bra and ket (if @c
  /// t.symmetry()!=Symmetry::Nonsymm ),
  ///       with the external indices appearing "before" (smaller particle
  ///       indices) than the internal indices
  ExprPtr apply(AbstractTensor& t) const override;

  /// Core of DefaultTensorCanonicalizer::apply, only does the canonicalization,
  /// i.e. no tagging/untagging
  template <typename IndexComp, typename IndexPairComp>
  ExprPtr apply(AbstractTensor& t, const IndexComp& idxcmp,
                const IndexPairComp& paircmp) const {
    // nothing to do for non-particle-symmetric tensors
    if (t._column_symmetry() == ColumnSymmetry::Nonsymm) return nullptr;

    auto s = symmetry(t);
    auto is_antisymm = (s == Symmetry::Antisymm);
    const auto _bra_rank = bra_rank(t);
    const auto _ket_rank = ket_rank(t);
    [[maybe_unused]] const auto _aux_rank = aux_rank(t);
    const auto _rank = std::min(_bra_rank, _ket_rank);

    // nothing to do for rank-1 tensors
    if (_bra_rank == 1 && _ket_rank == 1) return nullptr;

    using ranges::begin;
    using ranges::end;
    using ranges::views::counted;
    using ranges::views::take;
    using ranges::views::zip;

    bool even = true;
    switch (s) {
      case Symmetry::Antisymm:
      case Symmetry::Symm: {
        auto _bra = mutable_bra_range(t);
        auto _ket = mutable_ket_range(t);
        // std::{stable_}sort does not necessarily use swap! so must implement
        // sort ourselves .. thankfully ranks will be low so can stick with
        // bubble
        const int parity =
            bubble_sort_parity(_bra, idxcmp) * bubble_sort_parity(_ket, idxcmp);
        if (is_antisymm) even = parity == 1;
      } break;

      case Symmetry::Nonsymm: {
        // sort particles with bra and ket functions first,
        // then the particles with either bra or ket index
        auto _bra = mutable_bra_range(t);
        auto _ket = mutable_ket_range(t);
        auto _zip_braket = zip(take(_bra, _rank), take(_ket, _rank));
        bubble_sort(begin(_zip_braket), end(_zip_braket), paircmp);
        if (_bra_rank > _rank) {
          auto size_of_rest = _bra_rank - _rank;
          auto rest_of = counted(begin(_bra) + _rank, size_of_rest);
          bubble_sort(begin(rest_of), end(rest_of), idxcmp);
        } else if (_ket_rank > _rank) {
          auto size_of_rest = _ket_rank - _rank;
          auto rest_of = counted(begin(_ket) + _rank, size_of_rest);
          bubble_sort(begin(rest_of), end(rest_of), idxcmp);
        }
      } break;
    }

    // TODO: Handle auxiliary index symmetries once they are introduced
    // auto _aux = mutable_aux_range(t);
    // ranges::sort(_aux, comp);

    ExprPtr result =
        is_antisymm ? (even == false ? ex<Constant>(-1) : nullptr) : nullptr;
    return result;
  }

 private:
  container::set<Index> external_indices_;

 protected:
  void tag_indices(AbstractTensor& t) const;
};

class TensorBlockCanonicalizer : public DefaultTensorCanonicalizer {
 public:
  TensorBlockCanonicalizer() = default;
  ~TensorBlockCanonicalizer() = default;

  /// \param fold_signed_braket if false, canonicalize_braket leaves a tensor
  ///        whose bra<->ket exchange costs a sign (Antisymm) untouched, so
  ///        that only sign-free respellings (Symm) fold. Eval-boundary
  ///        bridge: the two orientations of such a tensor stay separate
  ///        spellings, each asked of a leaf provider as written, while the
  ///        phase of the respellings that do happen composes into the leaf's
  ///        retrieval transform and so reaches its value (normalize_leaf()
  ///        in eval_expr.cpp).
  explicit TensorBlockCanonicalizer(bool fold_signed_braket)
      : fold_signed_braket_(fold_signed_braket) {}

  template <typename IndexContainer>
  TensorBlockCanonicalizer(const IndexContainer& external_indices)
      : DefaultTensorCanonicalizer(external_indices) {}

  ExprPtr apply(AbstractTensor& t) const override;

  /// @param fold_kramers if true, canonicalize_kramers() is applied; OFF by
  ///        default: inside a network the per-tensor fold would flip one
  ///        tensor's dummies but not its partner's (the network fold owns
  ///        that decision); the eval leaf boundary opts in explicitly
  TensorBlockCanonicalizer& fold_kramers(bool fold_kramers) {
    fold_kramers_ = fold_kramers;
    return *this;
  }

 private:
  bool fold_signed_braket_ = true;
  bool fold_kramers_ = false;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_TENSOR_CANONICALIZER_HPP
