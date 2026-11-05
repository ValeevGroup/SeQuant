//
// Created by Eduard Valeyev on 8/25/25.
//

#ifndef SEQUANT_CORE_OPTIONS_HPP
#define SEQUANT_CORE_OPTIONS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer_fwd.hpp>
#include <SeQuant/core/utility/aggregate.hpp>

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace sequant {

/// canonicalization methods
enum class CanonicalizationMethod {
  /// Enables use of expression topology in the canonicalization, may be
  /// expensive.
  /// @note Canonicalization of tensor networks must take into account the
  ///       topology of the TN; this enables the use of canonical graph
  ///       sort of all TN elements. By itself it leaves the order of named
  ///       indices in (anti)symmetric slots to the graph, so its result
  ///       depends on how an expression is spelled (e.g. a lone tensor vs the
  ///       same tensor scaled) and may be aesthetically poor; Complete
  ///       canonicalization does not.
  Topological = 0b01,
  /// Enables the use of rapid canonicalization based on lexicographic sort.
  /// @note for TN this performs lexicographic sort of tensors, slots and
  /// slot bundles, and indices. Produces aesthetically pleasing result,
  /// but incomplete if some tensors are identical
  Lexicographic = 0b10,
  /// Enables use of Topological, then Lexicographic. Guaranteed to work even
  /// if some tensors are identical, and the result is aesthetically pleasing.
  Complete = Topological | Lexicographic,
  /// Lexicographic = quick-and-dirty
  Rapid = Lexicographic
};

CanonicalizationMethod operator&(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2);
CanonicalizationMethod operator|(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2);
std::wstring to_wstring(CanonicalizationMethod m);

/// tag that requests the named indices of an expression to be deduced as the
/// indices that occur once in each of its products
/// @warning deduction is a heuristic that fails for expressions in which an
/// external index occurs more than once (e.g. in tensor hypercontraction or
/// with quasi-diagonal tensors); the library itself never deduces
struct DeduceNamedIndices {
  explicit constexpr DeduceNamedIndices() = default;
  friend constexpr bool operator==(DeduceNamedIndices,
                                   DeduceNamedIndices) = default;
};

/// requests the deduction of named indices, see DeduceNamedIndices
inline constexpr DeduceNamedIndices deduce_named_indices{};

class NamedIndexSpec;

namespace detail {
/// @return the specification of the explicit named indices @p named_indices
/// @note for the library's own use: public functions take the named indices
/// as the indices of a ResultExpr, which also gives their symmetries
NamedIndexSpec explicit_named_indices(container::set<Index> named_indices);
}  // namespace detail

/// @brief specifies the named indices of a canonicalization, i.e. the indices
/// whose labels are meaningful and must not change (ordinarily the external
/// indices of an expression)
///
/// Either a request to deduce them (see DeduceNamedIndices), the only form a
/// user constructs, or an explicit set of indices, which the library derives
/// from a ResultExpr or from what it knows about an expression it builds
/// (see detail::explicit_named_indices()). There is deliberately no default:
/// every caller of a canonicalization states which indices are named.
class NamedIndexSpec {
 public:
  /// requests deduction
  NamedIndexSpec(DeduceNamedIndices) {}

  /// @return true if the named indices are to be deduced
  bool deduce() const { return !named_indices_.has_value(); }
  /// @return the named indices
  /// @pre `!deduce()`
  const container::set<Index>& indices() const;
  /// @return pointer to the named indices, or null if they are to be deduced
  const container::set<Index>* indices_ptr() const {
    return named_indices_ ? &*named_indices_ : nullptr;
  }

  friend bool operator==(const NamedIndexSpec&,
                         const NamedIndexSpec&) = default;

 private:
  explicit NamedIndexSpec(container::set<Index> named_indices)
      : named_indices_(std::move(named_indices)) {}
  friend NamedIndexSpec detail::explicit_named_indices(
      container::set<Index> named_indices);

  std::optional<container::set<Index>> named_indices_;
};

/// @brief settings of a canonicalization that only affect how its result
/// looks, i.e. which of the equivalent forms it picks; a member that is not
/// set takes the library default
struct CanonicalizationPresentation {
  SEQUANT_DESIGNATED_INIT_ONLY;

  /// the labels of the tensors to place first, in this order
  std::optional<container::vector<std::wstring>> cardinal_tensor_labels =
      std::nullopt;
  /// orders the indices of a tensor
  std::shared_ptr<const tensor_index_comparer_t> index_comparer = nullptr;
  /// orders the index pairs of a tensor
  std::shared_ptr<const tensor_index_pair_comparer_t> index_pair_comparer =
      nullptr;
  /// canonicalizers for the tensors with the given labels; the one for the
  /// empty label, if any, replaces the default canonicalizer of the tensors
  /// with any other label
  container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>
      tensor_canonicalizers = {};

  /// @return true if @p a and @p b agree in every member; comparers and
  /// canonicalizers are compared by identity
  friend bool operator==(const CanonicalizationPresentation& a,
                         const CanonicalizationPresentation& b);
};

/// @return the cardinal tensor labels of the library's canonical form: the
/// reserved labels of the antisymmetrizer, symmetrizer and transposition
/// operators
const container::vector<std::wstring>& default_cardinal_tensor_labels();

/// @brief options that control behavior of `canonicalize()`
struct CanonicalizeOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  enum class IgnoreNamedIndexLabel : bool { Yes = true, No = false };

  /// TN canonicalization method
  /// @internal
  /// TODO revert to CanonicalizationMethod::Topological once issues with
  /// external index handling exemplified by
  /// https://github.com/ValeevGroup/SeQuant/pull/406
  /// https://github.com/ValeevGroup/SeQuant/issues/426
  /// and https://github.com/ValeevGroup/SeQuant/issues/287
  /// is fixed
  /// @endinternal
  CanonicalizationMethod method = CanonicalizationMethod::Complete;
  /// whether the graph-based canonicalization ignores the labels of named
  /// indices, so that the structure of the result does not depend on them;
  /// named indices that an automorphism of the network can exchange are then
  /// placed in the order of their labels, as the context's index comparer
  /// orders them; the result is determined only if that comparer orders every
  /// two named indices of the network. A Sum canonicalizes its summands with
  /// labels distinguished regardless.
  IgnoreNamedIndexLabel ignore_named_index_labels = IgnoreNamedIndexLabel::Yes;
  /// presentation settings; null selects the library's canonical form, which
  /// is the only one recorded by the canonical mark of an Expr
  std::shared_ptr<const CanonicalizationPresentation> presentation = nullptr;

  /// @return true if these options produce the library's canonical form, i.e.
  /// have no presentation settings
  bool is_internal_form() const { return presentation == nullptr; }

  static CanonicalizeOptions default_options();
  CanonicalizeOptions copy_and_set(CanonicalizationMethod) const;
  CanonicalizeOptions copy_and_set(IgnoreNamedIndexLabel) const;
  CanonicalizeOptions copy_and_set(
      std::shared_ptr<const CanonicalizationPresentation>) const;

  /// @return true if @p a and @p b agree in every member; presentations are
  /// compared by value
  friend bool operator==(const CanonicalizeOptions& a,
                         const CanonicalizeOptions& b) {
    // the bindings name every member (guard is the SEQUANT_DESIGNATED_INIT_ONLY
    // tag), so that a member added to CanonicalizeOptions fails to compile here
    // until it is compared
    const auto& [guard_a, method_a, ignore_labels_a, presentation_a] = a;
    const auto& [guard_b, method_b, ignore_labels_b, presentation_b] = b;
    (void)guard_a;
    (void)guard_b;
    const bool same_presentation = presentation_a == presentation_b ||
                                   (presentation_a && presentation_b &&
                                    *presentation_a == *presentation_b);
    return method_a == method_b && ignore_labels_a == ignore_labels_b &&
           same_presentation;
  }
};

namespace detail {

/// the presentation settings that a canonicalization with given options uses;
/// each resolves a setting that @c opts.presentation does not set to the
/// library default
/// @{
container::vector<std::wstring> cardinal_tensor_labels(
    const CanonicalizeOptions& opts);
std::shared_ptr<const tensor_index_comparer_t> index_comparer(
    const CanonicalizeOptions& opts);
std::shared_ptr<const tensor_index_pair_comparer_t> index_pair_comparer(
    const CanonicalizeOptions& opts);
/// @return the canonicalizer for exactly @p label , or null
std::shared_ptr<TensorCanonicalizer> nondefault_tensor_canonicalizer(
    const CanonicalizeOptions& opts, std::wstring_view label);
/// @return the canonicalizer for @p label , else the one for the empty label,
/// or null if neither is set
std::shared_ptr<TensorCanonicalizer> tensor_canonicalizer(
    const CanonicalizeOptions& opts, std::wstring_view label);
/// @}

}  // namespace detail

}  // namespace sequant

#endif  // SEQUANT_CORE_OPTIONS_HPP
