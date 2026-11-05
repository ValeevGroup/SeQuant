//
// Created by Eduard Valeyev on 8/25/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <range/v3/algorithm/equal.hpp>

namespace sequant {

CanonicalizationMethod operator&(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2) {
  return static_cast<CanonicalizationMethod>(static_cast<int>(m1) &
                                             static_cast<int>(m2));
}

CanonicalizationMethod operator|(CanonicalizationMethod m1,
                                 CanonicalizationMethod m2) {
  return static_cast<CanonicalizationMethod>(static_cast<int>(m1) |
                                             static_cast<int>(m2));
}

std::wstring to_wstring(CanonicalizationMethod m) {
  switch (m) {
    case CanonicalizationMethod::Topological:
      return L"topological";
    case CanonicalizationMethod::Lexicographic:
      return L"lexicographic";
    case CanonicalizationMethod::Complete:
      return L"complete";
  }

  SEQUANT_UNREACHABLE;
}

bool operator==(const CanonicalizationPresentation& a,
                const CanonicalizationPresentation& b) {
  // the bindings name every member (guard is the SEQUANT_DESIGNATED_INIT_ONLY
  // tag), so that a member added to CanonicalizationPresentation fails to
  // compile here until it is compared
  const auto& [guard_a, labels_a, cmp_a, pair_cmp_a, canonicalizers_a] = a;
  const auto& [guard_b, labels_b, cmp_b, pair_cmp_b, canonicalizers_b] = b;
  (void)guard_a;
  (void)guard_b;
  return labels_a == labels_b && cmp_a.get() == cmp_b.get() &&
         pair_cmp_a.get() == pair_cmp_b.get() &&
         ranges::equal(canonicalizers_a, canonicalizers_b,
                       [](const auto& x, const auto& y) {
                         return x.first == y.first &&
                                x.second.get() == y.second.get();
                       });
}

const container::set<Index>& NamedIndexSpec::indices() const {
  if (!named_indices_)
    throw Exception(
        "NamedIndexSpec::indices: the named indices are to be deduced");
  return *named_indices_;
}

const container::vector<std::wstring>& default_cardinal_tensor_labels() {
  static const container::vector<std::wstring> result{
      reserved::antisymm_label(), reserved::symm_label(),
      reserved::transposition_label()};
  return result;
}

CanonicalizeOptions CanonicalizeOptions::default_options() {
  return sequant::get_default_context().canonicalization_options().value_or(
      CanonicalizeOptions{});
}

CanonicalizeOptions CanonicalizeOptions::copy_and_set(
    CanonicalizationMethod arg) const {
  auto result = *this;
  result.method = arg;
  return result;
}

CanonicalizeOptions CanonicalizeOptions::copy_and_set(
    IgnoreNamedIndexLabel arg) const {
  auto result = *this;
  result.ignore_named_index_labels = arg;
  return result;
}

CanonicalizeOptions CanonicalizeOptions::copy_and_set(
    std::shared_ptr<const CanonicalizationPresentation> arg) const {
  auto result = *this;
  result.presentation = std::move(arg);
  return result;
}

namespace detail {

// TODO: options without presentation settings still take them from the
// default context, until the context no longer holds them

container::vector<std::wstring> cardinal_tensor_labels(
    const CanonicalizeOptions& opts) {
  if (!opts.presentation)
    return get_default_context_snapshot().cardinal_tensor_labels();
  return opts.presentation->cardinal_tensor_labels.value_or(
      default_cardinal_tensor_labels());
}

std::shared_ptr<const tensor_index_comparer_t> index_comparer(
    const CanonicalizeOptions& opts) {
  if (!opts.presentation)
    return get_default_context_snapshot().index_comparer_ptr();
  if (opts.presentation->index_comparer)
    return opts.presentation->index_comparer;
  static const auto default_comparer =
      std::make_shared<const tensor_index_comparer_t>(
          TensorCanonicalizer::default_index_comparer());
  return default_comparer;
}

std::shared_ptr<const tensor_index_pair_comparer_t> index_pair_comparer(
    const CanonicalizeOptions& opts) {
  if (!opts.presentation)
    return get_default_context_snapshot().index_pair_comparer_ptr();
  if (opts.presentation->index_pair_comparer)
    return opts.presentation->index_pair_comparer;
  static const auto default_comparer =
      std::make_shared<const tensor_index_pair_comparer_t>(
          TensorCanonicalizer::default_index_pair_comparer());
  return default_comparer;
}

std::shared_ptr<TensorCanonicalizer> nondefault_tensor_canonicalizer(
    const CanonicalizeOptions& opts, std::wstring_view label) {
  if (!opts.presentation)
    return get_default_context_snapshot().nondefault_tensor_canonicalizer_ptr(
        label);
  const auto& map = opts.presentation->tensor_canonicalizers;
  auto it = map.find(std::wstring{label});
  return it != map.end() ? it->second : nullptr;
}

std::shared_ptr<TensorCanonicalizer> tensor_canonicalizer(
    const CanonicalizeOptions& opts, std::wstring_view label) {
  auto result = nondefault_tensor_canonicalizer(opts, label);
  if (!result) result = nondefault_tensor_canonicalizer(opts, L"");
  return result;
}

}  // namespace detail

}  // namespace sequant
