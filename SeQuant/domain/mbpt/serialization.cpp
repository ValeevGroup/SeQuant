#include <SeQuant/domain/mbpt/serialization.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/function_registry.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <typeindex>
#include <vector>

namespace sequant::mbpt {

namespace {

using io::serialization::SerializationError;
using Ranks = std::vector<std::int64_t>;
using BatchOrdinals = container::svector<std::size_t>;

/// whether \p label (without adjoint mark and perturbation order) is made by
/// a factory of its own rather than by the generic one for its OpClass
bool is_special(std::wstring_view label) {
  return label == reserved::antisymm_label() ||
         label == reserved::symm_label() || label == L"ã";
}

/// Builds the operator labeled \p base (without adjoint mark and perturbation
/// order) with the given \p ranks as the standard factories do. \p ranks must
/// be valid for the operator, see candidate_ranks(), which also checks that the
/// registry defines it.
template <Statistics S>
ExprPtr make(const std::wstring &base, const Ranks &ranks,
             const OpParams &params) {
  const bool default_params =
      params.order == 0 && params.batch_ordinals.empty();
  const std::int64_t first = ranks.front();
  const std::int64_t second = ranks.size() > 1 ? ranks[1] : first;

  if (is_special(base)) {
    if constexpr (S == Statistics::FermiDirac) {
      if (!default_params)
        throw SerializationError(
            "operator '" + toUtf8(base) +
            "' has no perturbation order or batch ordinals");
      if (base == reserved::antisymm_label())
        return mbpt::A(nₚ(first), nₕ(second));
      if (base == reserved::symm_label()) return mbpt::S(first);
      return mbpt::ã(first);
    } else {
      throw SerializationError("operator '" + toUtf8(base) + "' is fermionic");
    }
  }

  // mirrors the generic factories, e.g. Hʼ(), tʼ(), r() and l()
  qns_t qns;
  std::size_t nc = 0, na = 0;
  switch (to_op_class(base)) {
    case OpClass::Gen:
      qns = general_type_qns(first);
      nc = na = first;
      break;
    case OpClass::Ex:
      qns = generic_excitation_qns(first, second, get_particle_space(Spin::any),
                                   get_hole_space(Spin::any));
      nc = first, na = second;
      break;
    case OpClass::Deex:
      qns =
          generic_deexcitation_qns(first, second, get_particle_space(Spin::any),
                                   get_hole_space(Spin::any));
      nc = second, na = first;
      break;
  }
  return ex<Operator<qns_t, S>>(
      [base]() -> std::wstring_view { return base; },
      [base, nc, na, params]() -> ExprPtr {
        return OpMaker<S>(base, ncre(nc), nann(na), params)();
      },
      [qns](qns_t &q) { q = combine(qns, q); }, params);
}

/// Builds the operator labeled \p label, which may carry an adjoint mark and a
/// perturbation order
template <Statistics S>
ExprPtr make_labeled(std::wstring_view label, const Ranks &ranks,
                     const BatchOrdinals &batch) {
  const bool adjoint = is_adjoint_label(label);
  const std::wstring_view unmarked = strip_adjoint_label(label);
  const std::wstring_view base = detail::strip_pert_order(unmarked);
  const std::size_t order =
      base.size() == unmarked.size()
          ? 0
          : detail::pert_superscripts.find(unmarked.back());

  ExprPtr op = make<S>(std::wstring(base), ranks,
                       {.order = order, .batch_ordinals = batch});
  if (adjoint) op->adjoint();
  return op;
}

/// @return the rank lists that are valid for the operator labeled \p base
/// (without adjoint mark and perturbation order) with ranks up to \p max
/// @throw SerializationError if the operator registry does not define \p base
std::vector<Ranks> candidate_ranks(std::wstring_view base, std::int64_t max) {
  if (!is_special(base) &&
      !get_default_mbpt_context().op_registry()->contains(std::wstring(base)))
    throw SerializationError("operator '" + toUtf8(base) +
                             "' is not in the operator registry");

  std::vector<Ranks> result;
  const bool pair =
      base == reserved::antisymm_label() ||
      (!is_special(base) && to_op_class(std::wstring(base)) != OpClass::Gen);
  const bool signed_ranks =
      base == reserved::antisymm_label() || base == reserved::symm_label();
  const std::int64_t min = signed_ranks ? -max : 0;
  for (std::int64_t first = min; first <= max; ++first) {
    if (!pair) {
      if (first != 0) result.push_back({first});
      continue;
    }
    for (std::int64_t second = min; second <= max; ++second) {
      if (first == 0 && second == 0) continue;
      if (signed_ranks && first * second < 0) continue;
      result.push_back(first == second ? Ranks{first} : Ranks{first, second});
    }
  }
  return result;
}

template <Statistics S>
ExprPtr read(const io::serialization::FunctionCall &call,
             const io::serialization::DeserializationOptions &) {
  auto args = io::serialization::split_arguments(call.arguments);
  if (args.empty())
    throw SerializationError("expected an operator label and its ranks");

  BatchOrdinals batch;
  if (args.back().starts_with("batch=")) {
    std::string_view list = args.back().substr(6);
    if (!list.starts_with('[') || !list.ends_with(']'))
      throw SerializationError("expected batch=[ordinals]");
    for (std::string_view ordinal :
         io::serialization::split_arguments(list.substr(1, list.size() - 2)))
      batch.push_back(string_to<std::size_t>(ordinal));
    args.pop_back();
  }

  const std::wstring label = io::serialization::parse_label(args.front());
  Ranks ranks;
  for (std::size_t i = 1; i < args.size(); ++i)
    ranks.push_back(string_to<std::int64_t>(args[i]));
  if (ranks.size() == 2 && ranks[0] == ranks[1]) ranks.pop_back();

  // only accept what candidate_ranks() considers valid, which is what the
  // factories accept
  const std::wstring_view base =
      detail::strip_pert_order(strip_adjoint_label(label));
  std::int64_t max = 0;
  for (std::int64_t rank : ranks) max = std::max(max, rank < 0 ? -rank : rank);
  const std::vector<Ranks> valid = candidate_ranks(base, max);
  if (std::ranges::find(valid, ranks) == valid.end())
    throw SerializationError("invalid ranks for operator '" + toUtf8(label) +
                             "'");

  return make_labeled<S>(label, ranks, batch);
}

template <Statistics S>
std::wstring write(const Expr &expr,
                   const io::serialization::SerializationOptions &) {
  const auto &op = expr.as<Operator<qns_t, S>>();
  const std::wstring_view label = op.label();
  const BatchOrdinals batch = op.batch_ordinals().value_or(BatchOrdinals{});

  // no rank of the operator exceeds the largest change of a quantum number
  std::int64_t max = 0;
  for (const auto &interval : op()) max = std::max(max, interval.upper());

  // find the factory call that reproduces the operator
  const std::wstring_view base =
      detail::strip_pert_order(strip_adjoint_label(label));
  for (const Ranks &ranks : candidate_ranks(base, max)) {
    if (!(*make_labeled<S>(label, ranks, batch) == expr)) continue;

    std::wstring result = io::serialization::quote_label(label);
    for (std::int64_t rank : ranks) result += L", " + std::to_wstring(rank);
    if (!batch.empty()) {
      result += L", batch=[";
      result += join_strings<std::wstring>(
          batch, L",", [](std::size_t ord) { return std::to_wstring(ord); });
      result += L"]";
    }
    return result;
  }
  throw SerializationError("mbpt operator '" + toUtf8(label) +
                           "' cannot be reproduced by the standard factories");
}

}  // namespace

void register_serialization_functions() {
  using enum Statistics;
  io::serialization::register_function(
      L"fmbop", std::type_index(typeid(Operator<qns_t, FermiDirac>)),
      read<FermiDirac>, write<FermiDirac>);
  io::serialization::register_function(
      L"bmbop", std::type_index(typeid(Operator<qns_t, BoseEinstein>)),
      read<BoseEinstein>, write<BoseEinstein>);
}

}  // namespace sequant::mbpt
