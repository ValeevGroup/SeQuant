//
// Shared helpers for the basis-instance tests.
//

#ifndef SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP
#define SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/// test-local basis instances; SeQuant attaches no meaning to them
namespace test_csv {
inline constexpr sequant::IndexBasis::instance_type exchange = 1, coulomb = 2,
                                                    pno_pert = 10;
}  // namespace test_csv

namespace sequant::tests::csv {

/// calls @p f on every bra/ket/aux slot of every tensor of @p expr
template <typename F>
void for_each_index(ExprPtr const& expr, F&& f) {
  expr->visit(
      [&f](ExprPtr const& x) {
        if (!x->is<AbstractTensor>()) return;
        for (auto const& idx : x->as<AbstractTensor>()._slots()) f(idx);
      },
      /* atoms_only = */ true);
}

/// @return the number of tensor slots of @p expr per basis instance
/// (`std::nullopt` counts the slots without one)
inline std::map<IndexBasis::optional_instance, std::size_t> instance_histogram(
    ExprPtr const& expr) {
  std::map<IndexBasis::optional_instance, std::size_t> result;
  for_each_index(expr, [&result](Index const& idx) {
    ++result[idx.basis().basis_instance()];
  });
  return result;
}

/// the operator expression @p op with every pure-unoccupied (or, with
/// @p occupied, pure-occupied) index given @p inst -- what an OpRegistry
/// grant on that leg space makes OpMaker mint
inline ExprPtr with_leg_instance(ExprPtr const& op,
                                 IndexBasis::optional_instance inst,
                                 bool occupied = false) {
  auto const& isr = get_default_context().index_space_registry();
  container::map<Index, Index> m;
  for (auto const& idx : get_used_indices(op))
    if (occupied ? isr->is_pure_occupied(idx.space())
                 : isr->is_pure_unoccupied(idx.space()))
      m.emplace(idx, idx.replace_basis_instance(inst));
  return m.empty() ? op : transform_expr(op, m);
}

/// matches an exception whose message contains @p text
inline auto message_contains(std::string const& text) {
  return Catch::Matchers::MessageMatches(
      Catch::Matchers::ContainsSubstring(text));
}

/// the number of summands of @p expr if it is a Sum, else 1
inline std::size_t term_count(ExprPtr const& expr) {
  return expr->is<Sum>() ? expr->size() : 1;
}

/// the Context CSV-CCSD is derived under: min SR spaces plus the DF and PAO
/// spaces, Complete canonicalization, column-Symm deserialization
inline sequant::Context csv_cc_context() {
  auto isr = mbpt::make_min_sr_spaces();
  mbpt::add_df_spaces(isr);
  mbpt::add_pao_spaces(isr, IndexSpace::QuantumNumbers{mbpt::Spin::any});
  return Context({.index_space_registry_shared_ptr = std::move(isr),
                  .vacuum = Vacuum::SingleProduct,
                  .spbasis = SPBasis::Spinor,
                  .canonicalization_options =
                      CanonicalizeOptions::default_options().copy_and_set(
                          CanonicalizationMethod::Complete),
                  .deserialization_column_symmetry = ColumnSymmetry::Symm});
}

/// every tensor of @p expr labelled @p label
inline std::vector<AbstractTensor const*> tensors_labelled(
    ExprPtr const& expr, std::wstring_view label) {
  std::vector<AbstractTensor const*> result;
  expr->visit(
      [&](ExprPtr const& x) {
        if (!x->is<AbstractTensor>()) return;
        auto const& t = x->as<AbstractTensor>();
        if (std::wstring_view(t._label()) == label) result.push_back(&t);
      },
      /* atoms_only = */ true);
  return result;
}

/// the rank-1/1 overlaps of @p expr whose two ends carry different basis
/// instances
inline std::size_t standing_metrics(ExprPtr const& expr) {
  return std::ranges::count_if(
      tensors_labelled(expr, reserved::overlap_label()), [](auto const* t) {
        return t->_bra_rank() == 1 && t->_ket_rank() == 1 &&
               t->_bra()[0].basis().basis_instance() !=
                   t->_ket()[0].basis().basis_instance();
      });
}

/// the slots of the projector/symmetrizer tensors of @p term
inline std::vector<Index> projector_slots(ExprPtr const& term) {
  std::vector<Index> result;
  for (auto label : {reserved::antisymm_label(), reserved::symm_label()})
    for (AbstractTensor const* t : tensors_labelled(term, label))
      for (Index const& idx : t->_slots()) result.push_back(idx);
  return result;
}

}  // namespace sequant::tests::csv

#endif  // SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP
