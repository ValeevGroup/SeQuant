//
// Shared helpers for the basis-instance tests.
//

#ifndef SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP
#define SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
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

/// @p base with each of @p grants given on the particle space; a granted
/// perturbed amplitude (t¹, λ¹) first registers t¹, λ¹ and h¹ where absent
/// @pre the Context the registry will be used under is installed
inline std::shared_ptr<mbpt::OpRegistry> granted_registry(
    std::vector<std::pair<std::wstring, IndexBasis::instance_type>> const&
        grants,
    std::shared_ptr<mbpt::OpRegistry> base = mbpt::make_minimal_registry()) {
  const auto particles = get_particle_space(mbpt::Spin::any);
  for (auto const& [op, instance] : grants) {
    if ((op == L"t¹" || op == L"λ¹") && !base->contains(L"t¹"))
      base->add(L"t¹", mbpt::OpClass::Ex)
          .add(L"λ¹", mbpt::OpClass::Deex)
          .add(L"h¹", mbpt::OpClass::Gen);
    base->grant_basis(op, particles, instance);
  }
  return base;
}

/// installs csv_cc_context() and an mbpt Context with @p csv over
/// @p registry while this object lives
class ScopedCsvContext {
 public:
  explicit ScopedCsvContext(std::shared_ptr<mbpt::OpRegistry> registry =
                                mbpt::make_minimal_registry(),
                            mbpt::CSV csv = mbpt::CSV::Yes)
      : core_(sequant::set_scoped_default_context(csv_cc_context())),
        mbpt_(mbpt::set_scoped_default_mbpt_context(mbpt::Context(
            {.csv = csv, .op_registry_ptr = std::move(registry)}))) {}

 private:
  sequant::detail::ImplicitContextResetter<
      sequant::container::map<sequant::Statistics, sequant::Context>>
      core_;
  sequant::detail::ImplicitContextResetter<mbpt::Context> mbpt_;
};

/// @p equations, each closed-shell spin-traced
inline std::vector<ExprPtr> spintraced(std::vector<ExprPtr> const& equations) {
  std::vector<ExprPtr> result;
  for (auto const& eq : equations)
    result.push_back(mbpt::closed_shell_CC_spintrace_v2(eq));
  return result;
}

/// [E, R1, R2] of CC{2}.t(2, 0) under @p registry and the current CSV mode,
/// closed-shell spin-traced
inline std::vector<ExprPtr> derive_t(
    std::shared_ptr<mbpt::OpRegistry> registry) {
  Index::reset_tmp_index();
  auto mbpt_ctx = mbpt::set_scoped_default_mbpt_context(
      mbpt::Context({.csv = mbpt::get_default_mbpt_context().csv(),
                     .op_registry_ptr = std::move(registry)}));
  return spintraced(mbpt::CC{2}.t(2, 0));
}

/// CC{2}.λ() under @p registry and the current CSV mode, closed-shell
/// spin-traced
inline std::vector<ExprPtr> derive_λ(
    std::shared_ptr<mbpt::OpRegistry> registry) {
  Index::reset_tmp_index();
  auto mbpt_ctx = mbpt::set_scoped_default_mbpt_context(
      mbpt::Context({.csv = mbpt::get_default_mbpt_context().csv(),
                     .op_registry_ptr = std::move(registry)}));
  return spintraced(mbpt::CC{2}.λ());
}

}  // namespace sequant::tests::csv

#endif  // SEQUANT_TESTS_UNIT_CSV_TEST_UTILS_HPP
