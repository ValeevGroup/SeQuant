#ifndef SEQUANT_TESTS_UNIT_DERIV_CCSD_HELPERS_HPP
#define SEQUANT_TESTS_UNIT_DERIV_CCSD_HELPERS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace sequant::deriv::test {

/// @return every tensor leaf of @p equations but a projector or a Kronecker
/// delta, one per spelling, with indices of its own: its slot indices and
/// their protoindices renamed, per space and in the order they are written,
/// to ordinals from 51 on
inline container::svector<ExprPtr> leaves_as_variables(
    const std::vector<std::pair<std::wstring, ExprPtr>>& equations) {
  container::map<std::wstring, ExprPtr> variables;
  for (const auto& [name, e] : equations)
    e->visit(
        [&variables](const ExprPtr& leaf) {
          if (!leaf->is<Tensor>()) return;
          const auto& t = leaf->as<Tensor>();
          if (t.label() == reserved::antisymm_label() ||
              t.label() == reserved::symm_label() ||
              t.label() == reserved::kronecker_label())
            return;
          container::map<Index, Index> names;
          container::map<IndexSpace, std::size_t> count;
          auto name = [&](const Index& i) {
            if (!names.contains(i))
              names.emplace(i, Index(i.basis(), 50 + ++count[i.space()]));
            return names.at(i);
          };
          for (const auto& i : t.const_braketaux()) {
            if (!i.has_proto_indices() || names.contains(i)) {
              name(i);
              continue;
            }
            container::vector<Index> protos;
            for (const auto& p : i.proto_indices()) protos.push_back(name(p));
            names.emplace(i, Index(name(i.drop_proto_indices()), protos,
                                   i.symmetric_proto_indices()));
          }
          const auto x = transform_expr(leaf, names);
          variables.emplace(serialize(x), x);
        },
        /*atoms_only=*/true);
  container::svector<ExprPtr> result;
  for (const auto& [spelling, x] : variables) result.push_back(x);
  return result;
}

/// the equations SeQuant generates for a CC model, by name, and variables to
/// differentiate them by
struct GeneratedEquations {
  std::vector<std::pair<std::wstring, ExprPtr>> equations;
  container::svector<ExprPtr> variables;
};

/// @return under the current context, the equations of CC{2}: ℒ as ref_av
/// generates it, ℒ simplified, and every t()[k], λ()[k], tʼ(1, 1)[k],
/// λʼ(1, 1)[k]; closed-shell spin-traced if @p traced. The variables are t_1,
/// t_2, λ_1, λ_2 (spelled over a basis of their occupied indices if @p csv)
/// and, unless @p csv, the integrals f and g.
inline GeneratedEquations generated_equations(bool csv, bool traced) {
  mbpt::CC cc{2};
  GeneratedEquations result;
  auto& equations = result.equations;
  const auto raw = mbpt::op::ref_av((1 + mbpt::op::Λ(2)) * cc.hbar(),
                                    {.connect = cc.hbar_connections()});
  equations.emplace_back(L"ℒ raw", raw);
  equations.emplace_back(L"ℒ", simplify(raw->clone()));
  const auto add = [&equations](const wchar_t* name,
                                const std::vector<ExprPtr>& es) {
    for (std::size_t k = 0; k < es.size(); ++k)
      if (es[k])
        equations.emplace_back(
            std::wstring(name) + L"[" + std::to_wstring(k) + L"]", es[k]);
  };
  add(L"t", cc.t());
  add(L"λ", cc.λ());
  add(L"tʼ", cc.tʼ(1, 1));
  add(L"λʼ", cc.λʼ(1, 1));
  if (traced)
    for (auto& [name, e] : equations) {
      const bool scalar = name.starts_with(L"ℒ") || name.ends_with(L"[0]");
      e = scalar ? mbpt::closed_shell_spintrace(e, {})
                 : mbpt::closed_shell_CC_spintrace_v2(e);
      if (name != L"ℒ raw") e = simplify(e);
    }

  const std::wstring amplitude_symmetry = traced ? L":N-N-S" : L":A-N-S";
  const std::wstring integral_symmetry = traced ? L":N-C-S" : L":A-C-S";
  const auto amplitudes =
      csv ? std::vector<std::wstring>{L"t{a_51<i_51>;i_51}",
                                      L"t{a_51<i_51,i_52>,a_52<i_51,i_52>;"
                                      L"i_51,i_52}",
                                      L"λ{i_51;a_51<i_51>}",
                                      L"λ{i_51,i_52;a_51<i_51,i_52>,"
                                      L"a_52<i_51,i_52>}"}
          : std::vector<std::wstring>{
                L"t{a_51;i_51}", L"t{a_51,a_52;i_51,i_52}", L"λ{i_51;a_51}",
                L"λ{i_51,i_52;a_51,a_52}"};
  for (const auto& x : amplitudes)
    result.variables.push_back(deserialize(x + amplitude_symmetry));
  if (!csv)
    for (const auto* x : {L"f{i_51;a_51}", L"g{i_51,i_52;a_51,a_52}"})
      result.variables.push_back(deserialize(x + integral_symmetry));
  return result;
}

/// @return the canonical open-shell CC{2} equations: each spin case of ℒ and of
/// every t()[k], λ()[k], tʼ(1, 1)[k], λʼ(1, 1)[k] (k = 1, 2) traced by
/// open_shell_CC_spintrace and simplified, made under the Complete
/// canonicalization SeQuant's open-shell code requires. The variables are the
/// tensor leaves of the equations, one per spelling, with indices of their own
inline GeneratedEquations open_shell_equations() {
  auto context = get_default_context();
  context.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  const auto scope = set_scoped_default_context(context);
  Index::reset_tmp_index();
  mbpt::CC cc{2};
  GeneratedEquations result;
  const auto add = [&result](const std::wstring& name, const ExprPtr& e) {
    const auto cases = mbpt::open_shell_CC_spintrace(e);
    for (std::size_t c = 0; c != cases.size(); ++c)
      result.equations.emplace_back(name + L" case " + std::to_wstring(c),
                                    simplify(cases[c]->clone()));
  };
  add(L"ℒ", simplify(mbpt::op::ref_av((1 + mbpt::op::Λ(2)) * cc.hbar(),
                                      {.connect = cc.hbar_connections()})));
  const auto t = cc.t(), l = cc.λ(), tp = cc.tʼ(1, 1), lp = cc.λʼ(1, 1);
  for (std::size_t k = 1; k <= 2; ++k) {
    const auto rank = L"[" + std::to_wstring(k) + L"]";
    add(L"t" + rank, t[k]);
    add(L"λ" + rank, l[k]);
    add(L"tʼ" + rank, tp[k]);
    add(L"λʼ" + rank, lp[k]);
  }

  result.variables = leaves_as_variables(result.equations);
  return result;
}

}  // namespace sequant::deriv::test

#endif  // SEQUANT_TESTS_UNIT_DERIV_CCSD_HELPERS_HPP
