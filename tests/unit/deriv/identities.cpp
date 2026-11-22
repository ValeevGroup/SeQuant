#include "identities.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/series.hpp>
#include <catch2/catch_test_macros.hpp>
#include "ccsd_helpers.hpp"
#include "csv_test_utils.hpp"
#include "judge.hpp"

#include <algorithm>
#include <utility>

namespace sequant::deriv::test {

namespace {

using V = container::svector<ExprPtr>;

ExprPtr sum_of(Expressions terms) {
  if (terms.empty()) return ex<Constant>(0);
  if (terms.size() == 1) return terms.front();
  return ex<Sum>(terms.begin(), terms.end());
}

/// whether @p leaf holds of a leaf of @p e
bool holds(const ExprPtr& e, const std::function<bool(const Expr&)>& leaf) {
  bool result = false;
  e->visit([&](const ExprPtr& l) { result = result || leaf(*l); },
           /*atoms_only=*/true);
  return result;
}

bool is_tensor_labeled(const Expr& e, std::wstring_view label) {
  return e.is<Tensor>() && e.as<Tensor>().label() == label;
}

/// the most leaves of the symbol @p key a term of @p f holds
std::size_t degree_in(const ExprPtr& f, const std::wstring& key) {
  std::size_t result = 0;
  for (const auto& term : summands_of(f)) {
    std::size_t n = 0;
    term->visit(
        [&](const ExprPtr& l) {
          n += (l->is<Tensor>() || l->is<Variable>()) && symbol_of(*l) == key;
        },
        /*atoms_only=*/true);
    result = std::max(result, n);
  }
  return result;
}

/// @p e without its terms that hold a leaf labelled one of @p labels
ExprPtr without_terms_of(const ExprPtr& e,
                         std::initializer_list<std::wstring_view> labels) {
  Expressions kept;
  for (const auto& term : summands_of(e))
    if (!holds(term, [&](const Expr& l) {
          return std::ranges::any_of(labels, [&](std::wstring_view label) {
            return is_tensor_labeled(l, label) ||
                   (l.is<Variable>() && l.as<Variable>().label() == label);
          });
        }))
      kept.push_back(term);
  return sum_of(std::move(kept));
}

/// @p e with each leaf of the block of @p target replaced by @p value at its
/// indices (core's replace)
ExprPtr with_block_replaced(const ExprPtr& e, const ExprPtr& target,
                            const ExprPtr& value) {
  auto result = e->clone();
  replace(result, ExprMatcher(*target, {.tensor_cmp = TensorComparison::Block}),
          *value);
  return result;
}

/// the first factor Â or Ŝ of a term of @p e ; none if no term has one
const AbstractTensor* projector_of(const ExprPtr& e) {
  for (const auto& term : summands_of(e)) {
    Constant::scalar_type scalar = 1;
    for (const auto& f : factors_of(term, scalar))
      if (is_projector(*f)) return &f->as<AbstractTensor>();
  }
  return nullptr;
}

/// a leaf of @p e labelled @p label whose bra has the spaces of @p indices ,
/// in any order
const Tensor* leaf_over(const ExprPtr& e, std::wstring_view label,
                        const container::svector<Index>& indices) {
  auto spaces = [](const auto& range) {
    container::svector<std::wstring> result;
    for (const auto& i : range) result.emplace_back(i.space().base_key());
    std::ranges::sort(result);
    return result;
  };
  const Tensor* result = nullptr;
  e->visit(
      [&](const ExprPtr& l) {
        if (result || !is_tensor_labeled(*l, label)) return;
        const auto& t = l->as<Tensor>();
        if (spaces(t.bra()) == spaces(indices)) result = &t;
      },
      /*atoms_only=*/true);
  return result;
}

/// the occupied and the virtual free indices of the first product of @p e ,
/// each ordered by spin, then label: SeQuant writes the ↑ column of a
/// mixed-spin amplitude first
std::pair<container::svector<Index>, container::svector<Index>> free_indices(
    const ExprPtr& e) {
  const auto& registry = *get_default_context().index_basis_registry();
  container::svector<Index> occupied, virtuals;
  for (const auto& i :
       indices_of_product(products_of(summands_of(e).front()).front()).free)
    (registry.is_pure_occupied(i.space()) ? occupied : virtuals).push_back(i);
  for (auto* is : {&occupied, &virtuals})
    std::ranges::sort(*is, {}, [](const Index& i) {
      return std::pair(i.space().qns().to_int32(), std::wstring(i.label()));
    });
  return {std::move(occupied), std::move(virtuals)};
}

/// SeQuant's biorthogonal transform of the closed-shell moment @p r over the
/// external pairs (ket, bra) of the projector @p p , as
/// closed_shell_CC_spintrace_v2 applies it
ExprPtr biorthogonal(const ExprPtr& r, const AbstractTensor& p) {
  auto moment = r->clone();
  pop_tensor(moment, reserved::symm_label());
  container::svector<container::svector<Index>> pairs;
  for (std::size_t k = 0; k != p._bra_rank(); ++k)
    pairs.push_back({p._ket()[k], p._bra()[k]});
  return mbpt::biorthogonal_transform_pre_nnsproject(moment, pairs);
}

/// @p y spelled apart from x, as schwarz() differentiates by it
ExprPtr schwarz_spelling(const ExprPtr& y) { return renamed_apart(y, 20); }

Family family_of(const Model& model) {
  Index::reset_tmp_index();
  const auto generated = model.open_shell
                             ? open_shell_equations()
                             : generated_equations(model.csv, model.traced);
  Family family{.model = model, .equations = generated.equations};
  family.variables = leaves_as_variables(generated.equations);
  container::set<std::wstring> spellings;
  for (const auto& x : family.variables) spellings.insert(serialize(x));
  for (const auto& x : generated.variables)
    if (spellings.insert(serialize(x)).second) family.variables.push_back(x);
  return family;
}

}  // namespace

std::string name_of(const Model& model) {
  if (model.open_shell)
    return std::string("CC{2} ") + (model.csv ? "CSV " : "") + "open-shell";
  return std::string("CC{2} ") + (model.csv ? "CSV" : "canonical") +
         (model.traced ? " traced" : " spin-orbital");
}

container::svector<Model> family_models() {
#ifdef SEQUANT_SKIP_LONG_TESTS
  return {{}, {.csv = true}};
#else
  return {{},
          {.traced = true},
          {.csv = true},
          {.csv = true, .traced = true},
          {.open_shell = true}};
#endif
}

void with_family(const Model& model,
                 const std::function<void(const Family&)>& f) {
  sequant::tests::csv::ScopedCsvContext scoped{
      mbpt::make_legacy_registry(), model.csv ? mbpt::CSV::Yes : mbpt::CSV::No};
  static container::map<std::string, Family> made;
  auto it = made.find(name_of(model));
  if (it == made.end())
    it = made.emplace(name_of(model), family_of(model)).first;
  f(it->second);
}

std::uint64_t stable_hash(std::string_view text) {
  std::uint64_t h = 0xcbf29ce484222325ULL;
  for (const auto c : text) {
    h ^= static_cast<unsigned char>(c);
    h *= 0x100000001b3ULL;
  }
  return h;
}

ExprPtr equation_named(const Family& family, std::wstring_view name) {
  const auto it = std::ranges::find_if(
      family.equations, [&](const auto& e) { return e.first == name; });
  INFO(toUtf8(name));
  REQUIRE(it != family.equations.end());
  return it->second;
}

bool is_amplitude(const ExprPtr& x) {
  if (!x->is<Tensor>()) return false;
  const auto label = mbpt::detail::strip_pert_order(
      strip_adjoint_label(x->as<Tensor>().label()));
  return label == L"t" || label == L"λ";
}

CheckedVariables checked_variables(const Family& family) {
  static const std::string rejection = "deriv::diff: ";
  CheckedVariables result;
  Expressions others;
  for (const auto& x : family.variables) {
    std::string rejected;
    for (const auto& [name, e] : family.equations) {
      try {
        (void)diff(e, {x});
      } catch (const Exception& error) {
        const std::string what = error.what();
        if (!what.starts_with(rejection)) throw;
        rejected = what.substr(rejection.size());
        break;
      }
    }
    if (!rejected.empty())
      result.rejected[rejected].push_back(x);
    else
      (is_amplitude(x) ? result.variables : others).push_back(x);
  }
  std::ranges::sort(others, {},
                    [](const ExprPtr& x) { return stable_hash(text_of(x)); });
  others.resize(std::min(others.size(), sampled_variables));
  result.variables.insert(result.variables.end(), others.begin(), others.end());
  return result;
}

ExprPtr relabeled(const ExprPtr& x, std::wstring_view label) {
  const auto& t = x->as<Tensor>();
  if (is_adjoint_label(t.label()))
    return sequant::adjoint(relabeled(sequant::adjoint(x), label));
  return ex<Tensor>(label, t.bra(), t.ket(), t.aux(), t.symmetry(),
                    t.braket_symmetry(), t.column_symmetry());
}

ExprPtr renamed_apart(const ExprPtr& x, std::size_t by) {
  return transform_expr(
      x, names_for(indices_of_leaf(x),
                   [by](const Index& i) -> std::optional<Index> {
                     return Index(i.basis(), i.ordinal().value() + by);
                   }));
}

ExprPtr without_projectors_over(const ExprPtr& e,
                                const container::set<Index>& indices) {
  Expressions terms;
  for (const auto& term : summands_of(e)) {
    Constant::scalar_type scalar = 1;
    Expressions kept;
    for (const auto& f : factors_of(term, scalar)) {
      const bool over_indices =
          f->is<AbstractTensor>() && is_projector(*f) &&
          std::ranges::all_of(
              f->as<AbstractTensor>()._slots(),
              [&](const Index& i) { return indices.contains(i); });
      if (!over_indices) kept.push_back(f);
    }
    terms.push_back(ex<Product>(scalar, kept, Product::Flatten::No));
  }
  return sum_of(std::move(terms));
}

Comparison definition(const ExprPtr& f, const ExprPtr& x, const ExprPtr& r,
                      std::uint64_t seed) {
  auto result =
      compare_tangents(f, r, x, {}, {.count = samples_per_pair, .seed = seed});
  auto unread = [&] {
    return result.difference.empty() && result.nonzero == 0 && !r->is_zero();
  };
  // a nonzero r of which no nonzero element was drawn: drawn again over three
  // occupied values, where three distinct occupied labels meet; then compared
  // at every element
  if (unread())
    result = compare_tangents(
        f, r, x, {.occupied = 3},
        {.count = samples_per_pair, .seed = seed, .draws = retry_draws});
  if (unread()) result = compare_tangents(f, r, x);
  return result;
}

Comparison taylor(const ExprPtr& f, const ExprPtr& x, std::uint64_t seed,
                  std::size_t order) {
  const auto h = relabeled(x, L"ζ");
  const auto x_adjoint = sequant::adjoint(x);
  // where x⁺ holds the elements of x, the derivative by x holds them too
  const auto& t = x->as<Tensor>();
  const bool own_adjoint = t.braket_symmetry() == BraKetSymmetry::Symm ||
                           t.base_field() == Field::Real;
  const auto shifted = with_block_replaced(with_block_replaced(f, x, h),
                                           x_adjoint, sequant::adjoint(h));
  const auto n_max = degree_in(f, symbol_of(*x)),
             m_max = own_adjoint ? 0 : degree_in(f, symbol_of(*x_adjoint));
  const auto block = order == 0 ? block_order(x) : order;
  Expressions series{at_zero(f, {x})};
  for (std::size_t n = 0; n <= n_max; ++n)
    for (std::size_t m = 0; m <= m_max; ++m) {
      if (n + m == 0) continue;
      V xs;
      Expressions hs;
      container::set<Index> indices;
      // each derivative by a spelling of its own
      for (std::size_t k = 0; k != n + m; ++k) {
        const auto v = renamed_apart(k < n ? x : x_adjoint, 10 * (k + 1));
        xs.push_back(v);
        hs.push_back(relabeled(v, L"ζ"));
        indices.merge(indices_of_leaf(v));
      }
      const auto dn = diff(f, xs);
      if (dn->is_zero()) continue;
      // h has the symmetry of x: a projector over x's indices leaves it as it
      // is
      auto term = without_projectors_over(at_zero(dn, {x}), indices);
      if (term->is_zero()) continue;
      for (const auto& hk : hs) term = term * hk;
      sequant::intmax_t denominator = factorial(n) * factorial(m);
      for (std::size_t k = 0; k != n + m; ++k)
        denominator *= static_cast<sequant::intmax_t>(block);
      series.push_back(ex<Constant>(rational{1, denominator}) * term);
    }
  return compare_arrays(shifted, sum_of(std::move(series)), externals_of(f),
                        {.projectors = true},
                        {.count = samples_per_pair, .seed = seed});
}

Comparison schwarz(const ExprPtr& f, const ExprPtr& x, const ExprPtr& y,
                   std::uint64_t seed, Constant::scalar_type scale) {
  const auto z = schwarz_spelling(y);
  auto free = externals_of(f);
  free.merge(indices_of_leaf(x));
  free.merge(indices_of_leaf(z));
  const Sampling sampling{
      .count = samples_per_pair,
      .seed = seed,
      .groups = {externals_of(f), indices_of_leaf(x), indices_of_leaf(z)}};
  return compare_arrays(diff(diff(f, {x}), {z}),
                        ex<Constant>(scale) * diff(diff(f, {z}), {x}), free,
                        {.projectors = true}, sampling);
}

const ExprPtr& schwarz_partner(const ExprPtr& r, const ExprPtr& x,
                               const Expressions& amplitudes) {
  const auto n = amplitudes.size();
  const auto k = static_cast<std::size_t>(
      std::ranges::find_if(amplitudes,
                           [&](const ExprPtr& a) { return *a == *x; }) -
      amplitudes.begin());
  REQUIRE(k != n);
  for (std::size_t o = 1; o < n; ++o) {
    const auto& y = amplitudes[(k + o) % n];
    if (!diff(r, {schwarz_spelling(y)})->is_zero()) return y;
  }
  return amplitudes[(k + 1) % n];
}

container::svector<OwnEquation> own_equations(const Family& family) {
  const auto& model = family.model;
  const auto L = equation_named(family, model.open_shell ? L"ℒ case 0" : L"ℒ");
  container::svector<OwnEquation> result;
  for (const auto& [name, e] : family.equations) {
    const bool t_eq = name.starts_with(L"t["), λ_eq = name.starts_with(L"λ[");
    if ((!t_eq && !λ_eq) || name.starts_with(L"t[0]") ||
        name.starts_with(L"λ[0]"))
      continue;
    // the variable: the other amplitude, spelled with the external indices of
    // the equation, bra and ket as the projector holds them
    container::svector<Index> b, k;
    const auto* projector = projector_of(e);
    if (projector) {
      for (const auto& i : projector->_bra()) b.push_back(i);
      for (const auto& i : projector->_ket()) k.push_back(i);
    } else {
      const auto [occupied, virtuals] = free_indices(e);
      b = t_eq ? occupied : virtuals;
      k = t_eq ? virtuals : occupied;
    }
    const auto label = t_eq ? L"λ" : L"t";
    const auto* leaf = leaf_over(L, label, b);
    REQUIRE(leaf);
    const auto x = mbpt::detail::tensor_like(*leaf, label, b, k, {});
    const auto r = diff(L, {x});
    OwnEquation entry{.name = name, .generated = e, .expected = r};
    if (model.open_shell) {
      // the generated same-spin residuals carry no Â
      if (x->as<Tensor>().symmetry() == Symmetry::Antisymm && b.size() > 1)
        entry.generated =
            ex<Tensor>(reserved::antisymm_label(), bra(b), ket(k), aux{},
                       Symmetry::Antisymm, BraKetSymmetry::Nonsymm,
                       ColumnSymmetry::Symm) *
            e;
    } else if (model.traced) {
      // closed_shell_CC_spintrace_v2 writes the residual after SeQuant's
      // biorthogonal transform
      REQUIRE(projector);
      entry.expected = biorthogonal(r, *projector);
    }
    entry.free = externals_of(entry.generated);
    result.push_back(std::move(entry));
  }
  // the response pair: tʼ[k] (λʼ[k]) is ⟨μ|h̄¹ + H̄ T¹|0⟩ ∓ ω t¹ (λ¹), whose
  // T¹ (Λ¹) terms are d/dε t[k] (λ[k]) at t + ε t¹ (λ + ε λ¹)
  const auto eps = ex<Variable>(L"ε");
  for (const auto& [name, e] : family.equations) {
    const bool t_eq = name.starts_with(L"tʼ"), λ_eq = name.starts_with(L"λʼ");
    if (!t_eq && !λ_eq) continue;
    auto shifted = mbpt::series(
        equation_named(
            family,
            (t_eq ? std::wstring(L"t") : std::wstring(L"λ")) + name.substr(2)),
        L"t", eps, 1);
    if (λ_eq) shifted = mbpt::series(shifted, L"λ", eps, 1);
    auto generated = without_terms_of(e, {L"h¹", L"ω"});
    auto expected = at_zero(diff(shifted, {eps}), {eps});
    if (model.open_shell && !projector_of(e)) {
      // the same-spin residuals carry no Â
      const auto [occupied, virtuals] = free_indices(e);
      const bool same_spin =
          occupied.size() == 2 &&
          occupied[0].space().qns() == occupied[1].space().qns();
      if (same_spin) {
        const auto A = ex<Tensor>(
            reserved::antisymm_label(), bra(t_eq ? occupied : virtuals),
            ket(t_eq ? virtuals : occupied), aux{}, Symmetry::Antisymm,
            BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
        generated = A * generated;
        expected = A * expected;
      }
    }
    result.push_back({.name = name,
                      .generated = generated,
                      .expected = expected,
                      .free = externals_of(e)});
  }
  return result;
}

}  // namespace sequant::deriv::test
