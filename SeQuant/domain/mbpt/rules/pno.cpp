//
// Integral projection -- see pno.hpp
//

#include <SeQuant/domain/mbpt/rules/pno.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <algorithm>
#include <cstddef>
#include <string>

namespace sequant::mbpt {

bool is_occupied(Index const& idx) {
  return get_default_context().index_space_registry()->is_pure_occupied(
      idx.space());
}

bool is_ovov_integral(AbstractTensor const& tnsr) {
  if (tnsr._bra_rank() != 2 || tnsr._ket_rank() != 2) return false;
  for (std::size_t k = 0; k != 2; ++k)
    if (is_occupied(tnsr._bra()[k]) == is_occupied(tnsr._ket()[k]))
      return false;
  return true;
}

bool is_oovv_integral(AbstractTensor const& tnsr) {
  if (tnsr._bra_rank() != 2 || tnsr._ket_rank() != 2) return false;
  for (std::size_t k = 0; k != 2; ++k)
    if (is_occupied(tnsr._bra()[k]) != is_occupied(tnsr._ket()[k]))
      return false;
  // homogeneous columns must also disagree, else (oo|oo) or (vv|vv)
  return is_occupied(tnsr._bra()[0]) != is_occupied(tnsr._bra()[1]);
}

namespace {

/// what a term's integrals are rewritten against
struct projection_scope {
  ProjectionOptions const& opts;
  /// the number of amplitude tensors of the term
  std::size_t order = 0;
  /// the residual's external indices (the symmetrizer's slots)
  container::svector<Index> externals = {};
  /// moved leg -> its replacement, shared by sibling summands' integrals
  container::map<Index, Index> moved = {};
  /// the legs the enclosing products' integrals have moved so far: another
  /// integral of the same product holding one moves it to a fresh index
  container::set<Index> moved_here = {};
  /// the flat product being rewritten, null inside a Sum's summands
  Product const* product = nullptr;
};

/// @p sum with each summand replaced by @p f of it; @p sum itself if no
/// summand changed
template <typename F>
ExprPtr projection_per_summand(ExprPtr const& sum, F&& f) {
  container::svector<ExprPtr> summands;
  summands.reserve(sum->size());
  bool changed = false;
  for (auto const& summand : *sum) {
    summands.push_back(f(summand));
    changed = changed || summands.back() != summand;
  }
  return changed ? ex<Sum>(std::move(summands)) : sum;
}

/// @return the cell @p tnsr falls in within a term of order @p order, or
/// `None` if it is neither (ov|ov) nor (oo|vv)
ProjectionTerms projection_cell(AbstractTensor const& tnsr, std::size_t order) {
  if (is_ovov_integral(tnsr))
    return order <= 1 ? ProjectionTerms::ExchangeLinear
                      : ProjectionTerms::ExchangeNonlinear;
  if (is_oovv_integral(tnsr))
    return order <= 1 ? ProjectionTerms::CoulombLinear
                      : ProjectionTerms::CoulombNonlinear;
  return ProjectionTerms::None;
}

/// the distinct occupied indices of @p tnsr , sorted like proto indices
Index::index_vector projection_own_pair(AbstractTensor const& tnsr) {
  Index::index_vector occ;
  for (Index const& idx : tnsr._braket())
    if (is_occupied(idx) && !ranges::contains(occ, idx)) occ.push_back(idx);
  std::stable_sort(occ.begin(), occ.end());
  return occ;
}

/// @return the tensor @p node with its legs moved, times one overlap per moved
/// leg; @p node itself if no leg moves
ExprPtr projection_move_legs(ExprPtr const& node, projection_scope& scope) {
  auto const& opts = scope.opts;
  auto const& tnsr = node->as<AbstractTensor>();
  if (tnsr._label() != opts.integral_label) return node;
  if (tnsr._aux_rank() != 0)
    throw Exception(
        "mbpt::project_integral_domains: the integral " +
        toUtf8(tnsr._label()) +
        " carries auxiliary indices; run the projection BEFORE density "
        "fitting");
  const ProjectionTerms cell = projection_cell(tnsr, scope.order);
  if (!any(cell) || !contains(opts.terms, cell)) return node;

  const bool own_pair = opts.domain == ProjectionDomain::OwnPair;
  const Index::index_vector pair =
      own_pair ? projection_own_pair(tnsr) : Index::index_vector{};
  if (own_pair && pair.size() != 2) return node;
  const IndexBasis::optional_instance cell_instance =
      opts.basis == ProjectionBasis::Integral
          ? IndexBasis::optional_instance{opts.cell_instance(cell)}
          : std::nullopt;

  container::map<Index, Index> replacements;
  container::svector<ExprPtr> overlaps;
  std::size_t slot = 0;
  for (Index const& x : tnsr._braket()) {
    const bool bra_leg = slot++ < tnsr._bra_rank();
    if (is_occupied(x) || ranges::contains(scope.externals, x)) continue;
    auto const& target_pair = own_pair ? pair : x.proto_indices();
    auto const& target_instance = opts.basis == ProjectionBasis::Integral
                                      ? cell_instance
                                      : x.basis().basis_instance();
    auto on_target = [&](Index const& idx) {
      return ranges::equal(idx.proto_indices(), target_pair) &&
             idx.basis().basis_instance() == target_instance;
    };
    if (on_target(x)) continue;
    // sibling summands share a leg's replacement; two integrals of one
    // product do not (x' s{x';x} ... s{x;x''} x''), nor do two targets
    auto it = scope.moved.find(x);
    const bool reuse = it != scope.moved.end() && on_target(it->second) &&
                       !scope.moved_here.contains(x);
    const Index replacement =
        reuse ? it->second
              : Index::make_tmp_index(IndexBasis{x.space(), target_instance},
                                      target_pair,
                                      own_pair || x.symmetric_proto_indices());
    if (it == scope.moved.end()) scope.moved.emplace(x, replacement);
    if (scope.product) scope.moved_here.insert(x);
    replacements.emplace(x, replacement);
    // keeps x' and x each in one bra and one ket
    overlaps.push_back(bra_leg ? make_overlap(x, replacement)
                               : make_overlap(replacement, x));
  }
  if (replacements.empty()) return node;

  auto result = std::make_shared<Product>();
  result->append(1, transform_expr(node, replacements));
  for (auto const& s : overlaps) result->append(1, s);
  return result;
}

/// @p x with every integral of the term @p scope describes rewritten, nested
/// sums and products included; @p x itself if nothing changed
ExprPtr projection_rewrite(ExprPtr const& x, projection_scope& scope) {
  if (x->is<AbstractTensor>()) return projection_move_legs(x, scope);
  Product const* const enclosing = scope.product;
  if (x->is<Sum>()) {
    scope.product = nullptr;
    auto result = projection_per_summand(
        x, [&scope](ExprPtr const& s) { return projection_rewrite(s, scope); });
    scope.product = enclosing;
    return result;
  }
  if (!x->is<Product>()) return x;

  auto const& prod = x->as<Product>();
  scope.product = &prod;
  const auto moved_before = scope.moved_here;
  auto result = std::make_shared<Product>();
  result->scale(prod.scalar());
  bool changed = false;
  for (auto const& factor : prod) {
    auto rewritten = projection_rewrite(factor, scope);
    const bool moved = rewritten != factor;
    changed = changed || moved;
    // a moved tensor's overlaps join the enclosing product
    result->append(1, rewritten,
                   moved && factor->is<AbstractTensor>()
                       ? Product::Flatten::Once
                       : Product::Flatten::No);
  }
  scope.product = enclosing;
  scope.moved_here = moved_before;
  return changed ? result : x;
}

ExprPtr projection_term(ExprPtr const& term, ProjectionOptions const& opts) {
  projection_scope scope{.opts = opts};
  std::size_t rank = 0;
  auto const& registry = *get_default_mbpt_context().op_registry();
  term->visit(
      [&](ExprPtr const& x) {
        if (!x->is<AbstractTensor>()) return;
        auto const& t = x->as<AbstractTensor>();
        if (t._label() == reserved::symm_label() ||
            t._label() == reserved::antisymm_label()) {
          rank = t._bra_rank();
          for (Index const& idx : t._slots()) scope.externals.push_back(idx);
        } else if (is_amplitude_tensor(t, registry))
          ++scope.order;
      },
      /* atoms_only = */ true);
  if (rank != 2) return term;
  return projection_rewrite(term, scope);
}

}  // namespace

ExprPtr project_integral_domains(ExprPtr const& expr,
                                 ProjectionOptions const& opts) {
  if (!expr || !any(opts.terms)) return expr;
  if (opts.basis == ProjectionBasis::Integral && !opts.cell_instance)
    throw Exception(
        "mbpt::project_integral_domains: ProjectionBasis::Integral needs a "
        "cell_instance map; the caller owns what an instance means");
  if (opts.basis == ProjectionBasis::Amplitude &&
      opts.domain == ProjectionDomain::PartnerPair)
    throw Exception(
        "mbpt::project_integral_domains: ProjectionDomain::PartnerPair with "
        "ProjectionBasis::Amplitude moves nothing; pick OwnPair or Integral");
  if (!expr->is<Sum>()) return projection_term(expr, opts);
  return projection_per_summand(expr, [&opts](ExprPtr const& term) {
    return projection_term(term, opts);
  });
}

}  // namespace sequant::mbpt
