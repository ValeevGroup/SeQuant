#ifndef SEQUANT_CORE_IRREPS_PROPAGATION_HPP
#define SEQUANT_CORE_IRREPS_PROPAGATION_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/irreps/algebra.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <range/v3/algorithm/contains.hpp>

#include <concepts>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace sequant {

/// Combines an IRREP algebra with a resolver that determines the IRREPs of
/// leaf tensors. Objects of this type are the context needed to determine the
/// IRREPs of arbitrary expressions via irreps_of.
template <IrrepAlgebra A>
class IrrepContext {
 public:
  using algebra_type = A;
  using irrep_type = typename A::irrep_type;
  using set_type = IrrepSet<irrep_type>;
  using resolver_type = std::function<set_type(const Tensor &)>;
  using variable_resolver_type =
      std::function<std::optional<set_type>(const Variable &)>;
  using parser_type = std::function<irrep_type(std::string_view)>;

  /// @param parser Converts the name of an IRREP into the IRREP. Only needed
  /// for canonical_irrep.
  /// @param variable_resolver Determines the IRREPs of scalar variables.
  /// Variables for which it returns no value (or all, if not given) are
  /// totally symmetric.
  IrrepContext(A algebra, resolver_type resolver, parser_type parser = {},
               variable_resolver_type variable_resolver = {})
      : m_algebra(std::move(algebra)),
        m_resolver(std::move(resolver)),
        m_parser(std::move(parser)),
        m_variable_resolver(std::move(variable_resolver)) {}

  explicit IrrepContext(resolver_type resolver)
    requires std::default_initializable<A>
      : IrrepContext(A{}, std::move(resolver)) {}

  const A &algebra() const { return m_algebra; }
  const resolver_type &resolver() const { return m_resolver; }

  /// @returns The IRREPs of the given (leaf) tensor. Tensors with reserved
  /// labels (Kronecker delta, overlap, (anti)symmetrizers, reference
  /// densities, ...) are considered to be totally symmetric. All others are
  /// handed to the resolver.
  set_type irreps_of(const Tensor &tensor) const {
    if (ranges::contains(reserved::labels(), tensor.label())) {
      return totally_symmetric();
    }

    // Multiplying with the totally symmetric IRREP brings the resolved IRREPs
    // into the algebra's canonical form
    return product(totally_symmetric(), m_resolver(tensor));
  }

  /// @returns The IRREPs of the given variable. Unless the variable resolver
  /// says otherwise, variables are totally symmetric.
  set_type irreps_of(const Variable &variable) const {
    if (m_variable_resolver) {
      if (auto irreps = m_variable_resolver(variable)) {
        return product(totally_symmetric(), irreps.value());
      }
    }

    return totally_symmetric();
  }

  set_type totally_symmetric() const { return {m_algebra.totally_symmetric()}; }

  /// @returns Whether the given set contains nothing but the totally
  /// symmetric IRREP
  bool is_totally_symmetric(const set_type &set) const {
    return set == totally_symmetric();
  }

  /// @returns The IRREP of the given name in the algebra's canonical form
  /// (e.g. with relations between IRREPs applied)
  /// @throws Exception if this context has no parser or if the name does not
  /// refer to exactly one IRREP
  irrep_type canonical_irrep(std::string_view name) const {
    if (!m_parser) {
      throw Exception("IRREP context cannot parse IRREP names");
    }

    set_type set = product(totally_symmetric(), {m_parser(name)});
    if (set.is_any() || set.size() != 1) {
      throw Exception("'" + std::string(name) +
                      "' does not refer to exactly one IRREP");
    }

    return *set.begin();
  }

  set_type empty_set() const { return {}; }

  set_type product(const set_type &lhs, const set_type &rhs) const {
    return sequant::product(m_algebra, lhs, rhs);
  }

  set_type conjugate(const set_type &set) const {
    return sequant::conjugate(m_algebra, set);
  }

 private:
  A m_algebra;
  resolver_type m_resolver;
  parser_type m_parser;
  variable_resolver_type m_variable_resolver;
};

/// Requirements for a context that can be used to propagate IRREPs through
/// expressions (e.g. IrrepContext)
template <typename C>
concept IrrepPropagationContext =
    std::copyable<typename C::set_type> &&
    requires(const C &ctx, const Tensor &tensor, const Variable &variable,
             const typename C::set_type &set,
             typename C::set_type &mutable_set) {
      { ctx.irreps_of(tensor) } -> std::convertible_to<typename C::set_type>;
      { ctx.irreps_of(variable) } -> std::convertible_to<typename C::set_type>;
      { ctx.totally_symmetric() } -> std::convertible_to<typename C::set_type>;
      { ctx.empty_set() } -> std::convertible_to<typename C::set_type>;
      { ctx.product(set, set) } -> std::convertible_to<typename C::set_type>;
      { ctx.conjugate(set) } -> std::convertible_to<typename C::set_type>;
      mutable_set.unite(set);
    };

/// @returns A resolver that determines a tensor's IRREPs by its label.
/// Tensors whose label is not in the map are taken to span any IRREP.
template <typename Irrep>
std::function<IrrepSet<Irrep>(const Tensor &)> resolve_by_label(
    container::map<std::wstring, IrrepSet<Irrep>> map) {
  return [map = std::move(map)](const Tensor &tensor) {
    auto it = map.find(std::wstring(tensor.label()));
    return it == map.end() ? IrrepSet<Irrep>::any() : it->second;
  };
}

/// @returns A resolver that determines a tensor's IRREPs by its label and the
/// spaces of its indices (i.e. by tensor block).
/// Tensors not in the map are taken to span any IRREP.
template <typename Irrep>
std::function<IrrepSet<Irrep>(const Tensor &)> resolve_by_block(
    std::map<Tensor, IrrepSet<Irrep>, TensorBlockLessThanComparator> map) {
  return [map = std::move(map)](const Tensor &tensor) {
    auto it = map.find(tensor);
    return it == map.end() ? IrrepSet<Irrep>::any() : it->second;
  };
}

/// @returns A resolver that determines a tensor's IRREPs by matching it
/// against the given patterns. The first matching pattern wins. Tensors not
/// matching any pattern are taken to span any IRREP.
template <typename Irrep>
std::function<IrrepSet<Irrep>(const Tensor &)> resolve_by_pattern(
    std::vector<std::pair<ExprMatcher, IrrepSet<Irrep>>> patterns) {
  return [patterns = std::move(patterns)](const Tensor &tensor) {
    for (const auto &[matcher, irreps] : patterns) {
      if (matcher == tensor) {
        return irreps;
      }
    }
    return IrrepSet<Irrep>::any();
  };
}

/// The different kinds of nodes relevant for the propagation of IRREPs
enum class IrrepNodeKind { Leaf, Product, Sum, Adjoint };

/// Adapter that makes a given type of (expression) node accessible to the
/// generic IRREP propagation. Specializations provide
/// - kind(node) -> IrrepNodeKind
/// - for_each_child(node, f) calling f(const Node&) for every child
/// - leaf_irreps(node, ctx) -> the IRREPs of a leaf node
template <typename Node>
struct IrrepNodeTraits;

namespace detail {

struct NoOpIrrepCallback {
  void operator()(const auto &, const auto &) const {}
};

/// Generic bottom-up computation of the IRREPs of the given node
/// @param on_node Callback that is invoked as on_node(node, irreps) for every
/// visited node (post-order)
template <typename Node, IrrepPropagationContext Ctx, typename OnNode>
typename Ctx::set_type fold_irreps(const Node &node, const Ctx &ctx,
                                   OnNode &on_node) {
  using Traits = IrrepNodeTraits<Node>;
  using Set = typename Ctx::set_type;

  Set result;
  switch (Traits::kind(node)) {
    case IrrepNodeKind::Leaf:
      result = Traits::leaf_irreps(node, ctx);
      break;
    case IrrepNodeKind::Product:
      result = ctx.totally_symmetric();
      Traits::for_each_child(node, [&](const Node &child) {
        result = ctx.product(result, fold_irreps(child, ctx, on_node));
      });
      break;
    case IrrepNodeKind::Sum:
      result = ctx.empty_set();
      Traits::for_each_child(node, [&](const Node &child) {
        result.unite(fold_irreps(child, ctx, on_node));
      });
      break;
    case IrrepNodeKind::Adjoint:
      result = ctx.empty_set();
      Traits::for_each_child(node, [&](const Node &child) {
        result.unite(ctx.conjugate(fold_irreps(child, ctx, on_node)));
      });
      break;
  }

  on_node(node, result);

  return result;
}

}  // namespace detail

template <>
struct IrrepNodeTraits<Expr> {
  static IrrepNodeKind kind(const Expr &expr) {
    if (expr.is<Product>() && !expr.is_zero()) {
      return IrrepNodeKind::Product;
    }
    if (expr.is<Sum>()) {
      return IrrepNodeKind::Sum;
    }
    return IrrepNodeKind::Leaf;
  }

  template <typename F>
  static void for_each_child(const Expr &expr, F &&f) {
    if (expr.is<Product>()) {
      // Note: the scalar prefactor is irrelevant for symmetry
      for (const ExprPtr &factor : expr.as<Product>().factors()) {
        f(*factor);
      }
    } else {
      SEQUANT_ASSERT(expr.is<Sum>());
      for (const ExprPtr &summand : expr.as<Sum>().summands()) {
        f(*summand);
      }
    }
  }

  template <IrrepPropagationContext Ctx>
  static typename Ctx::set_type leaf_irreps(const Expr &expr, const Ctx &ctx) {
    if (expr.is<Tensor>()) {
      const Tensor &tensor = expr.as<Tensor>();
      if (is_adjoint_label(tensor.label())) {
        // Adjoint of a (non-symmetric) tensor -> resolve the bare tensor
        // (consistent with the explicit Adjoint nodes in evaluation trees)
        Tensor bare(tensor);
        bare.adjoint();
        return ctx.conjugate(ctx.irreps_of(bare));
      }
      return ctx.irreps_of(tensor);
    }
    if (expr.is<Constant>() || expr.is<Product>()) {
      // Only zero-valued products end up here (see kind())
      SEQUANT_ASSERT(expr.is<Constant>() || expr.is_zero());
      return expr.is_zero() ? ctx.empty_set() : ctx.totally_symmetric();
    }
    if (expr.is<Variable>()) {
      const Variable &variable = expr.as<Variable>();
      auto irreps = ctx.irreps_of(variable);
      return variable.conjugated() ? ctx.conjugate(irreps) : irreps;
    }
    if (expr.is<Power>()) {
      const Expr &base = *expr.as<Power>().base();
      if (base.is<Variable>() &&
          ctx.irreps_of(base.as<Variable>()) != ctx.totally_symmetric()) {
        throw Exception(
            "Can't determine IRREPs of a power of a variable that is not "
            "totally symmetric");
      }
      return ctx.totally_symmetric();
    }

    throw Exception("Can't determine IRREPs of expression of type " +
                    expr.type_name());
  }
};

/// @returns The IRREPs spanned by the given expression. Expressions that
/// vanish identically (zero constants or products, empty sums) span no IRREP
/// (empty set).
template <IrrepPropagationContext Ctx>
typename Ctx::set_type irreps_of(const Expr &expr, const Ctx &ctx) {
  detail::NoOpIrrepCallback no_op;
  return detail::fold_irreps<Expr>(expr, ctx, no_op);
}

/// @returns The IRREPs spanned by the given result expression
template <IrrepPropagationContext Ctx>
typename Ctx::set_type irreps_of(const ResultExpr &result, const Ctx &ctx) {
  return irreps_of(*result.expression(), ctx);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_IRREPS_PROPAGATION_HPP
