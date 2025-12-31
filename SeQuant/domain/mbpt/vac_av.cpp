//
// Created by Eduard Valeyev on 8/2/23.
//

#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>

#include <vector>

namespace sequant {
namespace mbpt {
inline namespace op {

ExprPtr expectation_value_impl(
    ExprPtr expr, const OpConnections<std::wstring>& op_connections,
    bool use_topology, bool screen, bool skip_clone, bool full_contractions) {
  // use cloned expr to avoid side effects
  if (!skip_clone) expr = expr->clone();

  auto vac_av_product = [&op_connections, use_topology, screen,
                         full_contractions](ExprPtr expr) {
    SEQUANT_ASSERT(expr.is<Product>());
    // extract scalar and factors
    const auto scalar = expr.as<Product>().scalar();
    auto factors = expr.as<Product>().factors();
    // remove Variable types from the Product temporarily
    auto variables = factors | ranges::views::filter([](const auto& factor) {
                       return factor.template is<Variable>();
                     }) |
                     ranges::to_vector;

    auto factors_filtered = factors |
                            ranges::views::filter([](const auto& factor) {
                              return !(factor.template is<Variable>());
                            }) |
                            ranges::to<decltype(factors)>;
    // construct Product with filtered factors
    auto product =
        ex<Product>(scalar, factors_filtered.begin(), factors_filtered.end());

    // compute connections
    std::vector<std::pair<int, int>> connections;
    {
      container::map<std::wstring, std::vector<int>>
          oplbl2pos;  // maps operator labels to the operator positions in the
                      // product
      int pos = 0;
      bool ops_only = true;
      for (const auto& factor : product.as<Product>()) {
        if (factor.is<op_t>()) {
          const auto& op = factor.as<op_t>();
          const std::wstring op_lbl = std::wstring(op.label());
          const auto it = oplbl2pos.find(op_lbl);
          if (it == oplbl2pos.end()) {  // new label
            oplbl2pos.emplace(op_lbl, std::vector<int>{pos});
          } else {
            it->second.emplace_back(pos);
          }
          ++pos;
        } else if (factor.is<FNOperator>() || factor.is<BNOperator>()) {
          ++pos;  // skip FNOperator and BNOperator
          ops_only = false;
        }
      }

      // if composed of ops only, screen out products with zero VEV
      if (ops_only && screen) {
        if (!can_change_qns(product, qns_t{})) {
          return ex<Constant>(0);
        }
      }

      for (const auto& [op1_lbl, op2_lbl] : op_connections) {
        auto it1 = oplbl2pos.find(op1_lbl);
        auto it2 = oplbl2pos.find(op2_lbl);
        if (it1 == oplbl2pos.end() || it2 == oplbl2pos.end())
          continue;  // one of the op labels is not present in the product
        const auto& [dummy1, op1_indices] = *it1;
        const auto& [dummy2, op2_indices] = *it2;
        for (const auto& op1_idx : op1_indices) {
          for (const auto& op2_idx : op2_indices) {
            if (op1_idx < op2_idx) {  // N.B. connections are directional
              connections.emplace_back(op1_idx, op2_idx);
            }
          }
        }
      }
    }

    // lower to tensor form
    lower_to_tensor_form(product);
    simplify(product);

    // compute expectation value
    // call the tensor-level impl function directly
    auto vev = tensor::expectation_value_impl(product, connections,
                                              use_topology, full_contractions);
    // restore Variable types to the Product
    if (!variables.empty())
      ranges::for_each(variables, [&vev](const auto& var) { vev *= var; });

    return simplify(vev);
  };

  ExprPtr result;
  if (expr.is<Product>()) {
    // expand sums in a product
    if (ranges::any_of(expr.as<Product>().factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      expr = expand(expr);
      simplify(expr);  // condense equivalent terms after expansion
      return expectation_value_impl(expr, op_connections, use_topology, screen,
                                    /* skip_clone = */ true, full_contractions);
    } else
      return vac_av_product(expr);
  } else if (expr.is<Sum>()) {
    result = sequant::transform_sum_expr(
        *expr, [&op_connections, use_topology, screen,
                full_contractions](const auto& op_product) {
          return expectation_value_impl(
              op_product, op_connections, use_topology, screen,
              /* skip_clone = */ true, full_contractions);
        });
    simplify(result);  // combine possible equivalent summands
    return result;
  } else if (expr.is<op_t>()) {
    return ex<Constant>(
        0);  // expectation value of a normal-ordered operator is 0
  } else if (expr.is<Constant>() || expr.is<Variable>()) {
    return expr;  // vacuum is normalized
  }
  throw std::invalid_argument(
      "mbpt::op::detail::expectation_value_impl(expr): unknown expression "
      "type");
}

ExprPtr ref_av(ExprPtr expr, const OpConnections<std::wstring>& op_connections,
               bool use_topology, bool screen, bool skip_clone) {
  auto isr = get_default_context().index_space_registry();
  const bool full_contractions =
      (isr->reference_occupied_space() == isr->vacuum_occupied_space()) ? true
                                                                        : false;
  return expectation_value_impl(expr, op_connections, use_topology, screen,
                                skip_clone, full_contractions);
}

ExprPtr vac_av(ExprPtr expr, const OpConnections<std::wstring>& op_connections,
               bool use_topology, bool screen, bool skip_clone) {
  return expectation_value_impl(expr, op_connections, use_topology, screen,
                                skip_clone,
                                /* full_contractions */ true);
}

}  // namespace op
}  // namespace mbpt
}  // namespace sequant
