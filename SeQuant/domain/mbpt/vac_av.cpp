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

namespace {

/// @brief Lowers operator label pairs to normal-operator position pairs
/// @param[in] oplbl2pos map from operator labels to their positions in the
///            product
/// @param[in] label_pairs the label pairs to lower
/// @return lowered position pairs
template <typename Container>
std::vector<std::pair<int, int>> lower_label_pairs(
    const container::map<std::wstring, std::vector<int>>& oplbl2pos,
    const Container& label_pairs) {
  std::vector<std::pair<int, int>> result;
  for (const auto& [op1_lbl, op2_lbl] : label_pairs) {
    auto it1 = oplbl2pos.find(op1_lbl);
    auto it2 = oplbl2pos.find(op2_lbl);
    if (it1 == oplbl2pos.end() || it2 == oplbl2pos.end())
      continue;  // one of the op labels is not present in the product
    const auto& [dummy1, op1_indices] = *it1;
    const auto& [dummy2, op2_indices] = *it2;
    for (const auto& op1_idx : op1_indices) {
      for (const auto& op2_idx : op2_indices) {
        if (op1_idx < op2_idx) {  // N.B. connections are directional
          result.emplace_back(op1_idx, op2_idx);
        }
      }
    }
  }
  return result;
}
}  // namespace

// fwd declare the tensor level impl function
namespace tensor {
ExprPtr expectation_value_impl(ExprPtr expr, OpConnections<int> connect,
                               OpConnections<int> avoid, bool use_top,
                               bool full_contractions);
}  // namespace tensor

ExprPtr expectation_value_impl(ExprPtr expr,
                               const OpConnections<std::wstring>& connect,
                               const OpConnections<std::wstring>& avoid,
                               bool use_topology, bool screen, bool skip_clone,
                               bool full_contractions) {
  // use cloned expr to avoid side effects
  if (!skip_clone) expr = expr->clone();

  auto vac_av_product = [&connect, &avoid, use_topology, screen,
                         full_contractions](ExprPtr expr) {
    SEQUANT_ASSERT(expr.is<Product>());
    // Note: Non-tensor factors are filtered out at the tensor level (see
    // tensor::expectation_value_impl) before reaching WickTheorem

    // compute connections (both required and avoided)
    OpConnections<int> t_connect;
    OpConnections<int> t_avoid;
    {
      container::map<std::wstring, std::vector<int>>
          oplbl2pos;  // maps operator labels to the operator positions in the
                      // product
      int pos = 0;
      bool ops_only = true;
      for (const auto& factor : expr.as<Product>()) {
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
        } else {
          // any other type
          ops_only = false;
        }
      }

      // if composed of ops only, screen out products with zero VEV
      if (ops_only && screen) {
        if (!can_change_qns(expr, qns_t{})) {
          return ex<Constant>(0);
        }
      }

      // assemble connectivity info for tensor level call
      t_connect = lower_label_pairs(oplbl2pos, connect);
      t_avoid = lower_label_pairs(oplbl2pos, avoid);
    }

    // lower to tensor form
    lower_to_tensor_form(expr);
    simplify(expr);

    // compute expectation value
    // call the tensor-level impl function directly
    auto vev = tensor::expectation_value_impl(expr, t_connect, t_avoid,
                                              use_topology, full_contractions);
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
      return expectation_value_impl(expr, connect, avoid, use_topology, screen,
                                    /* skip_clone = */ true, full_contractions);
    } else
      return vac_av_product(expr);
  } else if (expr.is<Sum>()) {
    result = sequant::transform_sum_expr(
        *expr, [&connect, &avoid, use_topology, screen,
                full_contractions](const auto& op_product) {
          return expectation_value_impl(
              op_product, connect, avoid, use_topology, screen,
              /* skip_clone = */ true, full_contractions);
        });
    simplify(result);  // combine possible equivalent summands
    return result;
  } else if (expr.is<op_t>()) {
    return ex<Constant>(
        0);  // expectation value of a normal-ordered operator is 0
  } else if (expr->is_scalar()) {
    return expr;  // vacuum is normalized
  }
  throw Exception(
      "mbpt::op::expectation_value_impl(expr): unknown expression "
      "type");
}

ExprPtr ref_av(ExprPtr expr, EVOptions<std::wstring> opts) {
  auto isr = get_default_context().index_space_registry();
  const bool full_contractions =
      isr->reference_occupied_space() == isr->vacuum_occupied_space();
  return expectation_value_impl(expr, opts.connect, opts.do_not_connect,
                                opts.use_topology, opts.screen, opts.skip_clone,
                                full_contractions);
}

ExprPtr vac_av(ExprPtr expr, EVOptions<std::wstring> opts) {
  return expectation_value_impl(expr, opts.connect, opts.do_not_connect,
                                opts.use_topology, opts.screen, opts.skip_clone,
                                /* full_contractions */ true);
}

}  // namespace op
}  // namespace mbpt
}  // namespace sequant
