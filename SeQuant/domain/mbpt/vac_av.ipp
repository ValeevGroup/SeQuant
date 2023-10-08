//
// Created by Eduard Valeyev on 8/2/23.
//

#ifndef SEQUANT_DOMAIN_MBPT_VAC_AV_IPP
#define SEQUANT_DOMAIN_MBPT_VAC_AV_IPP

// operator-level vac_av is same for SR and MR

ExprPtr vac_av(
    ExprPtr expr,
    std::vector<std::pair<std::wstring, std::wstring>> op_connections,
    bool skip_clone) {
  // use cloned expr to avoid side effects
  if (!skip_clone) expr = expr->clone();

  auto vac_av_product = [&op_connections](ExprPtr expr) {
    assert(expr.is<Product>());
    // compute connections
    std::vector<std::pair<int, int>> connections;
    {
      std::map<std::wstring, std::vector<int>>
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
        }
      }

      // if composed of ops only, screen out products with zero VEV
      if (ops_only) {
        if (!can_change_qns(expr, qns_t{})) {
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
            using std::min;
            using std::max;
            connections.emplace_back(min(op1_idx, op2_idx),
                                     max(op1_idx, op2_idx));
          }
        }
      }
    }

    // lower to tensor form
    auto lower_to_tensor_form = [](ExprPtr& expr) {
      auto op_lowerer = [](ExprPtr& leaf) {
        if (leaf.is<op_t>()) leaf = leaf.as<op_t>().tensor_form();
      };
      expr->visit(op_lowerer, /* atoms only = */ true);
    };
    lower_to_tensor_form(expr);
    expr = simplify(expr);

    // compute VEV
    auto vev = vac_av(expr, connections, /* use_topology = */ true);
    return simplify(vev); // simplify vev since vac_av does not
  };

  ExprPtr result;
  if (expr.is<Product>()) {
    // expand sums in a product
    if (ranges::any_of(expr.as<Product>().factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      expr = expand(expr);
      simplify(expr);  // condense equivalent terms after expansion
      return vac_av(expr, op_connections, /* skip_clone = */ true);
    } else
      return vac_av_product(expr);
  } else if (expr.is<Sum>()) {
    result = sequant::transform_reduce(
        *expr, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [&op_connections](const auto& op_product) {
          return vac_av(op_product, op_connections, /* skip_clone = */ true);
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
      "mpbt::*::op::vac_av(expr): unknown expression type");
}

#endif  // SEQUANT_DOMAIN_MBPT_VAC_AV_IPP
