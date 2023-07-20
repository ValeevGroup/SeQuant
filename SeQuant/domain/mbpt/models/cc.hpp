#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/timer.hpp>

namespace sequant::mbpt::sr {

/// derives equations of traditional coupled-cluster method
class cceqs {
  size_t N, P, PMIN;

 public:
  cceqs(size_t n, size_t p = std::numeric_limits<size_t>::max(),
        size_t pmin = 1);

  /// derives t amplitude equations
  std::vector<sequant::ExprPtr> t(bool screen = true, bool use_topology = true,
                                  bool use_connectivity = true,
                                  bool canonical_only = true);
  /// derives λ amplitude equations
  std::vector<sequant::ExprPtr> lambda(bool screen = false,
                                       bool use_topology = true,
                                       bool use_connectivity = true,
                                       bool canonical_only = true);

  std::vector<sequant::ExprPtr> t_pert();
  std::vector<sequant::ExprPtr> lambda_pert();
};  // class cceqs

}  // namespace sequant::mbpt::sr

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
