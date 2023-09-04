#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/timer.hpp>

namespace sequant::mbpt::sr {

/// derives equations of traditional coupled-cluster method
class cc {
  size_t N, P, PMIN;

 public:
  cc(size_t n, size_t p = std::numeric_limits<size_t>::max(),
        size_t pmin = 1);

  /// derives similarity-transformed expressions of mpbt::Operators
  /// @param expr expression to be transformed
  /// @param r order of truncation
  /// @return transformed expression
  ExprPtr sim_tr(ExprPtr expr, size_t r);

  /// derives t amplitude equations
  std::vector<sequant::ExprPtr> t(bool screen = true, bool use_topology = true,
                                  bool use_connectivity = true,
                                  bool canonical_only = true);
  /// derives λ amplitude equations
  std::vector<sequant::ExprPtr> λ(bool screen = false, bool use_topology = true,
                                  bool use_connectivity = true,
                                  bool canonical_only = true);

  /// derives first-order perturbed t amplitude equations
  std::vector<sequant::ExprPtr> pert_t1();

  /// derives first-order perturbed λ amplitude equations
  std::vector<sequant::ExprPtr> pert_λ1();
};  // class cc

}  // namespace sequant::mbpt::sr

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
