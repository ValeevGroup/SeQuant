#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/timer.hpp>

namespace sequant::mbpt::sr::so {

class cceqvec {
  size_t N, P, PMIN;
  bool antisymm;

 public:
  cceqvec(size_t n, bool antisymm = true,
          size_t p = std::numeric_limits<size_t>::max(), size_t pmin = 1);
  std::vector<sequant::ExprPtr> operator()(bool screen = true,
                                           bool use_topology = true,
                                           bool use_connectivity = true,
                                           bool canonical_only = true);
};  // class cceqvec

}  // namespace sequant::mbpt::sr::so

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
