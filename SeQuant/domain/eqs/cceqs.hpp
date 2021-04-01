#ifndef SEQUANT_EQS_CCEQS_HPP
#define SEQUANT_EQS_CCEQS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>

namespace sequant::eqs {

class cceqvec {
  size_t N, P, PMIN;

 public:
  cceqvec(size_t n, size_t p = std::numeric_limits<size_t>::max(),
          size_t pmin = 1);
  std::vector<sequant::ExprPtr> operator()(bool screen, bool use_topology,
                                           bool use_connectivity,
                                           bool canonical_only);
};  // class cceqvec

class compute_cceqvec {
  size_t P, PMIN, N;

 public:
  compute_cceqvec(size_t p, size_t pmin, size_t n);

  void operator()(bool print, bool screen, bool use_topology,
                  bool use_connectivity, bool canonical_only);
};  // class compute_cceqvec

class compute_all {
  size_t NMAX;

 public:
  compute_all(size_t nmax);

  void operator()(bool print = true, bool screen = true,
                  bool use_topology = true, bool use_connectivity = true,
                  bool canonical_only = true);
};  // class compute_all

}  // namespace sequant::eqs
#endif  // SEQUANT_EQS_CCEQS_HPP
