//
// Created by Bimal Gaudel on 7/17/21.
//

#ifndef SEQUANT_EVAL_SCF_HPP
#define SEQUANT_EVAL_SCF_HPP

#include <cstddef>
#include <iostream>

#include "examples/eval/calc_info.hpp"

namespace sequant::eval {

class SequantEvalScf {
 protected:
  CalcInfo info_;

  double energy_ = 0;

  bool converged_ = false;

  inline static size_t const double_precision =
      std::numeric_limits<decltype(energy_)>::max_digits10;

  inline static size_t const double_print_width = double_precision + 5;

  inline static size_t const time_print_width =
      std::to_wstring(std::numeric_limits<std::intmax_t>::max()).length();

  SequantEvalScf(CalcInfo const& srcc_info) : info_{srcc_info} {}

  [[nodiscard]] virtual double norm() const = 0;

  virtual void reset_cache() = 0;

  virtual double solve() = 0;

 public:
  virtual ~SequantEvalScf() = default;

  [[nodiscard]] bool converged() const;

  [[nodiscard]] double energy() const;

  void scf(std::basic_ostream<wchar_t>& log);
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_SCF_HPP
