//
// Created by Bimal Gaudel on 7/18/21.
//

#include "scf.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace sequant::eval {

bool SequantEvalScf::converged() const { return converged_; }

double SequantEvalScf::energy() const { return energy_; }

//
// todo separate logging functions
//
void SequantEvalScf::scf(std::basic_ostream<wchar_t>& log) {
  using HRC = std::chrono::high_resolution_clock;

  size_t const iter_print_width =
      std::to_wstring(info_.scf_opts.max_iter).length();

  size_t const norm_precision = 10;
  size_t const norm_print_width = double_print_width;

  static std::wstring const header_str = [iter_print_width]() {
    auto oss = std::wostringstream{};

    oss << std::setw(iter_print_width) << std::right << "Iter"
        << "  " << std::setw(double_print_width) << std::left << "Energy"
        << "  " << std::setw(norm_print_width) << std::left << "Diff(Energy)"
        << "  " << std::setw(norm_print_width) << std::left << "Diff(Norm)"
        << "  " << std::setw(time_print_width) << std::left
        << "Time(micro sec)";
    auto const dash_count = oss.str().length();
    oss << "\n" << std::wstring(dash_count, L'-');
    return oss.str();
  }();

  auto log_iter = [&log, iter_print_width](size_t iter, double energy,
                                           double ediff, double norm_diff,
                                           auto time_beg, auto time_end) {
    log << std::setw(iter_print_width) << std::right << iter << "  "
        << std::setprecision(double_precision) << std::setw(double_print_width)
        << std::left << energy << "  " << std::setprecision(norm_precision)
        << std::setw(norm_print_width) << std::left << std::scientific << ediff
        << "  " << std::setprecision(norm_precision)
        << std::setw(norm_print_width) << std::left << std::scientific
        << norm_diff << "  " << std::setw(time_print_width) << std::right
        << std::fixed
        << std::chrono::duration_cast<std::chrono::microseconds>(time_end -
                                                                 time_beg)
               .count()
        << std::endl;
  };

  if (info_.log_opts.level > 0) log << header_str << std::endl;

  double ediff = 0;
  double norm_diff = 0;
  size_t iter = 0;
  do {
    auto t_beg_iter = HRC::now();

    ++iter;
    reset_cache();

    auto norm_last = norm();
    auto energy_last = energy_;

    energy_ = solve();

    norm_diff = std::fabs(norm_last - norm());
    ediff = std::fabs(energy_last - energy_);

    converged_ =
        !(norm_diff > info_.scf_opts.conv || ediff > info_.scf_opts.conv);

    auto t_end_iter = HRC::now();
    if (info_.log_opts.level > 0)
      log_iter(iter, energy_, ediff, norm_diff, t_beg_iter, t_end_iter);
  } while (!converged_ && iter < info_.scf_opts.max_iter);
  //  reset_cache_all();
  if (info_.log_opts.level > 0)
    log << std::wstring(20, L'-') << "\n"
        << "Delta(CC) = " << std::setprecision(double_precision)
        << std::setw(double_print_width) << energy_ << "\n"
        << "Num iters = " << iter << "\n"
        << std::endl;
}

}  // namespace sequant::eval
