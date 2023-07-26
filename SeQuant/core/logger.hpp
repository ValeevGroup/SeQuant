//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_LOGGER_HPP
#define SEQUANT_LOGGER_HPP

#include "../core/utility/singleton.hpp"

#include <iostream>

namespace sequant {

/// controls logging within SeQuant components, only useful for
/// troubleshooting/learning
struct Logger : public Singleton<Logger> {
  bool wick_harness = false;
  bool wick_topology = false;
  bool wick_contract = false;
  bool wick_reduce = false;
  bool wick_stats = false;
  bool expand = false;
  bool canonicalize = false;
  bool canonicalize_dot = false;
  bool simplify = false;
  bool tensor_network = false;
  size_t log_level_eval = 1;
  std::ostream& stream = std::cout;

 private:
  friend class Singleton<Logger>;
  Logger(int log_level = 0) {
    if (log_level > 0) {
      wick_topology = true;
      wick_contract = true;
      wick_reduce = true;
      wick_stats = true;
      expand = true;
      canonicalize = true;
      simplify = true;
      tensor_network = true;
    }
  }
};

template <typename... Args>
void write_log(Logger& l, Args const&... args) noexcept {
  ((l.stream << args), ...);
}

}  // namespace sequant

#endif  // SEQUANT_LOGGER_HPP
