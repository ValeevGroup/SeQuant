//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_LOGGER_HPP
#define SEQUANT_LOGGER_HPP

#include <SeQuant/core/utility/singleton.hpp>

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
  bool canonicalize_input_graph = false;
  bool canonicalize_dot = false;
  bool simplify = false;
  bool tensor_network = false;

  struct {
    ///
    /// Evaluation log verbosity level
    ///   0: No log.
    ///   1: Log what is being evaluated in the evaluation tree (independent
    ///      of tensor algebra backend (TA/BTAS etc.)
    ///   2: Also log invocation of tensor algebra backend within sequant.
    ///   3: Also log TiledArray memory use.
    ///
    size_t level;

    /// the stream for logging; can be set to nullptr
    std::ostream* stream;
  } eval = {0, nullptr};

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
  if (l.eval.stream) {
    ((*l.eval.stream << args), ...);
    (*l.eval.stream).flush();
  }
}

}  // namespace sequant

#endif  // SEQUANT_LOGGER_HPP
