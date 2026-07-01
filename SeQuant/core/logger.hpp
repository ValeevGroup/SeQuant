//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_LOGGER_HPP
#define SEQUANT_LOGGER_HPP

#include <SeQuant/core/utility/singleton.hpp>

#include <cstddef>
#include <functional>
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
  bool export_equations = false;

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

    /// Optional reducer for the per-op RSS reported in the eval trace: maps
    /// this rank's local RSS (bytes) to the value to report (e.g. the sum over
    /// all ranks = true total app memory, instead of a misleading single-rank
    /// RSS). Injected by the tensor-algebra backend, which holds the World
    /// needed for a collective reduction. The eval log path runs on EVERY rank
    /// (printing() is level>0, identical across ranks), so a collective reducer
    /// here is matched across ranks. Empty = report this rank's RSS unchanged.
    std::function<std::size_t(std::size_t)> rss_reduce = {};

    /// Optional supplier of an extra, backend-defined suffix appended to each
    /// eval-trace line (after the rss field). Injected by the tensor-algebra
    /// backend to report allocator-level memory that RSS alone cannot
    /// distinguish -- e.g. glibc all-arena in-use vs system bytes, so one can
    /// tell live heap from retained-free heap right at an RSS jump. Runs on
    /// EVERY rank on the eval log path (printing() is level>0, identical across
    /// ranks), so an injected collective reduction here is matched across
    /// ranks. Returns a preformatted, already-reduced string (e.g.
    /// "heap_inuse=...B | heap_sys=...B"); empty function = omit the suffix.
    std::function<std::string()> heap_stats = {};
  } eval = {0, nullptr, {}, {}};

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
      export_equations = true;
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
