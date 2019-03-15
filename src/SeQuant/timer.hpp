//
// Created by Eduard Valeyev on 3/28/18.
//

#ifndef SEQUANT_TIMER_HPP
#define SEQUANT_TIMER_HPP

#include <chrono>
#include <iostream>

namespace sequant {

/// TimerPool aggregates \c N C++11 "timers"; used to profile
/// stages of a computation with high resolution
/// @tparam N the number of timers
/// @note member functions are not reentrant, use one Timers object per thread
template<size_t N = 1>
class TimerPool {
 public:
  typedef std::chrono::duration<double> dur_t;
  typedef std::chrono::high_resolution_clock clock_t;
  typedef std::chrono::time_point<clock_t> time_point_t;

  TimerPool() {
    clear();
    set_now_overhead(0);
  }

  /// returns the current time point
  static time_point_t now() { return clock_t::now(); }

  /// use this to report the overhead of now() call; if set, the reported
  /// timings will be adjusted for this overhead
  /// @param ns overhead in nanoseconds
  /// @note this is clearly compiler and system dependent, please measure
  /// carefully (turn off turboboost, etc.)
  void set_now_overhead(size_t ns) { overhead_ = std::chrono::nanoseconds(ns); }

  /// starts timer \c t
  void start(size_t t = 0) { tstart_[t] = now(); }
  /// stops timer \c t
  /// @return the duration, corrected for overhead, elapsed since the last call
  /// to \c start(t)
  dur_t stop(size_t t = 0) {
    const auto tstop = now();
    const dur_t result = (tstop - tstart_[t]) - overhead_;
    timers_[t] += result;
    return result;
  }
  /// reads value (in seconds) of timer \c t , converted to \c double
  double read(size_t t = 0) const { return timers_[t].count(); }
  /// resets timers to zero
  void clear() {
    for (auto t = 0; t != ntimers; ++t) {
      timers_[t] = dur_t::zero();
      tstart_[t] = time_point_t();
    }
  }

 private:
  constexpr static auto ntimers = N;
  dur_t timers_[ntimers];
  time_point_t tstart_[ntimers];
  dur_t overhead_;  // the duration of now() call ... use this to automatically
  // adjust reported timings is you need fine-grained timing
};

}   // namespace sequant

#define SEQUANT_PROFILE_SINGLE(id, call)                                    \
    {                                                                        \
    sequant::TimerPool<> timer;                                             \
    timer.start();                                                           \
    { (call); }                                                              \
    timer.stop();                                                            \
    auto elapsed_seconds = timer.read();                                     \
    std::wcout << id ": elapsed_time = " << std::scientific                  \
               << elapsed_seconds << " seconds" << std::endl;                \
    }


#endif //SEQUANT_TIMER_HPP
