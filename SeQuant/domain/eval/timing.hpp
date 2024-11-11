//
// Created by Bimal Gaudel on 11/11/24.
//

#ifndef SEQUANT_EVAL_TIMING_HPP
#define SEQUANT_EVAL_TIMING_HPP

#include <SeQuant/core/wstring.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace sequant {

class EvaluationTiming {
 public:
  using Duration = double;

  ///
  /// \brief Clears the timing data.
  ///
  inline void clear() { timing_.clear(); }

  ///
  /// \brief Append a timing info.
  ///
  /// \param label An identifier for an evaluation.
  /// \param time  The timing associated with the evaluation.
  ///
  inline void append(std::string_view label, Duration time) {
    timing_.emplace_back(label.data(), time);
  }

  ///
  /// \brief Append a timing info.
  ///
  /// \param label An identifier for an evaluation.
  /// \param time  The timing associated with the evaluation.
  ///
  inline void append(std::wstring_view label, Duration time) {
    append(sequant::to_string(label), time);
  }

  ///
  /// \brief Sort the timing entries based by comparing durations.
  ///        Default is to put the entries with longer duration first.
  ///
  /// \tparam Cmp A duration comparer. Default: std::greater<Duration>.
  ///
  /// \param cmp Duration comparer value. Default: std::greater<Duration>{}.
  ///
  template <typename Cmp = std::greater<Duration>>
  inline void time_sort(Cmp cmp = Cmp{}) {
    auto sort_fn = [cmp](Entry const& lhs, Entry const& rhs) {
      return cmp(lhs.time(), rhs.time());
    };
    ranges::sort(timing_, sort_fn);
  }

  ///
  /// \return The timing durations.
  ///
  [[nodiscard]] inline auto time() const {
    return timing_ | ranges::views::transform(&Entry::time);
  }

  ///
  /// \return The labels used to construct timings.
  ///
  [[nodiscard]] inline auto label() const {
    return timing_ | ranges::views::transform(&Entry::label);
  }

  ///
  /// \return A zipped view of time() and label().
  ///
  [[nodiscard]] inline auto time_to_label() const {
    return ranges::views::zip(time(), label());
  }

 private:
  class Entry {
    std::string label_;
    Duration time_;

   public:
    Entry(std::string_view l, Duration t) : label_(l.data()), time_(t) {}

    [[nodiscard]] inline std::string const& label() const { return label_; }

    [[nodiscard]] inline Duration time() const { return time_; }
  };

  std::vector<Entry> timing_;
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_TIMING_HPP
