//
// Created by Eduard Valeyev on 3/24/18.
//

#ifndef SEQUANT2_ITERATOR_HPP
#define SEQUANT2_ITERATOR_HPP

#include <range/v3/all.hpp>

namespace sequant2 {

/// flattened_rangenest is a flattened view over a nest of ranges
/// @tparam RangeNest the type of a nest of ranges; currently only a range of ranges is supported
template <typename RangeNest>
class flattened_rangenest : public ranges::view_facade<flattened_rangenest<RangeNest>> {
 public:
  flattened_rangenest() = default;

  explicit flattened_rangenest(RangeNest* r) : range_(r) {}

 private:
  using Range = typename RangeNest::value_type;

  RangeNest* range_;

  friend ranges::range_access;

  /// cursor (=iterator?) type
  struct cursor
  {
   private:
    RangeNest* range_;
    typename RangeNest::iterator range_iter_;  // sequence iterator pointing to the current element's range in the sequence
    typename Range::iterator elem_iter_;  // iterator pointing to the current element
    mutable int64_t elem_index_ = -1;  // index of the current element within the sequence

    void compute_elem_index() const {
      // accumulate all elements before this range_iter_
      elem_index_ = std::accumulate(begin(*range_), range_iter_, 0, [](std::size_t v, const Range& r) { return v + std::size(r); });
      // accumulate all elements before this elem_iter_
      if (range_iter_ != end(*range_))
        elem_index_ += elem_iter_ - begin(*range_iter_);
    }

   public:
    /// constructs an uninitialized cursor
    cursor() = default;
    /// constructs a cursor pointing to the begin, if range is not empty
    /// @note has O(1) complexity
    cursor(RangeNest* range)
        : range_(range), range_iter_(begin(*range_)), elem_iter_(!ranges::empty(*range_) ? begin(*range_iter_) : decltype(elem_iter_){}), elem_index_{!ranges::empty(*range_) ? 0 : -1}
    {}
    /// constructs a cursor pointing to the end
    /// @note has O(1) complexity
    cursor(RangeNest* range, ranges::default_sentinel) : range_(range), range_iter_(end(*range_)), decltype(elem_iter_){} {}

    /// constructs a cursor pointing to particular @c range_iter and @c elem_iter
    /// @note has O(N) complexity due to the need to compute @c elem_index_
    cursor(RangeNest* range, typename RangeNest::iterator range_iter, typename Range::iterator elem_iter)
        : range_(range), range_iter_(range_iter), elem_iter_(elem_iter)
    {
      compute_elem_index();
    }
    const auto& read() const {
      return *elem_iter_;
    }
    bool equal(const cursor& that) const {
      if (ranges::empty(*range_) || range_iter_ == end(*range_))
        return true;
      else
        return elem_iter_ == that.elem_iter_;
    }
    void next() {
      ++elem_index_;
      ++elem_iter_;
      if (elem_iter_ == end(*range_iter_)) {
        ++range_iter_;
        if (range_iter_ != end(*range_))
          elem_iter_ = begin(*range_iter_);
      }
    }
    const auto index() const {
      if (elem_index_ < 0)
        compute_elem_index();
      return elem_index_;
    }
  };
  cursor begin_cursor() const
  {
    return {range_};
  }
  cursor end_cursor() const
  {
    return {range_, end(*range_), typename Range::iterator{}};
  }
};

}

#endif //SEQUANT2_ITERATOR_HPP
