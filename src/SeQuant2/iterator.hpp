//
// Created by Eduard Valeyev on 3/24/18.
//

#ifndef SEQUANT2_ITERATOR_HPP
#define SEQUANT2_ITERATOR_HPP

#include <range/v3/all.hpp>

namespace sequant2 {

/// flattened_rangenest is a flattened view over a nest of ranges
/// @tparam RangeNest the type of a nest of ranges; currently only a range of
/// ranges is supported
template <typename RangeNest>
class flattened_rangenest
    : public ranges::view_facade<flattened_rangenest<RangeNest>> {
public:
  flattened_rangenest() = default;

  explicit flattened_rangenest(RangeNest *r) : range_(r) {}

  flattened_rangenest(const flattened_rangenest&) = default;
  flattened_rangenest(flattened_rangenest&&) = default;
  flattened_rangenest& operator=(const flattened_rangenest&) = default;
  flattened_rangenest& operator=(flattened_rangenest&&) = default;

  RangeNest* range() const{ return range_; }

 private:
  using Range = typename RangeNest::value_type;

  RangeNest *range_;

  friend ranges::range_access;

  /// cursor (=iterator?) type
  struct cursor {
  private:
    RangeNest *range_;
    typename RangeNest::iterator
        range_iter_; // sequence iterator pointing to the current element's
                     // range in the sequence
    typename Range::iterator
        elem_iter_; // iterator pointing to the current element
    mutable int64_t elem_index_ =
        -1; // index of the current element within the sequence

    void compute_elem_index() const {
      // accumulate all elements before this range_iter_
      elem_index_ = std::accumulate(
          _begin(*range_), range_iter_, 0,
          [](std::size_t v, const Range &r) { return v + std::size(r); });
      // accumulate all elements before this elem_iter_
      if (range_iter_ != _end(*range_))
        elem_index_ += elem_iter_ - _begin(*range_iter_);
    }

    template <typename Range> static auto _begin(Range& rng) {
      using std::begin;
      return begin(rng);
    }
    template <typename Range> static auto _end(Range& rng) {
      using std::end;
      return end(rng);
    }

  public:
    /// constructs an uninitialized cursor
    cursor() = default;
    /// constructs a cursor pointing to the begin, if range is not empty
    /// @note has O(1) complexity
    cursor(RangeNest *range)
        : range_(range), range_iter_(find_if(_begin(*range_), _end(*range_), [](const auto& e) { using std::empty; return !empty(e); } )),
          elem_iter_(range_iter_ != _end(*range_) ? _begin(*range_iter_) : decltype(elem_iter_){}),
          elem_index_{range_iter_ != _end(*range_) ? 0 : -1} {}
    /// constructs a cursor pointing to the end
    /// @note has O(1) complexity
    cursor(RangeNest *range, ranges::default_sentinel)
        : range_(range), range_iter_(_end(*range_)) {}

    /// constructs a cursor pointing to particular @c range_iter and @c
    /// elem_iter
    /// @note has O(N) complexity due to the need to compute @c elem_index_
    cursor(RangeNest *range, typename RangeNest::iterator range_iter,
           typename Range::iterator elem_iter)
        : range_(range), range_iter_(range_iter), elem_iter_(elem_iter) {
      compute_elem_index();
    }
    const auto &read() const { return *elem_iter_; }
    bool equal(const cursor &that) const {
      if (range_ == that.range_) {  // make sure these point to same range
        const auto end_range_iter = _end(*range_);
        const auto this_is_the_end = range_iter_ == end_range_iter;
        const auto that_is_the_end = that.range_iter_ == end_range_iter;
        if (this_is_the_end && that_is_the_end)
          return true;
        else if (this_is_the_end || that_is_the_end)
          return false;
        else
          return elem_iter_ == that.elem_iter_;
      }
      else  // this points to a different range from that
        return false;
    }
    void next() {
      ++elem_index_;
      ++elem_iter_;
      if (elem_iter_ == _end(*range_iter_)) {
        ++range_iter_;
        // skip empty ranges
        const auto this_is_the_end = _end(*range_);
        while (range_iter_ != this_is_the_end && ranges::empty(*range_iter_))
          ++range_iter_;
        if (range_iter_ != this_is_the_end)
          elem_iter_ = _begin(*range_iter_);
      }
    }

    const auto range_iter() const { return range_iter_; }
    const auto elem_iter() const { return elem_iter_; }
    const auto index() const {
      if (elem_index_ < 0)
        compute_elem_index();
      return elem_index_;
    }

    /// calls erase on the current iterator
    void erase() {
      assert(range_iter_ != _end(*range_));
      // TODO resolve the compilation issue
      //      ranges::erase(*range_iter_, elem_iter_);
      // verify that capacity does not change
      const auto capacity = range_iter_->capacity();
      range_iter_->erase(elem_iter_);
      assert(capacity == range_iter_->capacity());
    }
    /// calls erase on the current iterator
    template <typename T> void insert(T &&elem) const {
      assert(range_iter_ != _end(*range_));
      // TODO resolve the compilation issue
      //      ranges::insert(*range_iter_, elem_iter_, std::forward<T>(elem));
      // verify that capacity does not change
      const auto capacity = range_iter_->capacity();
      range_iter_->insert(elem_iter_, std::forward<T>(elem));
      assert(capacity == range_iter_->capacity());
    }
  };
  cursor begin_cursor() const { return {range_}; }
  cursor end_cursor() const { return {range_, ranges::default_sentinel{}}; }
};

} // namespace sequant2

#endif // SEQUANT2_ITERATOR_HPP
