//
// Created by Eduard Valeyev on 3/24/18.
//

#ifndef SEQUANT_RANGES_HPP
#define SEQUANT_RANGES_HPP

#include <range/v3/all.hpp>

namespace sequant {

/// @brief a flattened view over a nest of ranges
/// @note this is just like view::join, but its iterator provides not only
/// elements but also their indices as well as the host ranges.
///       this is needed to be able to iterate over pairs of elements while
///       skipping pairs of elements from the same subrange
/// @tparam RangeNest the type of a nest of ranges; use flattened_rangenest for
/// RangeNest if you want recursive flattening.
template <typename RangeNest>
class flattened_rangenest
    : public ranges::view_facade<flattened_rangenest<RangeNest>> {
 public:
  using base_type = ranges::view_facade<flattened_rangenest<RangeNest>>;

  flattened_rangenest() = default;

  explicit flattened_rangenest(RangeNest *r) : range_(r) {}

  flattened_rangenest(const flattened_rangenest &) = default;
  flattened_rangenest(flattened_rangenest &&) = default;
  flattened_rangenest &operator=(const flattened_rangenest &) = default;
  flattened_rangenest &operator=(flattened_rangenest &&) = default;

  RangeNest *range() const { return range_; }

  using value_type = typename RangeNest::value_type::value_type;

 private:
  using Range = typename RangeNest::value_type;

  RangeNest *range_;

  friend ranges::range_access;

  /// the cursor type
  struct cursor {
   private:
    RangeNest *range_;
    std::conditional_t<std::is_const_v<RangeNest>,
                       typename RangeNest::const_iterator,
                       typename RangeNest::iterator>
        range_iter_;  // sequence iterator pointing to the current element's
                      // range in the sequence
    std::conditional_t<std::is_const_v<Range>, typename Range::const_iterator,
                       typename Range::iterator>
        elem_iter_;  // iterator pointing to the current element
    mutable int64_t ordinal_ =
        -1;  // index of the current element within the sequence

    // *range_iter_ produces a const lvalue ref, this return nonconst lvalue ref
    Range &current_range() const { return const_cast<Range &>(*range_iter_); }

    void compute_ordinal() const {
      // accumulate all elements before this range_iter_
      ordinal_ = std::accumulate(
          _begin(*range_), range_iter_, 0,
          [](std::size_t v, const Range &r) { return v + std::size(r); });
      // accumulate all elements before this elem_iter_
      if (range_iter_ != _end(*range_))
        ordinal_ += elem_iter_ - this->_begin(current_range());
    }

    static auto _begin(Range &rng) {
      using std::begin;
      return begin(rng);
    }
    static auto _end(Range &rng) {
      using std::end;
      return end(rng);
    }

    static auto _begin(RangeNest &rng) {
      using std::begin;
      return begin(rng);
    }
    static auto _end(RangeNest &rng) {
      using std::end;
      return end(rng);
    }

   public:
    /// constructs an uninitialized cursor
    cursor() = default;
    /// constructs a cursor pointing to the begin, if range is not empty
    /// @note has O(1) complexity
    cursor(RangeNest *range)
        : range_(range),
          range_iter_(std::find_if(this->_begin(*range), this->_end(*range),
                                   [](const auto &e) {
                                     using std::empty;
                                     return !empty(e);
                                   })),
          elem_iter_(range_iter_ != this->_end(*range)
                         ? this->_begin(current_range())
                         : decltype(elem_iter_){}),
          ordinal_{range_iter_ != this->_end(*range) ? 0 : -1} {}
    /// constructs a cursor pointing to the end
    /// @note has O(1) complexity
    cursor(RangeNest *range, ranges::default_sentinel_t)
        : range_(range), range_iter_(this->_end(*range)) {}

    /// constructs a cursor pointing to particular @c range_iter and @c
    /// elem_iter
    /// @note has O(N) complexity due to the need to compute the ordinal
    cursor(RangeNest *range, typename RangeNest::iterator range_iter,
           typename Range::iterator elem_iter)
        : range_(range), range_iter_(range_iter), elem_iter_(elem_iter) {
      compute_ordinal();
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
      } else  // this points to a different range from that
        return false;
    }
    void next() {
      ++ordinal_;
      ++elem_iter_;
      if (elem_iter_ == this->_end(current_range())) {
        ++range_iter_;
        // skip empty ranges
        const auto this_is_the_end = this->_end(*range_);
        while (range_iter_ != this_is_the_end && ranges::empty(*range_iter_))
          ++range_iter_;
        if (range_iter_ != this_is_the_end)
          elem_iter_ = this->_begin(current_range());
      }
    }

    /// @return the iterator pointing to the range in which this is located
    const auto range_iter() const { return range_iter_; }
    /// @return ordinal index of the range in which this is located
    const auto range_ordinal() const { return range_iter_ - _begin(*range_); }
    /// @return the iterator pointing to the element at which this is located
    const auto elem_iter() const { return elem_iter_; }
    /// @return ordinal index of the element at which this is located (i.e.
    /// distance to the first element of the first range)
    const auto ordinal() const {
      if (ordinal_ < 0) compute_ordinal();
      return ordinal_;
    }

    /// calls erase on the current iterator
    void erase() {
      assert(range_iter_ != this->_end(*range_));
      // TODO resolve the compilation issue
      //      ranges::erase(*range_iter_, elem_iter_);
      // verify that capacity does not change
      const auto capacity = range_iter_->capacity();
      range_iter_->erase(elem_iter_);
      assert(capacity == range_iter_->capacity());
    }
    /// calls insert on the current iterator
    template <typename T>
    void insert(T &&elem) {
      assert(range_iter_ != this->_end(*range_));
      // TODO resolve the compilation issue
      //      ranges::insert(*range_iter_, elem_iter_, std::forward<T>(elem));
      // verify that capacity does not change
      const auto capacity = range_iter_->capacity();
      range_iter_->insert(elem_iter_, std::forward<T>(elem));
      assert(capacity == range_iter_->capacity());
    }
  };
  cursor begin_cursor() const { return {range_}; }
  cursor end_cursor() const { return {range_, ranges::default_sentinel_t{}}; }

 public:
  using iterator = ranges::basic_iterator<cursor>;
  using const_iterator = ranges::basic_iterator<const cursor>;
};

}  // namespace sequant

#endif  // SEQUANT_RANGES_HPP
