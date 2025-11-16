#ifndef SEQUANT_EXPRESSIONS_EXPR_RANGE_HPP
#define SEQUANT_EXPRESSIONS_EXPR_RANGE_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/facade.hpp>

namespace sequant {

/// @brief a view of the leaves/atoms of an Expr tree
/// @note this is just like view::join, except fully recursive and its iterator
/// provides not only elements but also their indices as well as the host
/// ranges.
/// @note to traverse all nodes, not just the leaves, use Expr::visit @sa
/// Expr::visit
class ExprRange : public ranges::view_facade<ExprRange> {
 public:
  using base_type = ranges::view_facade<ExprRange>;

  ExprRange() = default;

  explicit ExprRange(ExprPtr top) : top_(std::move(top)) {}

  ExprRange(const ExprRange&) = default;
  ExprRange(ExprRange&&) = default;
  ExprRange& operator=(const ExprRange&) = default;
  ExprRange& operator=(ExprRange&&) = default;

  ExprPtr top() const { return top_; }

 private:
  ExprPtr top_;

  friend ranges::range_access;

  /// the cursor type
  struct cursor {
   private:
    ExprPtr* top_;
    // current element is encoded by a sequence of {pointer, index} pairs to its
    // parents e.g. consider the a -> b -> {c,d} tree; when the cursor points to
    // d address_ will contain {{&a, 0}, {&b, 1}}, i.e. the cursor is pointing
    // in the second element within b, which is in turn is the first element in
    // a
    container::svector<std::pair<ExprPtr*, int64_t>> address_;
    ExprPtr* element_ptr_ =
        nullptr;  // pointer to the element, for a valid cursor should be equal
                  // to address_.back().first[address_.back().second]
    int64_t ordinal_ = -1;  // scalar "index" of the current element within the
                            // sequence, -1 marks the end

    // recursively seek the first atom under
    void next_atom(ExprPtr& expr) {
      if (!expr->is_atom()) {
        address_.push_back(std::make_pair(&expr, 0));
        next_atom((*expr)[0]);
      } else {
        SEQUANT_ASSERT(!address_.empty());
        const auto& parent_plus_child = address_.back();
        element_ptr_ = &((**parent_plus_child.first)[parent_plus_child.second]);
      }
    }

   public:
    /// constructs an uninitialized cursor
    cursor() = default;
    /// constructs a cursor pointing to the begin, if range is not empty
    /// @note has O(d) complexity (where d = tree depth)
    cursor(ExprPtr& top) : top_(&top) {
      // if top is nonnull, initialize the address of the first element
      if (*top_) next_atom(*top_);
      if (element_ptr_)  // if have at least one atom
        ordinal_ = 0;
    }
    /// constructs a cursor pointing to the end
    /// @note has O(1) complexity
    cursor(ExprPtr& top, ranges::default_sentinel_t) : top_(&top) {}

    ExprPtr& read() const {
      SEQUANT_ASSERT(ordinal_ != -1);
      return *element_ptr_;
    }
    bool equal(const cursor& that) const {
      return this->element_ptr_ == that.element_ptr_;
    }
    void next() {
      SEQUANT_ASSERT(element_ptr_);
      // first the next parent with children
      auto* parent_plus_child = &(address_.back());
      while ((size_t)parent_plus_child->second + 1 ==
             ranges::size(**(parent_plus_child->first))) {
        address_.pop_back();     // step up, look for next atom
        if (address_.empty()) {  // we might be done if address is empty (i.e.
                                 // we are back at the top
          break;
        } else {
          parent_plus_child = &(address_.back());
        }
      }

      // update element_ptr_, ordinal, and address_ (if not done yet)
      if (address_.empty()) {
        ordinal_ = -1;
        element_ptr_ = nullptr;
      } else {
        ++parent_plus_child->second;
        next_atom((**parent_plus_child->first)[parent_plus_child->second]);
        ++ordinal_;
      }
    }

    auto address() const { return address_; }
    auto ordinal() const {
      SEQUANT_ASSERT(ordinal_ >= 0);
      return ordinal_;
    }

    //    /// calls erase on the current iterator
    //    void erase() {
    //      SEQUANT_ASSERT(range_iter_ != this->_end(*range_));
    //      // TODO resolve the compilation issue
    //      //      ranges::erase(*range_iter_, elem_iter_);
    //      // verify that capacity does not change
    //      const auto capacity = range_iter_->capacity();
    //      range_iter_->erase(elem_iter_);
    //      SEQUANT_ASSERT(capacity == range_iter_->capacity());
    //    }
    //    /// calls erase on the current iterator
    //    template <typename T> void insert(T &&elem) {
    //      SEQUANT_ASSERT(range_iter_ != this->_end(*range_));
    //      // TODO resolve the compilation issue
    //      //      ranges::insert(*range_iter_, elem_iter_,
    //      std::forward<T>(elem));
    //      // verify that capacity does not change
    //      const auto capacity = range_iter_->capacity();
    //      range_iter_->insert(elem_iter_, std::forward<T>(elem));
    //      SEQUANT_ASSERT(capacity == range_iter_->capacity());
    //    }
  };
  cursor begin_cursor() { return {top_}; }
  cursor end_cursor() { return {top_, ranges::default_sentinel_t{}}; }

  // public:
  //  using iterator = ranges::basic_iterator<cursor>;
};

// For backwards compatibility
using expr_range = ExprRange;

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_RANGE_HPP
