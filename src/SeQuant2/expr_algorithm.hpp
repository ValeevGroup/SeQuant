//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT2_EXPR_ALGORITHM_HPP
#define SEQUANT2_EXPR_ALGORITHM_HPP

namespace sequant2 {

/// Recursively canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; will be replaced if canonicalization is impure
inline void canonicalize(ExprPtr& expr) {
  const auto biproduct = expr->canonicalize();
  if (biproduct && biproduct->is<Constant>()) {
    expr = biproduct * expr;
  }
}

namespace detail {
struct simplify_visitor {
  void operator()(ExprPtr &expr) {
    if (debug) std::wcout << "simplify_visitor received " << to_latex(expr) << std::endl;
    // apply simplify() iteratively until done
    while (simplify(expr)) {
      if (debug) std::wcout << "after 1 round of simplification have " << to_latex(expr) << std::endl;
    }
    if (debug) std::wcout << "simplification result = " << to_latex(expr) << std::endl;
    // canonicalization to be done by other visitors
  }

  /// simplifies a Product by:
  /// - flattening subproducts
  /// - factoring in constants
  /// @param[in,out] expr (shared_ptr to ) a Product
  bool simplify_product(ExprPtr &expr) {
    auto &expr_ref = *expr;
    std::shared_ptr<Product> result;
    const auto nsubexpr = ranges::size(*expr);
    for (std::size_t i = 0; i != nsubexpr; ++i) {
      const auto expr_i_is_product = expr_ref[i]->is<Product>();
      const auto expr_i_is_constant = expr_ref[i]->is<Constant>();
      if (expr_i_is_product || expr_i_is_constant) {
        // allocate the result, if not done yet
        if (!result) {
          auto expr_product = std::static_pointer_cast<Product>(expr);
          // copy preceding factors
          result = std::make_shared<Product>(expr_product->scalar(), begin(expr->expr()), begin(expr->expr()) + i);
        }
        // append correctly handles products and constants
        result->append(1, std::move(expr_ref[i]));
      }
    }
    bool expr_changed = false;
    if (result) {
      expr = std::static_pointer_cast<Expr>(result);
      expr_changed = true;
    }
    const auto expr_size = ranges::size(*expr);
    auto expr_product = std::static_pointer_cast<Product>(expr);
    if (expr_product->scalar() == 0.) {  // if scalar = 0, make it 0 (too aggressive?)
      expr = ex<Constant>(0);
      expr_changed = true;
    }
    else if (expr_size == 0) {  // if product reduced to a constant make it a constant
      expr = ex<Constant>(expr_product->scalar());
      expr_changed = true;
    }
    else if (expr_size == 1 && expr_product->scalar() == 1.) {  // if product has 1 term and the scalar is 1, lift the factor
      expr = (*expr)[0];
      expr_changed = true;
    }
    return expr_changed;
  }

  /// simplifies a Sum
  bool simplify_sum(ExprPtr &expr) {
    auto &expr_ref = *expr;
    std::shared_ptr<Sum> result;  // will keep the result if one or more summands is expanded
    const auto nsubexpr = ranges::size(*expr);
    if (result) { // if anything changed
      expr = std::static_pointer_cast<Expr>(result);
      return true;
    } else
      return false;
  }

  // @return true if any simplifications were performed
  bool simplify(ExprPtr &expr) {
    if (expr->is<Product>()) {
      return simplify_product(expr);
    } else if (expr->is<Sum>()) {
      return simplify_sum(expr);
    } else
      return false;
  }

  bool debug = false;
};
};  // namespace detail

/// Simplifies an Expr by applying all known transformations (e.g. eliminating trivial math, CSE, etc.)
/// @param[in,out] expr expression to be simplified
inline void simplify(ExprPtr& expr) {
  detail::simplify_visitor simplifier{};
  expr->visit(simplifier);
  simplifier(expr);
}

namespace detail {
struct expand_visitor {
  void operator()(ExprPtr& expr) {
    if (debug) std::wcout << "expand_visitor received " << to_latex(expr) << std::endl;
    // apply expand() iteratively until done
    while(expand(expr)) {
      if (debug)  std::wcout << "after 1 round of expansion have " << to_latex(expr) << std::endl;
    }
    if (debug) std::wcout << "expansion result = " << to_latex(expr) << std::endl;
    // simplification and canonicalization are to be done by other visitors
  }

  /// expands the first Sum in a Product
  /// @param[in,out] expr (shared_ptr to ) a Product whose first Sum gets expanded; on return @c expr contains the result
  bool expand_product(ExprPtr& expr) {
    auto& expr_ref = *expr;
    std::shared_ptr<Sum> result;
    const auto nsubexpr = ranges::size(*expr);
    for(std::size_t i=0; i != nsubexpr; ++i) {
      if (expr_ref[i]->is<Sum>()) {
        // make template for expr cloning to avoid cloning the Sum we are about to expand
        auto scalar = std::static_pointer_cast<Product>(expr)->scalar();
        auto exprseq_clone_template = container::svector<ExprPtr>(ranges::begin(*expr), ranges::end(*expr));
        exprseq_clone_template[i].reset();
        // allocate the result, if not done yet
        if (!result)
          result = std::make_shared<Sum>();
        ExprPtr subexpr_to_expand = expr_ref[i];
        for(auto& subsubexpr: *subexpr_to_expand) {
          auto exprseq_clone = clone(exprseq_clone_template); // clone the product factors without the expanded sum
          exprseq_clone[i] = subsubexpr;  // scavenging summands here
          using std::begin;
          using std::end;
          result->append(std::move(ex<Product>(scalar, begin(exprseq_clone), end(exprseq_clone))));
        }
        expr = std::static_pointer_cast<Expr>(result); // expanded one Sum, return
        return true;
      }
    }
    return false;
  }

  /// expands a Sum
  bool expand_sum(ExprPtr& expr) {
    auto& expr_ref = *expr;
    std::shared_ptr<Sum> result;  // will keep the result if one or more summands is expanded
    const auto nsubexpr = ranges::size(*expr);
    if (debug) std::wcout << "in expand_sum: expr = " << to_latex(expr) << std::endl;
    for(std::size_t i=0; i != nsubexpr; ++i) {
      if (expr_ref[i]->is<Product>()) {
        const auto this_term_expanded = expand_product(expr_ref[i]);
        // if this is the first term that was expanded, create a result and copy all preceeding subexpressions into it
        if (!result && this_term_expanded) {
          result = std::make_shared<Sum>();
          for(std::size_t j=0; j != i; ++j)
            result->append(expr_ref[j]);
        }
        // if expr != expanded result append current subexpr
        if (result)
          result->append(expr_ref[i]);
        if (debug) std::wcout << "in expand_sum: after expand_product(" << (this_term_expanded ? "true)" : "false)") << " result = " << to_latex(result ? result : expr) << std::endl;
      }
    }
    bool expr_changed = false;
    if (result) { // if any summand was expanded, copy result into expr
      expr = std::static_pointer_cast<Expr>(result);
      expr_changed = true;
    }
    if (ranges::size(*expr) == 1) {  // if sum contains 1 element, raise it
      expr = (*expr)[0];
      expr_changed = true;
    }
    else if (ranges::size(*expr) == 0) {  // if sum contains 0 elements, turn to 0
      expr = ex<Constant>(0);
      expr_changed = true;
    }
    return expr_changed;
  }

  // @return true if expanded Product of Sum into Sum of Product
  bool expand(ExprPtr& expr) {
    if (expr->is<Product>()) {
      return expand_product(expr);
    } else if (expr->is<Sum>()) {
      return expand_sum(expr);
    } else
      return false;
  }

  bool debug = false;
};
};  // namespace detail

/// Recursively expands products of sums
inline void expand(ExprPtr& expr) {
  detail::expand_visitor expander{};
  expr->visit(expander);
  expander(expr);
}

/// @brief a view of the leafs/atoms of an Expr tree
/// @note this is just like view::join, except fully recursive and its iterator provides not only elements but
/// also their indices as well as the host ranges.
/// @note to traverse all nodes, not just the leafs, use Expr::visit @sa Expr::visit
class expr_range
    : public ranges::view_facade<expr_range> {
 public:
  using base_type = ranges::view_facade<expr_range>;

  expr_range() = default;

  explicit expr_range(ExprPtr top) : top_(std::move(top)) {}

  expr_range(const expr_range&) = default;
  expr_range(expr_range&&) = default;
  expr_range& operator=(const expr_range&) = default;
  expr_range& operator=(expr_range&&) = default;

  ExprPtr top() const{ return top_; }

 private:
  ExprPtr top_;

  friend ranges::range_access;

  /// the cursor type
  struct cursor {
   private:
    ExprPtr* top_;
    // current element is encoded by a sequence of {pointer, index} pairs to its parents
    // e.g. consider the a -> b -> {c,d} tree; when the cursor points to d address_ will contain {{&a, 0}, {&b, 1}}, i.e.
    // the cursor is pointing in the second element within b, which is in turn is the first element in a
    container::svector<std::pair<ExprPtr*, int64_t>> address_;
    ExprPtr* element_ptr_ = nullptr;  // pointer to the element, for a valid cursor should be equal to address_.back().first[address_.back().second]
    int64_t ordinal_ =
        -1; // scalar "index" of the current element within the sequence, -1 marks the end

    // recursively seek the first atom under
    void next_atom(ExprPtr& expr) {
      if (!expr->is_atom()) {
        address_.push_back(std::make_pair(&expr,0));
        next_atom((*expr)[0]);
      }
      else {
        assert(!address_.empty());
        const auto& parent_plus_child = address_.back();
        element_ptr_ = &((**parent_plus_child.first)[parent_plus_child.second]);
      }
    }

   public:
    /// constructs an uninitialized cursor
    cursor() = default;
    /// constructs a cursor pointing to the begin, if range is not empty
    /// @note has O(d) complexity (where d = tree depth)
    cursor(ExprPtr& top)
        : top_(&top) {
      // if top is nonnull, initialize the address of the first element
      if (*top_)
        next_atom(*top_);
      if (element_ptr_)  // if have at least one atom
        ordinal_ = 0;
    }
    /// constructs a cursor pointing to the end
    /// @note has O(1) complexity
    cursor(ExprPtr& top, ranges::default_sentinel)
        : top_(&top) {}

    ExprPtr& read() const {
      assert(ordinal_ != -1);
      return *element_ptr_;
    }
    bool equal(const cursor &that) const {
      return this->element_ptr_ == that.element_ptr_;
    }
    void next() {
      assert(element_ptr_);
      // first the next parent with children
      auto* parent_plus_child = &(address_.back());
      while (parent_plus_child->second + 1 == ranges::size(**(parent_plus_child->first))) {
        address_.pop_back();  // step up, look for next atom
        if (address_.empty()) { // we might be done if address is empty (i.e. we are back at the top
          break;
        }
        else {
          parent_plus_child = &(address_.back());
        }
      }

      // update element_ptr_, ordinal, and address_ (if not done yet)
      if (address_.empty()) {
        ordinal_ = -1;
        element_ptr_ = nullptr;
      }
      else {
        ++parent_plus_child->second;
        next_atom((**parent_plus_child->first)[parent_plus_child->second]);
        ++ordinal_;
      }
    }

    const auto address() const { return address_; }
    const auto ordinal() const {
      assert(ordinal_ >= 0);
      return ordinal_;
    }

//    /// calls erase on the current iterator
//    void erase() {
//      assert(range_iter_ != this->_end(*range_));
//      // TODO resolve the compilation issue
//      //      ranges::erase(*range_iter_, elem_iter_);
//      // verify that capacity does not change
//      const auto capacity = range_iter_->capacity();
//      range_iter_->erase(elem_iter_);
//      assert(capacity == range_iter_->capacity());
//    }
//    /// calls erase on the current iterator
//    template <typename T> void insert(T &&elem) {
//      assert(range_iter_ != this->_end(*range_));
//      // TODO resolve the compilation issue
//      //      ranges::insert(*range_iter_, elem_iter_, std::forward<T>(elem));
//      // verify that capacity does not change
//      const auto capacity = range_iter_->capacity();
//      range_iter_->insert(elem_iter_, std::forward<T>(elem));
//      assert(capacity == range_iter_->capacity());
//    }
  };
  cursor begin_cursor() { return {top_}; }
  cursor end_cursor() { return {top_, ranges::default_sentinel{}}; }

// public:
//  using iterator = ranges::basic_iterator<cursor>;
};


}  // namespace sequant2

#endif //SEQUANT2_EXPR_ALGORITHM_HPP
