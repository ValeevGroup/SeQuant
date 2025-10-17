#ifndef SEQUANT_EXPRESSIONS_EXPR_HPP
#define SEQUANT_EXPRESSIONS_EXPR_HPP

#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <boost/core/demangle.hpp>

#include <atomic>
#include <memory>
#include <optional>

#include <range/v3/all.hpp>

namespace sequant {

/// @brief the wchar used for labeling adjoints, i.e. the superscript + sign
static const wchar_t adjoint_label = L'\u207A';

/// @brief Base expression class

/// Expr represents the interface needed to form expression trees. Classes that
/// represent expressions should publicly derive from this class. Each Expr on a
/// tree has links to its children Expr objects. The lifetime of Expr objects is
/// expected to be managed by std::shared_ptr . Expr is an Iterable over
/// subexpressions (each of which is an Expr itself). More precisely, Expr meets
/// the SizedIterable concept (see
/// https://raw.githubusercontent.com/ericniebler/range-v3/master/doc/std/D4128.md).
/// Specifically, iterators to subexpressions
/// dereference to ExprPtr. Since Expr is a range, it provides begin/end/etc.
/// that can participate in overloads
///       with other functions in the derived class. Consider a Container class
///       derived from a BaseContainer class:
/// @code
///   template <typename T> class Container : public BaseContainer, public Expr
///   {
///     // WARNING: BaseContainer::begin clashes with Expr::begin
///     // WARNING: BaseContainer::end clashes with Expr::end
///     // etc.
///   };
/// @endcode
/// There are two possible scenarios:
///   - if @c Container is a container of Expr objects, BaseContainer will
///   iterate over ExprPtr objects already
///     and both ranges will be equivalent; it is sufficient to add `using
///     BaseContainer::begin`, etc. to Container's public API.
///   - if @c Container is a container of non-Expr objects, iteration over
///   BaseContainer is likely to be more commonly used
///     in practice, hence again adding `using BaseContainer::begin`, etc. will
///     suffice. To be able to iterate over subexpression range (in this case it
///     is empty) Expr provides Expr::expr member to cast to Expr:
/// @code
///    Container c(...);
///    for(const auto& e: c) {  // iterates over elements of BaseContainer
///    }
///    for(const auto& e: c.expr()) {  // iterates over subexpressions
///    }
/// @endcode
class Expr : public std::enable_shared_from_this<Expr>,
             public ranges::view_facade<Expr> {
 public:
  using range_type = ranges::view_facade<Expr>;
  using hash_type = std::size_t;
  using type_id_type = int;  // to speed up comparisons

  Expr() = default;
  virtual ~Expr() = default;

  /// @return true if this is a leaf
  bool is_atom() const { return ranges::empty(*this); }

  /// @return the string representation of @c this in the LaTeX format
  virtual std::wstring to_latex() const;

  /// @return the string representation of @c this in the Wolfram Language
  /// format
  virtual std::wstring to_wolfram() const;

  /// @return a clone of this object, i.e. an object that is equal to @c this
  /// @note - must be overridden in the derived class.
  ///       - the default implementation throws an exception
  virtual ExprPtr clone() const;

  /// like Expr::shared_from_this, but returns ExprPtr
  /// @return a shared_ptr to this object wrapped into ExprPtr
  /// @throw std::bad_weak_ptr if this object is not managed by a shared_ptr
  ExprPtr exprptr_from_this() {
    return static_cast<ExprPtr>(this->shared_from_this());
  }

  /// like Expr::shared_from_this, but returns ExprPtr
  /// @return a shared_ptr to this object wrapped into ExprPtr
  /// @throw std::bad_weak_ptr if this object is not managed by a shared_ptr
  ExprPtr exprptr_from_this() const {
    return static_cast<const ExprPtr>(
        std::const_pointer_cast<Expr>(this->shared_from_this()));
  }

  /// Canonicalizes @c this and returns the byproduct of canonicalization (e.g.
  /// phase)
  /// @return the byproduct of canonicalization, or @c nullptr if no byproduct
  /// generated
  virtual ExprPtr canonicalize(
      CanonicalizeOptions = CanonicalizeOptions::default_options()) {
    return {};  // by default do nothing and return nullptr
  }

  /// Performs approximate, but fast, canonicalization of @c this and returns
  /// the byproduct of canonicalization (e.g. phase) The default is to use
  /// canonicalize(), unless overridden in the derived class.
  /// @return the byproduct of canonicalization, or @c nullptr if no byproduct
  /// generated
  virtual ExprPtr rapid_canonicalize(
      CanonicalizeOptions = CanonicalizeOptions::default_options().copy_and_set(
          CanonicalizationMethod::Rapid)) {
    return this->canonicalize({.method = CanonicalizationMethod::Rapid});
  }

  // clang-format off
  /// recursively visit this expression, i.e. call visitor on each subexpression
  /// in depth-first fashion.
  /// @warning this will only work for tree expressions; no checking is
  /// performed that each subexpression has only been visited once
  /// TODO make work for graphs
  /// @tparam Visitor a callable of type void(ExprPtr&) or void(const ExprPtr&)
  /// @param visitor the visitor object
  /// @param atoms_only if true, will visit only the leaves; the default is to
  /// visit all nodes
  /// @return true if this object was visited
  /// @sa expr_range
  // clang-format on
  template <typename Visitor>
  bool visit(Visitor &&visitor, const bool atoms_only = false) {
    return visit_impl(*this, std::forward<Visitor>(visitor), atoms_only);
  }

  /// const version of visit
  template <typename Visitor>
  bool visit(Visitor &&visitor, const bool atoms_only = false) const {
    return visit_impl(*this, std::forward<Visitor>(visitor), atoms_only);
  }

  auto begin_subexpr() { return range_type::begin(); }

  auto end_subexpr() { return range_type::end(); }

  auto begin_subexpr() const { return range_type::begin(); }

  auto end_subexpr() const { return range_type::end(); }

  Expr &expr() { return *this; }
  const Expr &expr() const { return *this; }

  template <typename T, typename Enabler = void>
  struct is_shared_ptr_of_expr : std::false_type {};
  template <typename T>
  struct is_shared_ptr_of_expr<std::shared_ptr<T>,
                               std::enable_if_t<is_expr_v<T>>>
      : std::true_type {};
  template <typename T, typename Enabler = void>
  struct is_shared_ptr_of_expr_or_derived : std::false_type {};
  template <typename T>
  struct is_shared_ptr_of_expr_or_derived<std::shared_ptr<T>,
                                          std::enable_if_t<is_an_expr_v<T>>>
      : std::true_type {};

  /// @brief Reports if this is a c-number
  /// (https://en.wikipedia.org/wiki/C-number), i.e. it commutes
  /// multiplicatively with c-numbers and q-numbers
  /// @return true if this is a c-number
  /// @warning this returns true for all leaves, hence must be overridden for
  /// leaf q-numbers
  /// @note for leaves this has O(1) cost, for non-leaves this involves checking
  /// subexpressions
  virtual bool is_cnumber() const {
    if (is_atom())
      return true;
    else {
      bool result = true;
      for (auto it = begin_subexpr(); result && it != end_subexpr(); ++it) {
        result &= (*it)->is_cnumber();
      }
      return result;
    }
  }

  /// @brief Checks if this commutes (wrt multiplication) with @c that
  /// @return true if this commutes with @c that
  /// @note the default implementation checks if either is c-number; if both are
  /// q-numbers
  ///       this checks commutativity of each subexpression with (each
  ///       subexpression of) @c that
  /// @note expressions are assumed to always commute with respect to additions
  /// since
  ///       it does not appear that +nonabelian near-rings
  ///       (https://en.wikipedia.org/wiki/Near-ring) are commonly needed.
  /// @note commutativity of leaves is checked by commutes_with_atom()
  bool commutes_with(const Expr &that) const {
    auto this_is_atom = is_atom();
    auto that_is_atom = that.is_atom();
    bool result = true;
    if (this_is_atom && that_is_atom) {
      result =
          this->is_cnumber() || that.is_cnumber() || commutes_with_atom(that);
    } else if (this_is_atom) {
      if (!this->is_cnumber()) {
        for (auto it = that.begin_subexpr(); result && it != that.end_subexpr();
             ++it) {
          result &= this->commutes_with(**it);
        }
      }
    } else {
      for (auto it = this->begin_subexpr(); result && it != this->end_subexpr();
           ++it) {
        result &= (*it)->commutes_with(that);
      }
    }
    return result;
  }

  /// @brief changes this to its adjoint
  /// @note base implementation throws, must be reimplemented in the derived
  /// class
  virtual void adjoint();

  /// Computes and returns the hash value. If default @p hasher is used then the
  /// value will be memoized, otherwise @p hasher will be used to compute the
  /// hash every time.
  /// @param hasher the hasher object, if omitted, the default is used (@sa
  /// Expr::memoizing_hash )
  /// @note always returns 0 unless this derived class overrides
  /// Expr::memoizing_hash
  /// @return the hash value for this Expr
  hash_type hash_value(
      std::function<hash_type(const std::shared_ptr<const Expr> &)> hasher = {})
      const {
    return hasher ? hasher(shared_from_this()) : memoizing_hash();
  }

  /// Computes and returns the derived type identifier
  /// @note this function must be overridden in the derived class
  /// @sa Expr::get_type_id
  /// @return the hash value for this Expr
  virtual type_id_type type_id() const
#if __GNUG__
  {
    abort();
  }
#else
      = 0;
#endif

  friend inline bool operator==(const Expr &a, const Expr &b);

  /// @tparam T Expr or a class derived from Expr
  /// @return true if @c *this is less than @c that
  /// @note the derived class must implement Expr::static_less_than
  template <typename T, typename = std::enable_if<is_an_expr_v<T>>>
  bool operator<(const T &that) const {
    if (type_id() ==
        that.type_id()) {  // if same type, use generic (or type-specific, if
                           // available) comparison
      return static_less_than(static_cast<const Expr &>(that));
    } else {  // order types by type id
      return type_id() < that.type_id();
    }
  }

  /// @return (unique) type id of class T
  template <typename T>
  static type_id_type get_type_id() {
    return type_id_accessor<T>();
  }

  /// sets (unique) type id of class T
  /// @param id the value of type id of class T
  /// @note since get_type_id does not check for duplicates, it's user's
  /// responsiblity to make sure that there are no collisions between type ids
  template <typename T>
  static void set_type_id(type_id_type id) {
    type_id_accessor<T>() = id;
  }

  /// @tparam T an Expr type
  /// @return true if this object is of type @c T
  template <typename T>
  bool is() const {
    if constexpr (is_expr_v<T>)
      return true;
    else if constexpr (std::is_base_of_v<Expr, T>)
      return this->type_id() == get_type_id<meta::remove_cvref_t<T>>();
    else
      return dynamic_cast<const T *>(this) != nullptr;
  }

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  const T &as() const {
    SEQUANT_ASSERT(this->is<T>());
    if constexpr (std::is_base_of_v<Expr, T>) {
      return static_cast<const T &>(*this);
    } else
      return dynamic_cast<const T &>(*this);
  }

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  T &as() {
    SEQUANT_ASSERT(this->is<T>());
    if constexpr (std::is_base_of_v<Expr, T>) {
      return static_cast<T &>(*this);
    } else
      return dynamic_cast<T &>(*this);
  }

  /// @return the (demangled) name of this type
  /// @note uses RTTI
  std::string type_name() const {
    return boost::core::demangle(typeid(*this).name());
  }

  /** @name in-place arithmetic operators
   *  Virtual in-place arithmetic operators to be overridden in expressions for
   * which these make sense.
   */
  ///@{

  /// @brief in-place multiply @c *this by @c that
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator*=(const Expr &that);

  /// @brief in-place non-commutatively-multiply @c *this by @c that
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator^=(const Expr &that);

  /// @brief in-place add @c that to @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator+=(const Expr &that);

  /// @brief in-place subtract @c that from @c *this
  /// @return reference to @c *this
  /// @throw std::logic_error if not implemented for this class, or cannot be
  /// implemented for the particular @c that
  virtual Expr &operator-=(const Expr &that);

  ///@}

 private:
  friend ranges::range_access;

  template <typename E, typename Visitor,
            typename =
                std::enable_if_t<std::is_same_v<meta::remove_cvref_t<E>, Expr>>>
  static bool visit_impl(E &&expr, Visitor &&visitor, const bool atoms_only) {
    if (expr.weak_from_this().use_count() == 0)
      throw std::invalid_argument(
          "Expr::visit: cannot visit expressions not managed by shared_ptr");
    for (auto &subexpr_ptr : expr.expr()) {
      const auto subexpr_is_an_atom = subexpr_ptr->is_atom();
      const auto need_to_visit_subexpr = !atoms_only || subexpr_is_an_atom;
      bool visited = false;
      if (!subexpr_is_an_atom)  // if not a leaf, recur into it
        visited = visit_impl(*subexpr_ptr, std::forward<Visitor>(visitor),
                             atoms_only);
      // call on the subexpression itself, if not yet done so
      if (need_to_visit_subexpr && !visited) visitor(subexpr_ptr);
    }
    // N.B. can only visit itself if visitor is nonmutating!
    bool this_visited = false;
    if constexpr (std::is_invocable_r_v<void, std::remove_reference_t<Visitor>,
                                        const ExprPtr &>) {
      if (!atoms_only || expr.is_atom()) {
        const ExprPtr this_exprptr = expr.exprptr_from_this();
        visitor(this_exprptr);
        this_visited = true;
      }
    }
    return this_visited;
  }

 protected:
  Expr(Expr &&) = default;
  Expr(const Expr &) = default;
  Expr &operator=(Expr &&) = default;
  Expr &operator=(const Expr &) = default;

  struct cursor {
    using value_type = ExprPtr;

    cursor() = default;
    constexpr explicit cursor(ExprPtr *subexpr_ptr) noexcept
        : ptr_{subexpr_ptr} {}
    /// when take const ptr note runtime const flag
    constexpr explicit cursor(const ExprPtr *subexpr_ptr) noexcept
        : ptr_{const_cast<ExprPtr *>(subexpr_ptr)}, const_{true} {}
    bool equal(const cursor &that) const { return ptr_ == that.ptr_; }
    void next() { ++ptr_; }
    void prev() { --ptr_; }
    // TODO figure out why can't return const here if want to be able to assign
    // to *begin(Expr&)
    ExprPtr &read() const {
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    ExprPtr &read() {
      RANGES_EXPECT(const_ == false);
      RANGES_EXPECT(ptr_);
      return *ptr_;
    }
    void assign(const ExprPtr &that_ptr) {
      RANGES_EXPECT(ptr_);
      *ptr_ = that_ptr;
    }
    std::ptrdiff_t distance_to(cursor const &that) const {
      return that.ptr_ - ptr_;
    }
    void advance(std::ptrdiff_t n) { ptr_ += n; }

   private:
    ExprPtr *ptr_ =
        nullptr;  // both begin and end will be represented by this, so Expr
                  // without subexpressions begin() equals end() automatically
    bool const_ = false;  // assert in nonconst ops
  };

  /// @return the cursor for the beginning of the range (must override in a
  /// derived Expr that has subexpressions)
  virtual cursor begin_cursor() { return cursor{}; }
  /// @return the cursor for the end of the range (must override in a derived
  /// Expr that has subexpressions)
  virtual cursor end_cursor() { return cursor{}; }

  /// @return the cursor for the beginning of the range (must override in a
  /// derived Expr that has subexpressions)
  virtual cursor begin_cursor() const { return cursor{}; }
  /// @return the cursor for the end of the range (must override in a derived
  /// Expr that has subexpressions)
  virtual cursor end_cursor() const { return cursor{}; }

  mutable std::optional<hash_type> hash_value_;  // not initialized by default
  virtual hash_type memoizing_hash() const {
    static const hash_type default_hash_value = 0;
    if (hash_value_)
      return *hash_value_;
    else
      return default_hash_value;
  }
  virtual void reset_hash_value() const { hash_value_.reset(); }

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be
  /// statically cast
  /// @return true if @c that is equivalent to *this
  virtual bool static_equal([[maybe_unused]] const Expr &that) const
#if __GNUG__
  {
    abort();
  }
#else
      = 0;
#endif

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be
  /// statically cast
  /// @note base comparison compares Expr::hash_value() , specialize to each
  /// type as needed
  /// @return true if @c *this is less than @c that
  virtual bool static_less_than(const Expr &that) const {
    return this->hash_value() < that.hash_value();
  }

  /// @param that an Expr object
  /// @note @c *this and @c that are guaranteed to be leaves, and neither is a
  /// c-number , hence honest checking is needed
  /// @return true if @c *this multiplicatively commutes with @c that
  /// @note this returns true unless overridden in derived class
  virtual bool commutes_with_atom([[maybe_unused]] const Expr &that) const {
    return true;
  }

 private:
  /// @return returns next type id in the grand class list
  static type_id_type get_next_type_id() {
    static std::atomic<type_id_type> grand_type_id = 0;
    return ++grand_type_id;
  }

  /// sets (unique) type id of class T
  /// @param id the value of type id of class T
  template <typename T>
  static type_id_type &type_id_accessor() {
    static type_id_type type_id = get_next_type_id();
    return type_id;
  }

 private:
  /// @input[in] fn the name of function that is missing in this class
  /// @return a std::logic_error object containing a message describing that @p
  /// fn is missing from this type
  std::logic_error not_implemented(const char *fn) const;
};  // class Expr

template <>
struct Expr::is_shared_ptr_of_expr<ExprPtr, void> : std::true_type {};
template <>
struct Expr::is_shared_ptr_of_expr_or_derived<ExprPtr, void> : std::true_type {
};

/// @return true if @c a is equal to @c b
inline bool operator==(const Expr &a, const Expr &b) {
  if (a.type_id() != b.type_id())
    return false;
  else
    return a.static_equal(b);
}

#if __cplusplus < 202002L
/// @return true if @c a is not equal to @c b
inline bool operator!=(const Expr &a, const Expr &b) { return !(a == b); }
#endif  // __cplusplus < 202002L

/// binary predicate that returns true is 2 expressions differ by a factor
struct proportional_to {
  /// @param[in] expr1
  /// @param[in] expr2
  /// @return true if @p expr1 is proportional to @p expr2
  bool operator()(const ExprPtr &expr1, const ExprPtr &expr2) const;
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_HPP
