#ifndef SEQUANT_EXPRESSIONS_EXPR_PTR_HPP
#define SEQUANT_EXPRESSIONS_EXPR_PTR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/traits.hpp>

#include <memory>

namespace sequant {

/// @brief ExprPtr is a multiple-owner smart pointer to Expr

/// It can be used mostly interchangeably with `std::shared_ptr<Expr>`, but
/// also provides convenient mathematical operators (`+=`, etc.)
class ExprPtr : public std::shared_ptr<Expr> {
 public:
  using base_type = std::shared_ptr<Expr>;
  using base_type::operator->;
  using base_type::base_type;

  ExprPtr() = default;
  ExprPtr(const ExprPtr &) = default;
  ExprPtr(ExprPtr &&) = default;
  template <typename E, typename = std::enable_if_t<
                            std::is_same_v<std::remove_const_t<E>, Expr> ||
                            std::is_base_of_v<Expr, std::remove_const_t<E>>>>
  ExprPtr(const std::shared_ptr<E> &other_sptr) : base_type(other_sptr) {}
  template <typename E, typename = std::enable_if_t<is_an_expr_v<E>>>
  ExprPtr(std::shared_ptr<E> &&other_sptr) : base_type(std::move(other_sptr)) {}
  template <typename E, typename = std::enable_if_t<is_an_expr_v<E>>>
  ExprPtr &operator=(const std::shared_ptr<E> &other_sptr) {
    as_shared_ptr() = other_sptr;
    return *this;
  }
  template <typename E, typename = std::enable_if_t<is_an_expr_v<E>>>
  ExprPtr &operator=(std::shared_ptr<E> &&other_sptr) {
    as_shared_ptr() = std::move(other_sptr);
    return *this;
  }

  ExprPtr &operator=(const ExprPtr &) = default;
  ExprPtr &operator=(ExprPtr &&) = default;

  ~ExprPtr() = default;

  /// @return a copy of this object
  /// @sa Expr::clone()
  [[nodiscard]] ExprPtr clone() const &;
  /// @return a moved copy of this object
  /// @note this object is null after the call
  [[nodiscard]] ExprPtr clone() && noexcept;

  base_type &as_shared_ptr() &;
  const base_type &as_shared_ptr() const &;
  base_type &&as_shared_ptr() &&;

  template <typename E, typename = std::enable_if_t<!is_expr_v<E>>>
  std::shared_ptr<E> as_shared_ptr() const {
    assert(this->is<E>());
    return std::static_pointer_cast<E>(this->as_shared_ptr());
  }

  /// dereference operator
  /// @return non-const lvalue reference to the contained Expr object
  /// @pre `this->operator bool()`
  Expr &operator*() &;

  /// dereference operator
  /// @return const lvalue reference to the contained Expr object
  /// @pre `this->operator bool()`
  const Expr &operator*() const &;

  /// dereference operator
  /// @return non-const rvalue reference to the contained Expr object
  /// @pre `this->operator bool()`
  Expr &&operator*() &&;

  // n.b. equality comparison is deep

  /// ExprPtr is equal to a null pointer if it's uninitialized
  friend inline bool operator==(const ExprPtr &x, std::nullptr_t) {
    return x.get() == nullptr;
  }

  friend bool operator==(const ExprPtr &x, const ExprPtr &y);

  /// in-place addition operator

  /// if this is non-null, adds @c other to the contained expressions, otherwise
  /// will make this point to a clone of @c other (see Expr::clone())
  /// @param other expression to add to this
  /// @return reference to @c *this
  ExprPtr &operator+=(const ExprPtr &other);

  /// in-place subtraction operator

  /// if this is non-null, subtracts @c other from the contained expressions,
  /// otherwise will make this point to the negative of a clone of @c other (see
  /// Expr::clone())
  /// @param other expression to add to this
  /// @return reference to @c *this
  ExprPtr &operator-=(const ExprPtr &);

  /// in-place multiplication operator

  /// if this is non-null, multiplies the contained expressions by @c other ,
  /// otherwise will make this point to a clone of @c other (see Expr::clone())
  /// @param other expression to add to this
  /// @return reference to @c *this
  ExprPtr &operator*=(const ExprPtr &);

  /// @tparam T an Expr type
  /// @return true if this object is of type @c T
  template <typename T>
  bool is() const;

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  const T &as() const;

  /// @tparam T an Expr type
  /// @return this object cast to type @c T
  template <typename T>
  T &as();

  /// @return the range size of the contained expression
  std::size_t size() const;

  std::wstring to_latex() const;
};  // class ExprPtr

/// make an ExprPtr to a new object of type T
/// @tparam T a class derived from Expr
/// @tparam Args a parameter pack type such that T(std::forward<Args>...) is
/// well-formed
/// @param args a parameter pack such that T(args...) is well-formed
template <typename T, typename... Args>
ExprPtr ex(Args &&...args) {
  return std::make_shared<T>(std::forward<Args>(args)...);
}

// this is needed when using std::make_shared<X>({ExprPtr,ExprPtr}), i.e. must
// std::make_shared<X>(ExprPtrList{ExprPtr,ExprPtr})
using ExprPtrList = std::initializer_list<ExprPtr>;

static auto expr_ptr_comparer = [](const auto &ptr1, const auto &ptr2) {
  return *ptr1 == *ptr2;
};

using ExprPtrVector = container::svector<ExprPtr>;

/// @brief computes the adjoint of @p expr
/// @param[in] expr an Expr object
/// @return the adjoint of @p expr
ExprPtr adjoint(const ExprPtr &expr);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_PTR_HPP
