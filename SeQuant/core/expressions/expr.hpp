#ifndef SEQUANT_EXPRESSIONS_EXPR_HPP
#define SEQUANT_EXPRESSIONS_EXPR_HPP

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/tree_index.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <boost/core/demangle.hpp>

#include <atomic>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <utility>

namespace sequant {

/// the trailing label mark of the adjointed state, `t⁺{...}`
inline constexpr wchar_t adjoint_label = L'⁺';
/// the trailing label mark of the K-conjugated state, `t꙳{...}`, and of a
/// conjugated Variable, `x꙳` (U+A673 SLAVONIC ASTERISK: Unicode has no
/// spacing superscript asterisk; LaTeX renders it as `^{*}`)
inline constexpr wchar_t conjugate_label = L'꙳';

/// @return true if @p label ends with the adjoint marker ::adjoint_label
inline bool is_adjoint_label(std::wstring_view label) {
  return !label.empty() && label.back() == adjoint_label;
}

/// @return @p label without its trailing adjoint marker, if present
inline std::wstring_view strip_adjoint_label(std::wstring_view label) {
  return is_adjoint_label(label) ? label.substr(0, label.size() - 1) : label;
}

/// appends the adjoint marker to @p label, or removes it if already present
inline void toggle_adjoint_label(std::wstring &label) {
  if (is_adjoint_label(label))
    label.pop_back();
  else
    label.push_back(adjoint_label);
}

/// strips the trailing state marks off @p label: ::adjoint_label and
/// ::conjugate_label, in either order and at most one of each
/// @param[in,out] label the name; on return without its marks
/// @return whether an adjoint mark and whether a conjugation mark was found
/// @throw Exception if a mark is repeated
inline std::pair<bool, bool> split_state_marks(std::wstring &label) {
  bool adjointed = false, kconjugated = false;
  while (!label.empty()) {
    const wchar_t c = label.back();
    if (c == adjoint_label) {
      if (adjointed) throw Exception("repeated adjoint mark in the name");
      adjointed = true;
    } else if (c == conjugate_label) {
      if (kconjugated) throw Exception("repeated conjugation mark in the name");
      kconjugated = true;
    } else {
      break;
    }
    label.pop_back();
  }
  return {adjointed, kconjugated};
}

namespace detail {

/// @return the 64-bit FNV-1a hash of @p str
constexpr std::uint64_t fnv1a_64(std::string_view str) {
  std::uint64_t hash = 0xcbf29ce484222325ull;
  for (const char c : str) {
    hash ^= static_cast<unsigned char>(c);
    hash *= 0x100000001b3ull;
  }
  return hash;
}

/// @return the decimal representation of @p n , for composing Expr type names
constexpr std::string uint_to_string(std::uint64_t n) {
  std::string result;
  do {
    result.insert(result.begin(), static_cast<char>('0' + n % 10));
    n /= 10;
  } while (n != 0);
  return result;
}

/// @return the signature of this function as spelled by the compiler, which
/// names @c T ; deterministic for a given compiler, but not portable across
/// compilers
template <typename T>
constexpr std::string compiler_type_name() {
#if defined(_MSC_VER)
  return __FUNCSIG__;
#else
  return __PRETTY_FUNCTION__;
#endif
}

/// the canonical mark of an Expr (see Expr::is_canonical()); copying copies
/// it, moving moves it out of the source
struct CanonicalMark {
  /// identifies the state of the own data of an Expr (not of its
  /// subexpressions) since its last mutation; 0 if not assigned
  std::atomic<std::uint64_t> stamp{0};
  /// if nonzero, the digest of the validity key of the canonicalization that
  /// produced the Expr and of the stamps of its subtree
  std::atomic<std::size_t> seal{0};

  CanonicalMark() = default;
  CanonicalMark(const CanonicalMark &other) noexcept { *this = other; }
  CanonicalMark(CanonicalMark &&other) noexcept { *this = std::move(other); }
  CanonicalMark &operator=(const CanonicalMark &other) noexcept {
    stamp.store(other.stamp.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
    seal.store(other.seal.load(std::memory_order_relaxed),
               std::memory_order_relaxed);
    return *this;
  }
  CanonicalMark &operator=(CanonicalMark &&other) noexcept {
    if (this != &other) {
      *this = static_cast<const CanonicalMark &>(other);
      other.reset();
    }
    return *this;
  }

  void reset() noexcept {
    stamp.store(0, std::memory_order_relaxed);
    seal.store(0, std::memory_order_relaxed);
  }
};

}  // namespace detail

/// @return `T::static_type_name(std::type_identity<T>{})` if @c T declares
/// it, else a name derived from the compiler (see detail::compiler_type_name)
/// @note this is the name used by Expr::get_type_id ; a class template that
/// declares `static_type_name` composes it from the names of its template
/// arguments via this function
/// @note the `std::type_identity<T>` parameter, which does not convert to
/// `std::type_identity` of a base, keeps a derived class from picking up the
/// name its base declares
template <typename T>
constexpr std::string type_name_of() {
  if constexpr (requires {
                  {
                    T::static_type_name(std::type_identity<T>{})
                  } -> std::convertible_to<std::string>;
                })
    return std::string(T::static_type_name(std::type_identity<T>{}));
  else
    return detail::compiler_type_name<T>();
}

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
class Expr : public std::enable_shared_from_this<Expr> {
 public:
  using hash_type = std::size_t;
  using type_id_type = std::uint64_t;
  using type_rank_type = std::uint8_t;

  /// rank of an Expr type that does not declare `type_rank`
  /// @sa Expr::get_type_id
  static constexpr type_rank_type default_type_rank = 128;

  Expr() = default;
  virtual ~Expr() = default;

  /// @return true if this is a leaf
  bool is_atom() const { return empty(); }

  /// @return true if this is zero
  virtual bool is_zero() const { return false; }

  /// @return the string representation of @c this in the LaTeX format
  virtual std::wstring to_latex() const;

  /// @return a clone of this object, i.e. an object that is equal to @c this
  virtual ExprPtr clone() const = 0;

  /// like Expr::shared_from_this, but returns ExprPtr
  /// @return a shared_ptr to this object wrapped into ExprPtr, if this object
  /// is already managed by a shared_ptr, else returns a shared_ptr to a clone
  /// of this object wrapped into ExprPtr
  ExprPtr exprptr_from_this() {
    if (weak_from_this().use_count() == 0)
      return this->clone();
    else
      return static_cast<ExprPtr>(this->shared_from_this());
  }

  /// like Expr::shared_from_this, but returns ExprPtr
  /// @return a shared_ptr to this object wrapped into ExprPtr, if this object
  /// is already managed by a shared_ptr, else returns a shared_ptr to a clone
  /// of this object wrapped into ExprPtr
  ExprPtr exprptr_from_this() const {
    if (weak_from_this().use_count() == 0)
      return this->clone();
    else
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
      CanonicalizeOptions opts =
          CanonicalizeOptions::default_options().copy_and_set(
              CanonicalizationMethod::Rapid)) {
    return this->canonicalize(opts.copy_and_set(CanonicalizationMethod::Rapid));
  }

  /// @param opts the canonicalization options
  /// @return true if this is the full canonical form that canonicalization
  /// with @p opts produced, and since then neither this nor any of its
  /// subexpressions has been mutated or replaced, nor have the contexts in
  /// effect changed (see current_contexts_version())
  /// @note always false if @p opts does not request Complete canonicalization,
  /// the only method whose result is the same for every spelling of an
  /// expression
  bool is_canonical(const CanonicalizeOptions &opts =
                        CanonicalizeOptions::default_options()) const;

  /// records that this, with its subexpressions as they are now, is the full
  /// canonical form for @p opts, so that is_canonical(opts) holds until this
  /// or any of its subexpressions is mutated or replaced, or the contexts in
  /// effect change
  /// @param opts the canonicalization options
  /// @param contexts_version the value of current_contexts_version() when the
  /// canonicalization started; if it has changed since, the mark is not valid
  /// @note no-op if @p opts does not request Complete canonicalization
  /// @warning only to be called with the result of canonicalization with
  /// @p opts: canonicalization leaves an expression marked as canonical alone
  void mark_canonical(
      const CanonicalizeOptions &opts,
      std::uint64_t contexts_version = current_contexts_version()) const;

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

  /// @brief Reports if this is a pure scalar (number-like) expression
  /// @return true if this is a scalar
  /// @note This is distinct from is_cnumber()
  /// @warning this returns false for all leaves by default, hence must be
  /// overridden for scalar leaf types.
  virtual bool is_scalar() const {
    if (is_atom()) return false;
    for (auto it = begin(); it != end(); ++it) {
      if (!(*it)->is_scalar()) return false;
    }
    return true;
  }

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
      for (auto it = begin(); result && it != end(); ++it) {
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
        for (auto it = that.begin(); result && it != that.end(); ++it) {
          result &= this->commutes_with(**it);
        }
      }
    } else {
      for (auto it = this->begin(); result && it != this->end(); ++it) {
        result &= (*it)->commutes_with(that);
      }
    }
    return result;
  }

  /// @brief changes this to its adjoint and returns the sign byproduct: +1,
  /// or −1 when the adjoint is minus the resulting object (an anti-Hermitian
  /// Tensor); like canonicalize()'s byproduct it must be applied by the
  /// caller, see sequant::adjoint(const ExprPtr&)
  [[nodiscard]] virtual std::int8_t adjoint() = 0;

  /// @brief complex conjugation of the represented operator, `K O K⁻¹`, with
  /// `K` complex conjugation in the coordinate representation: on a scalar
  /// the complex conjugate, on a Tensor the K-conjugated state (see
  /// Tensor::kconjugate), on an operator string the identity (the
  /// conjugation acts through the coefficients). Not the conjugate of a
  /// matrix element's value, which is adjoint() with the slots exchanged.
  /// @return the sign to apply, as for adjoint()
  [[nodiscard]] virtual std::int8_t kconjugate() = 0;

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
  /// @sa Expr::get_type_id
  /// @return the hash value for this Expr
  virtual type_id_type type_id() const = 0;

  friend inline bool operator==(const Expr &a, const Expr &b);

  /// @tparam T Expr or a class derived from Expr
  /// @return true if @c *this is less than @c that
  /// @note the derived class must implement Expr::static_less_than
  template <typename T>
    requires(is_an_expr_v<T>)
  bool operator<(const T &that) const {
    if (type_id() ==
        that.type_id()) {  // if same type, use generic (or type-specific, if
                           // available) comparison
      return static_less_than(static_cast<const Expr &>(that));
    } else {  // order types by type id
      return type_id() < that.type_id();
    }
  }

  /// @return the type id of a type of rank @p rank named @p name : @p rank in
  /// the top 8 bits, the top 56 bits of `fnv1a_64(name)` below them
  static constexpr type_id_type make_type_id(type_rank_type rank,
                                             std::string_view name) {
    return (static_cast<type_id_type>(rank) << 56) |
           (detail::fnv1a_64(name) >> 8);
  }

  /// @return the rank encoded in type id @p id
  /// @sa Expr::make_type_id
  static constexpr type_rank_type type_rank_of(type_id_type id) {
    return static_cast<type_rank_type>(id >> 56);
  }

  /// @return the (unique) type id of class T
  /// @details The id is `make_type_id(rank, name)`, where `rank` is
  /// `T::type_rank` (a `static constexpr Expr::type_rank_type`) if
  /// @c T declares it, else Expr::default_type_rank , and `name` is
  /// `sequant::type_name_of<T>()`. Hence Expr::operator< orders unlike types by
  /// rank first. A type whose relative order must not depend on the compiler
  /// declares `static constexpr std::string
  /// static_type_name(std::type_identity<Self> = {})`; `constexpr` makes it
  /// usable in constant expressions, but is not required. A derived class
  /// inherits its base's `type_rank` but not its name, so unless it declares
  /// its own name it gets the compiler-derived one. A type in an unnamed
  /// namespace declares a name if another translation unit may define one of
  /// the same name, since types are told apart by mangled name and such a pair
  /// would silently share an id.
  /// @throw sequant::Exception if a type of another mangled name already has
  /// this id
  template <typename T>
  static type_id_type get_type_id() {
    static const type_id_type id = [] {
      const std::string name = sequant::type_name_of<T>();
      type_rank_type rank = default_type_rank;
      if constexpr (requires { T::type_rank; }) {
        static_assert(std::in_range<type_rank_type>(T::type_rank),
                      "type_rank must fit in Expr::type_rank_type");
        rank = static_cast<type_rank_type>(T::type_rank);
      }
      const type_id_type result = make_type_id(rank, name);
      register_type_id(result, name, typeid(T));
      return result;
    }();
    return id;
  }

  /// @tparam T an Expr type
  /// @return true if this object is of type @c T
  template <typename T>
  bool is() const {
    if constexpr (is_expr_v<T>)
      return true;
    else if constexpr (std::is_base_of_v<Expr, T>)
      return this->type_id() == get_type_id<std::remove_cvref_t<T>>();
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

  ExprIterator begin();
  ExprIterator end();
  ConstExprIterator begin() const;
  ConstExprIterator end() const;
  ConstExprIterator cbegin() const;
  ConstExprIterator cend() const;

  virtual ExprIterator begin_subexpr();
  virtual ExprIterator end_subexpr();
  virtual ConstExprIterator begin_subexpr() const;
  virtual ConstExprIterator end_subexpr() const;

  std::size_t size() const;

  bool empty() const;

  /// unchecked element access
  /// @note the bounds check is only performed if `SEQUANT_ASSERT_ENABLED` is
  ///       #defined; use at() for a bounds check that is always performed
  ExprPtr &operator[](std::size_t idx);
  /// @copydoc operator[](std::size_t)
  const ExprPtr &operator[](std::size_t idx) const;

  /// @return The subexpression identified by the given index
  Expr &operator[](const TreeIndex &idx);
  const Expr &operator[](const TreeIndex &idx) const;

  /// checked element access
  /// @throw Exception if @p idx is not less than size()
  ExprPtr &at(std::size_t idx);
  /// @copydoc at(std::size_t)
  const ExprPtr &at(std::size_t idx) const;

  /// @throw Exception if this Expr is empty (e.g. is an atom)
  ExprPtr &front();
  /// @copydoc front()
  const ExprPtr &front() const;

  /// @throw Exception if this Expr is empty (e.g. is an atom)
  ExprPtr &back();
  /// @copydoc back()
  const ExprPtr &back() const;

 private:
  /// reports an out-of-range access by at()/front()/back()
  [[noreturn]] void throw_out_of_range(std::size_t idx) const;

  template <
      typename E, typename Visitor,
      typename = std::enable_if_t<std::is_same_v<std::remove_cvref_t<E>, Expr>>>
  static bool visit_impl(E &&expr, Visitor &&visitor, const bool atoms_only) {
    if (expr.weak_from_this().use_count() == 0)
      throw Exception(
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

  mutable std::optional<hash_type> hash_value_;  // not initialized by default
  virtual hash_type memoizing_hash() const {
    static const hash_type default_hash_value = 0;
    if (hash_value_)
      return *hash_value_;
    else
      return default_hash_value;
  }
  /// invalidates the memoized hash and the canonical mark
  /// @note to be called by every mutation of this object's own data
  virtual void reset_hash_value() const {
    hash_value_.reset();
    reset_canonical_mark();
  }

  /// invalidates the canonical mark, see is_canonical()
  /// @note to be called by every mutation of this object's own data that does
  /// not call reset_hash_value() (mutations of subexpressions, and their
  /// replacement, are detected without it)
  void reset_canonical_mark() const { canonical_mark_.reset(); }

  /// copies the canonical mark of @p other
  /// @pre the own data of @c *this is identical to that of @p other , and its
  /// subexpressions are copies (or clones) of those of @p other
  void copy_canonical_mark(const Expr &other) const {
    canonical_mark_ = other.canonical_mark_;
  }

  /// @param that an Expr object
  /// @note @c that is guaranteed to be of same type as @c *this, hence can be
  /// statically cast
  /// @return true if @c that is equivalent to *this
  virtual bool static_equal(const Expr &that) const = 0;

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
  /// records that type @p type , named @p name , has type id @p id
  /// @throw sequant::Exception if @p id is already recorded for a type of
  /// another mangled name
  static void register_type_id(type_id_type id, const std::string &name,
                               std::type_index type);
  mutable detail::CanonicalMark canonical_mark_;

  /// assigns a stamp to every node of this subtree that lacks one
  void stamp_canonical_subtree() const;

  /// @return the digest of the stamps of this subtree, or null if any of its
  /// nodes lacks a stamp
  std::optional<std::size_t> canonical_subtree_digest() const;
};  // class Expr

/// ranks (`type_rank`) of SeQuant's own Expr types; Expr::operator< orders
/// unlike types by rank first, so this is their order. Spaced so that a new
/// type can be inserted between them; types of Expr::default_type_rank sort
/// between `power` and `boperator`.
namespace expr_type_rank {
inline constexpr Expr::type_rank_type tensor = 10;
inline constexpr Expr::type_rank_type product = 20;
inline constexpr Expr::type_rank_type constant = 30;
inline constexpr Expr::type_rank_type sum = 40;
inline constexpr Expr::type_rank_type variable = 50;
inline constexpr Expr::type_rank_type power = 60;
inline constexpr Expr::type_rank_type boperator = 250;
inline constexpr Expr::type_rank_type foperator = 251;
inline constexpr Expr::type_rank_type bnoperator = 252;
inline constexpr Expr::type_rank_type fnoperator = 253;
}  // namespace expr_type_rank

static_assert(std::ranges::sized_range<Expr>);
static_assert(std::ranges::bidirectional_range<Expr>);
static_assert(std::ranges::random_access_range<Expr>);

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

/// binary predicate that returns true is 2 expressions differ by a factor
struct proportional_to {
  /// @param[in] expr1
  /// @param[in] expr2
  /// @return true if @p expr1 is proportional to @p expr2
  bool operator()(const ExprPtr &expr1, const ExprPtr &expr2) const;
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_HPP
