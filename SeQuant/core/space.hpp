//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT_SPACE_H
#define SEQUANT_SPACE_H

#include <bitset>
#include <cassert>
#include <cmath>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/wstring.hpp>

#include <range/v3/algorithm/any_of.hpp>

namespace sequant {

class QuantumNumbersAttr;  // used to constrain TypeAttr ctor

class Index;  // friend of TypeAttr

/// @brief TypeAttr denotes the type of index space.
///
/// This class models a host (complete) space partitioned into disjoint
/// subspaces. To simplify implementation of set operations
/// (intersection, union, etc.) it is encoded as a fixed-width (32) bitset.
class TypeAttr {
 public:
  using bitset_t = int32_t;

  /// default ctor creates a null TypeAttr
  constexpr TypeAttr() noexcept = default;

  /// the null TypeAttr
  const static TypeAttr null;

  /// @brief Constructs from a bitset representation

  /// @warning first (most significant) bit is reserved for internal use
  /// @param bitset bitset representation of this Type
  /// @pre `(bitset & make_reserved().bitset) == null.bitset`
  explicit constexpr TypeAttr(bitset_t bitset) noexcept : bitset(bitset) {
    assert((this->bitset & reserved.bitset) == null.bitset);
  }

  /// construct TypeAddr from things that can be cast to bitset_t, but exclude
  /// bool and QuantumNumbersAttr
  template <typename T,
            typename = std::enable_if_t<
                meta::is_statically_castable_v<std::decay_t<T>, bitset_t> &&
                !std::is_same_v<std::decay_t<T>, bool> &&
                !std::is_same_v<std::decay_t<T>, QuantumNumbersAttr> &&
                !std::is_same_v<std::decay_t<T>, TypeAttr>>>
  constexpr TypeAttr(T &&value) noexcept
      : TypeAttr(static_cast<bitset_t>(std::forward<T>(value))) {}

  constexpr explicit operator int64_t() const {
    return static_cast<int64_t>(bitset);
  }
  constexpr explicit operator bitset_t() const { return bitset; }
  constexpr int32_t to_int32() const { return bitset; }

  /// @return true if this object is non-null (i.e. has any bits set)
  constexpr explicit operator bool() const { return bitset != 0; }

  constexpr TypeAttr(const TypeAttr &other) { bitset = other.to_int32(); }

  /// @return union of `*this` and @p other, i.e. `*this` AND @p other
  /// @note equivalent to `this->to_int32() | other.to_int32()`
  constexpr TypeAttr unIon(TypeAttr other) const {
    return TypeAttr(this->to_int32() | other.to_int32());
  }

  /// @return union of @p a and @p b, i.e. @p a AND @p b
  friend constexpr TypeAttr operator|(const TypeAttr a, const TypeAttr b) {
    return a.unIon(b);
  }

  /// @return `*this` XOR @p other
  /// @note equivalent to `this->to_int32() ^ other.to_int32()`
  constexpr const TypeAttr xOr(TypeAttr other) const {
    return TypeAttr(this->to_int32() ^ other.to_int32());
  }

  /// @return @p a XOR @p b
  friend constexpr TypeAttr operator^(const TypeAttr a, const TypeAttr b) {
    return a.xOr(b);
  }

  /// @return intersection of `*this` AND @p other
  /// @note equivalent to `this->to_int32() & other.to_int32()`
  constexpr const TypeAttr intersection(TypeAttr other) const {
    return TypeAttr(this->to_int32() & other.to_int32());
  }

  /// @return intersection of @p a AND @p b
  friend constexpr TypeAttr operator&(const TypeAttr a, const TypeAttr b) {
    return a.intersection(b);
  }

  /// @return complement of `*this`
  /// @note equivalent to `~this->to_int32()`
  constexpr TypeAttr operator~() const { return ~this->to_int32(); }

  friend constexpr bool operator==(const TypeAttr lhs, const TypeAttr rhs) {
    return lhs.to_int32() == rhs.to_int32();
  }
  friend constexpr bool operator!=(const TypeAttr lhs, const TypeAttr rhs) {
    return !(lhs == rhs);
  }

  /// @return true if \c other is included in this object
  constexpr bool includes(TypeAttr other) const {
    return intersection(other) == other;
  }
  /// @return true if in canonical order this object preceeds \c other
  friend constexpr bool operator<(TypeAttr a, TypeAttr b) {
    return a.to_int32() < b.to_int32();
  }

 private:
  bitset_t bitset = 0;

  friend class Index;

  /// first (most significant) bit reserved for creating default space used by
  /// Index that is distinct from the null space
  const static TypeAttr reserved;

  /// makes reserved object
  static TypeAttr make_reserved() {
    TypeAttr result;
    result.bitset = 0x80000000;
    return result;
  }
};  // struct TypeAttr

inline const TypeAttr TypeAttr::null;
inline const TypeAttr TypeAttr::reserved = TypeAttr::make_reserved();

/// denotes other quantum numbers (particle type, spin, etc.)
class QuantumNumbersAttr {
 public:
  using bitset_t = int32_t;

  /// default ctor creates a null QuantumNumbersAttr
  /// @post `static_cast<bool>(*this) == false`
  constexpr QuantumNumbersAttr() noexcept = default;

  /// the null TypeAttr
  const static QuantumNumbersAttr null;

  /// @brief Constructs from a bitset representation

  /// @warning first (most significant) bit is reserved for internal use
  /// @param bitset bitset representation of this Type
  /// @pre `(bitset & make_reserved().bitset()) == null.bitset()`
  explicit constexpr QuantumNumbersAttr(bitset_t bitset) noexcept
      : bitset(bitset) {
    assert((this->bitset & reserved.bitset) == null.bitset);
  }

  template <typename QN,
            typename = std::enable_if_t<
                meta::is_statically_castable_v<std::decay_t<QN>, bitset_t> &&
                !std::is_same_v<std::decay_t<QN>, bool> &&
                !std::is_same_v<std::decay_t<QN>, TypeAttr> &&
                !std::is_same_v<std::decay_t<QN>, QuantumNumbersAttr>>>
  constexpr QuantumNumbersAttr(QN &&value) noexcept
      : bitset(static_cast<bitset_t>(std::forward<QN>(value))) {}

  constexpr explicit operator int64_t() const {
    return static_cast<int64_t>(bitset);
  }
  constexpr explicit operator bitset_t() const { return bitset; }
  constexpr int32_t to_int32() const { return bitset; }

  /// @return true if this object is non-null (i.e. has any bits set)
  constexpr explicit operator bool() const { return bitset != 0; }

  /// @return `*this` XOR @p other
  /// @note equivalent to `this->to_int32() ^ other.to_int32()`
  constexpr QuantumNumbersAttr xOr(QuantumNumbersAttr other) const {
    return QuantumNumbersAttr(this->to_int32() ^ other.to_int32());
  }

  /// @return @p a XOR @p b
  friend constexpr QuantumNumbersAttr operator^(const QuantumNumbersAttr a,
                                                const QuantumNumbersAttr b) {
    return a.xOr(b);
  }

  /// @return union of `*this` and @p other, i.e. `*this` AND @p other
  /// @note equivalent to `this->to_int32() | other.to_int32()`
  constexpr QuantumNumbersAttr unIon(QuantumNumbersAttr other) const {
    return QuantumNumbersAttr(this->to_int32() | other.to_int32());
  }

  /// @return union of @p a and @p b, i.e. @p a AND @p b
  friend constexpr QuantumNumbersAttr operator|(const QuantumNumbersAttr a,
                                                const QuantumNumbersAttr b) {
    return a.unIon(b);
  }

  /// @return intersection of `*this` AND @p other
  /// @note equivalent to `this->to_int32() & other.to_int32()`
  constexpr QuantumNumbersAttr intersection(QuantumNumbersAttr other) const {
    return QuantumNumbersAttr(this->to_int32() & other.to_int32());
  }

  /// @return intersection of @p a AND @p b
  friend constexpr QuantumNumbersAttr operator&(const QuantumNumbersAttr a,
                                                const QuantumNumbersAttr b) {
    return a.intersection(b);
  }

  /// @return complement of `*this`
  /// @note equivalent to `~this->to_int32()`
  constexpr QuantumNumbersAttr operator~() const { return ~this->to_int32(); }

  friend constexpr bool operator==(QuantumNumbersAttr lhs,
                                   QuantumNumbersAttr rhs) {
    return lhs.to_int32() == rhs.to_int32();
  }
  friend constexpr bool operator!=(QuantumNumbersAttr lhs,
                                   QuantumNumbersAttr rhs) {
    return !(lhs == rhs);
  }

  /// @return true if \c other is included in this object
  bool includes(QuantumNumbersAttr other) const {
    return intersection(other) == other;
  }
  /// @return true if in canonical order this object preceeds \c other
  friend constexpr bool operator<(QuantumNumbersAttr lhs,
                                  QuantumNumbersAttr rhs) {
    return lhs.to_int32() < rhs.to_int32();
  }

 private:
  bitset_t bitset = 0;

  friend class Index;

  /// first (most significant) bit reserved for creating default space used by
  /// Index that is distinct from the null space
  const static QuantumNumbersAttr reserved;

  /// makes reserved object
  static QuantumNumbersAttr make_reserved() {
    QuantumNumbersAttr result;
    result.bitset = 0x80000000;
    return result;
  }
};  // struct QuantumNumbersAttr

inline const QuantumNumbersAttr QuantumNumbersAttr::null;
inline const QuantumNumbersAttr QuantumNumbersAttr::reserved =
    QuantumNumbersAttr::make_reserved();

/// @brief a collection of attributes which define a space of (1-particle)
/// states
///
/// IndexSpace has a base_label, TypeAttr or interpretable bitset in the
/// context of other spaces. IndexSpace may also have QuantumNumberAttr which
/// differentiate spaces with different quanta, such as spin projection quantum
/// numbers. IndexSpace may additionally have an approximate extent which is
/// useful in symbolic manipulation of indexed expressions, such as
/// tensor network evaluation.
class IndexSpace {
 public:
  /// @brief Attr describes all attributes of a space (occupancy + quantum
  /// numbers)
  struct Attr : TypeAttr, QuantumNumbersAttr {
    Attr(TypeAttr type, QuantumNumbersAttr qns) noexcept
        : TypeAttr(type), QuantumNumbersAttr(qns){};
    Attr(int32_t type, int32_t qns) noexcept
        : TypeAttr(type), QuantumNumbersAttr(qns){};

    /// @brief default ctor creates a null Attr
    /// @post `static_cast<bool>(*this) == false`
    Attr() = default;
    Attr(const Attr &) = default;
    Attr(Attr &&) = default;
    Attr &operator=(const Attr &) = default;
    Attr &operator=(Attr &&) = default;

    const static Attr null;

    constexpr const TypeAttr &type() const {
      return static_cast<const TypeAttr &>(*this);
    }
    constexpr TypeAttr &type() { return static_cast<TypeAttr &>(*this); }
    constexpr const QuantumNumbersAttr &qns() const {
      return static_cast<const QuantumNumbersAttr &>(*this);
    }
    constexpr QuantumNumbersAttr &qns() {
      return static_cast<QuantumNumbersAttr &>(*this);
    }

    constexpr explicit operator int64_t() const {
      return (static_cast<int64_t>(this->type()) << 32) +
             static_cast<int64_t>(this->qns());
    }

    /// @return true if either `type()` or `qns()` is non-null
    explicit operator bool() const {
      return static_cast<bool>(this->type()) || static_cast<bool>(this->qns());
    }

    /// union of Attr = union of TypeAttr and union of QuantumNumbersAttr
    /// @return union of `*this` and @p other, i.e. `*this` AND @p other
    Attr unIon(Attr other) const {
      return {this->type().unIon(other.type()), this->qns().unIon(other.qns())};
    }

    /// @return union of @p a and @p b
    /// @sa Attr::unIon
    friend Attr operator|(Attr a, Attr b) { return a.unIon(b); }

    /// intersection of Attr = intersection of TypeAttr and intersection of
    /// QuantumNumbersAttr
    /// @return intersection of `*this` AND @p other
    Attr intersection(Attr other) const {
      return Attr(this->type().intersection(other.type()),
                  this->qns().intersection(other.qns()));
    }
    /// @return intersection of @p a and @p b
    /// @sa Attr::intersection
    friend Attr operator&(Attr a, Attr b) { return a.intersection(b); }

    /// @return true if \p other is included in this object
    bool includes(Attr other) const {
      return this->type().includes(other.type()) &&
             this->qns().includes(other.qns());
    }

    constexpr bool operator==(Attr other) const {
      return this->type() == other.type() && this->qns() == other.qns();
    }
    constexpr bool operator!=(Attr other) const { return !(*this == other); }

    /// Attr objects are ordered by quantum numbers, then by type
    constexpr bool operator<(Attr other) const {
      if (this->qns() == other.qns()) {
        return this->type() < other.type();
      } else {
        return this->qns() < other.qns();
      }
    }
  };  // struct Attr

  using Type = TypeAttr;
  using QuantumNumbers = QuantumNumbersAttr;

  /// exception type thrown when ancountered unknown/invalid
  /// IndexSpace::base_key() or Index::label()
  struct bad_key : std::invalid_argument {
    bad_key() : std::invalid_argument("bad key") {}
    template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
    bad_key(S &&key)
        : std::invalid_argument(std::string("bad key: ") +
                                sequant::to_string(key)) {}
  };

  struct KeyCompare {
    using is_transparent = void;
    bool operator()(const IndexSpace &a, const IndexSpace &b) const {
      return a.base_key() < b.base_key();
    }
    bool operator()(const std::wstring &a, const IndexSpace &b) const {
      return a < b.base_key();
    }
    bool operator()(const std::wstring_view &a, const IndexSpace &b) const {
      return a < b.base_key();
    }
    bool operator()(const IndexSpace &a, const std::wstring &b) const {
      return a.base_key() < b;
    }
    bool operator()(const IndexSpace &a, const std::wstring_view &b) const {
      return a.base_key() < b;
    }
    bool operator()(const std::wstring &a, const std::wstring &b) const {
      return a < b;
    }
    bool operator()(const std::wstring &a, const std::wstring_view &b) const {
      return a < b;
    }
    bool operator()(const std::wstring_view &a, const std::wstring &b) const {
      return a < b;
    }
  };

  friend constexpr bool operator==(IndexSpace const &,
                                   IndexSpace const &) noexcept;
  friend constexpr bool operator!=(IndexSpace const &,
                                   IndexSpace const &) noexcept;
  friend constexpr bool operator<(IndexSpace const &,
                                  IndexSpace const &) noexcept;

  constexpr Attr attr() const noexcept { return attr_; }
  constexpr Type type() const noexcept { return attr().type(); }
  QuantumNumbers qns() const noexcept { return attr().qns(); }

  /// Default ctor creates null space (with null label, type and quantum
  /// numbers)
  IndexSpace() noexcept {}

  const static IndexSpace null;

  explicit operator bool() const { return *this != null; }

  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpace(S &&type_label, TypeAttr typeattr,
             QuantumNumbersAttr qnattr = QuantumNumbersAttr{0},
             unsigned long approximate_size = 10)
      : attr_(typeattr, qnattr),
        base_key_(sequant::to_wstring(std::forward<S>(type_label))),
        approximate_size_(approximate_size) {}

  IndexSpace(const IndexSpace &other) = default;
  IndexSpace(IndexSpace &&other) = default;
  IndexSpace &operator=(const IndexSpace &other) = default;
  IndexSpace &operator=(IndexSpace &&other) = default;

  const std::wstring &base_key() const { return base_key_; }
  static std::wstring_view reduce_key(std::wstring_view key) {
    const auto underscore_position = key.rfind(L'_');
    if (underscore_position != std::wstring::npos) {  // key can be reduced
      return key.substr(0, underscore_position);
    } else {
      return key;
    }
  }
  static std::wstring reduce_key(std::string_view key) {
    const auto underscore_position = key.rfind(L'_');
    if (underscore_position != std::string::npos) {  // key can be reduced
      return sequant::to_wstring(key.substr(0, underscore_position));
    } else {
      return sequant::to_wstring(key);
    }
  }

  /// @return approximate size of a space
  std::size_t approximate_size() const { return approximate_size_; }

  /// Set the approximate size of a space.
  void approximate_size(size_t n) { approximate_size_ = n; }

 private:
  Attr attr_;
  std::wstring base_key_;
  std::size_t approximate_size_;

  static std::wstring to_wstring(std::wstring_view key) {
    return std::wstring(key.begin(), key.end());
  }
};  // class IndexSpace

inline const IndexSpace IndexSpace::null;
inline const IndexSpace::Attr IndexSpace::Attr::null;

/// @return true if type2 is included in type1, i.e. intersection(type1, type2)
/// == type2
inline bool includes(IndexSpace::Type type1, IndexSpace::Type type2) {
  return type1.includes(type2);
}
/// @return true if qns2 is included in qns1, i.e. \code intersection(qns1,
/// qns2) == qns2 \endcode is true
inline bool includes(IndexSpace::QuantumNumbers qns1,
                     IndexSpace::QuantumNumbers qns2) {
  return qns1.includes(qns2);
}
/// @return true if space2 is included in space1, i.e. intersection(space1,
/// space2) == space2
inline bool includes(const IndexSpace &space1, const IndexSpace &space2) {
  return space1.attr().includes(space2.attr());
}

/// IndexSpace are ordered by their attributes (i.e. labels do not matter one
/// bit)
[[nodiscard]] inline constexpr bool operator<(
    const IndexSpace &space1, const IndexSpace &space2) noexcept {
  return space1.attr() < space2.attr();
}

///
/// IndexSpace are equal if they have equal @c IndexSpace::type(),
/// @c IndexSpace::qns(), and @c IndexSpace::base_key().
///
[[nodiscard]] inline constexpr bool operator==(
    IndexSpace const &space1, IndexSpace const &space2) noexcept {
  return space1.type() == space2.type() && space1.qns() == space2.qns() &&
         space1.base_key() == space2.base_key();
}

///
/// IndexSpace are equal if they have equal @c IndexSpace::type(),
/// @c IndexSpace::qns(), and @c IndexSpace::base_key().
///
[[nodiscard]] inline constexpr bool operator!=(
    IndexSpace const &space1, IndexSpace const &space2) noexcept {
  return !(space1 == space2);
}

}  // namespace sequant

#endif  // SEQUANT_SPACE_H
