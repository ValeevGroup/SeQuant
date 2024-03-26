//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT_SPACE_H
#define SEQUANT_SPACE_H

#include <cassert>
#include <codecvt>
#include <cstdint>
#include <cstdlib>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include <SeQuant/core/container.hpp>

#include <range/v3/algorithm/any_of.hpp>

namespace sequant {

/// @brief TypeAttr denotes the type of index space.
///
/// This class models a host (complete) space partitioned into disjoint
/// subspaces. To simplify implementation of set operations
/// (intersection, union, etc.) it is encoded as a fixed-width (32) bitset.
/// By default, all bits are used, but user can disable some bits from being
/// used for comparison
struct TypeAttr {
  /// @brief used_bits is a bitset of bits that are used for comparison
  static int32_t used_bits;
  int32_t bitset = 0;

  constexpr explicit TypeAttr(int32_t value) noexcept : bitset(value) {}

  [[nodiscard]] constexpr explicit operator int64_t() const {
    return static_cast<int64_t>(bitset);
  }
  [[nodiscard]] constexpr int32_t to_int32() const { return bitset; }
  [[nodiscard]] int32_t to_int32_used() const { return used_bits & bitset; }
  [[nodiscard]] constexpr TypeAttr intersection(TypeAttr other) const {
    return TypeAttr(this->to_int32() & other.to_int32());
  }
  [[nodiscard]] constexpr TypeAttr unIon(TypeAttr other) const {
    return TypeAttr(this->to_int32() | other.to_int32());
  }

  [[nodiscard]] friend bool operator==(TypeAttr lhs, TypeAttr rhs) {
    return lhs.to_int32_used() == rhs.to_int32_used();
  }

  [[nodiscard]] friend bool operator!=(TypeAttr lhs, TypeAttr rhs) {
    return !(lhs == rhs);
  }

  /// @return true if \c other is included in this object
  [[nodiscard]] bool includes(TypeAttr other) const {
    return intersection(other) == other;
  }
  /// @return true if in canonical order this object preceeds \c other
  [[nodiscard]] bool operator<(TypeAttr other) const {
    return this->to_int32_used() < other.to_int32_used();
  }

  /// @return an invalid TypeAttr
  [[nodiscard]] static constexpr TypeAttr invalid() noexcept {
    return TypeAttr(0xffff);
  }
};

/// denotes other quantum numbers (particle type, spin, etc.)
struct QuantumNumbersAttr {
  int32_t bitset = 0;

  constexpr explicit QuantumNumbersAttr(int32_t value) noexcept
      : bitset(value) {}
  constexpr explicit operator int64_t() const {
    return static_cast<int64_t>(bitset);
  }
  constexpr int32_t to_int32() const { return bitset; }
  constexpr QuantumNumbersAttr intersection(QuantumNumbersAttr other) const {
    return QuantumNumbersAttr(this->to_int32() & other.to_int32());
  }
  constexpr QuantumNumbersAttr unIon(QuantumNumbersAttr other) const {
    return QuantumNumbersAttr(this->to_int32() | other.to_int32());
  }
  constexpr QuantumNumbersAttr operator~() const {
    return QuantumNumbersAttr(~this->to_int32());
  }

  friend constexpr bool operator==(QuantumNumbersAttr lhs,
                                   QuantumNumbersAttr rhs) {
    return lhs.to_int32() == rhs.to_int32();
  }
  friend constexpr bool operator!=(QuantumNumbersAttr lhs,
                                   QuantumNumbersAttr rhs) {
    return !(lhs == rhs);
  }

  /// @return true if \c other is included in this object
  constexpr bool includes(QuantumNumbersAttr other) const {
    return intersection(other) == other;
  }
  /// @return true if in canonical order this object preceeds \c other
  constexpr bool operator<(QuantumNumbersAttr other) const {
    return this->to_int32() < other.to_int32();
  }

  /// @return an invalid TypeAttr
  static constexpr QuantumNumbersAttr invalid() noexcept {
    return QuantumNumbersAttr(-0);
  }
};

std::wstring to_wstring(TypeAttr type);
std::wstring to_wstring(QuantumNumbersAttr qns);

/// @brief space of Index objects
///
/// IndexSpace is a set of attributes associated 1-to-1 with keys
class IndexSpace {
 public:
  using TypeAttr = sequant::TypeAttr;
  using QuantumNumbersAttr = sequant::QuantumNumbersAttr;

  /// @brief Attr describes all attributes of a space (occupancy + quantum
  /// numbers)
  struct Attr : TypeAttr, QuantumNumbersAttr {
    Attr(TypeAttr type, QuantumNumbersAttr qns) noexcept
        : TypeAttr(type), QuantumNumbersAttr(qns){};
    Attr(int32_t type, int32_t qns) noexcept
        : TypeAttr(type), QuantumNumbersAttr(qns){};
    //    explicit Attr(int64_t value) : TypeAttr((value & 0xffffffff00000000)
    //    >> 32), QuantumNumbersAttr(value & 0x00000000ffffffff) {}
    Attr(const Attr &) = default;
    Attr(Attr &&) = default;
    Attr &operator=(const Attr &) = default;
    Attr &operator=(Attr &&) = default;

    const TypeAttr &type() const {
      return static_cast<const TypeAttr &>(*this);
    }
    TypeAttr &type() { return static_cast<TypeAttr &>(*this); }
    const QuantumNumbersAttr &qns() const {
      return static_cast<const QuantumNumbersAttr &>(*this);
    }
    QuantumNumbersAttr &qns() {
      return static_cast<QuantumNumbersAttr &>(*this);
    }

    explicit operator int64_t() const {
      return (static_cast<int64_t>(this->type()) << 32) +
             static_cast<int64_t>(this->qns());
    }

    Attr intersection(Attr other) const {
      return Attr(this->type().intersection(other.type()),
                  this->qns().intersection(other.qns()));
    }
    Attr unIon(Attr other) const {
      return Attr(this->type().unIon(other.type()),
                  this->qns().unIon(other.qns()));
    }

    /// @return true if \p other is included in this object
    bool includes(Attr other) const {
      return this->type().includes(other.type()) &&
             this->qns().includes(other.qns());
    }

    friend bool operator==(Attr a, Attr b) {
      return a.type() == b.type() && a.qns() == b.qns();
    }
    friend bool operator!=(Attr a, Attr b) { return !(a == b); }

    static Attr null() noexcept { return Attr{nulltype, nullqns}; }
    static Attr invalid() noexcept {
      return Attr{TypeAttr::invalid(), QuantumNumbersAttr::invalid()};
    }
    bool is_valid() const noexcept { return *this != Attr::invalid(); }

    /// Attr objects are ordered by quantum numbers, then by type
    bool operator<(Attr other) const {
      if (this->qns() == other.qns()) {
        return this->type() < other.type();
      } else {
        return this->qns() < other.qns();
      }
    }
  };
  using Type = TypeAttr;
  using QuantumNumbers = QuantumNumbersAttr;

  /// \name default space tags

  /// @{
  // clang-format off
  /// null space (empty subset), needed to define intersection operation
  static constexpr Type nulltype{0};
  /// represents any space, standard (see below) or otherwise
  static constexpr Type nonnulltype{0x7fffffff};
  /// @}

  /// \name standard space tags

  /// standard space tags are predefined that helps implement set theory of
  /// standard spaces as binary ops on bitsets
  /// @{
  // clang-format off
  /// space of sp states that are fully occupied (i.e., non-correlated) in the reference (vacuum) state and are "frozen" in their reference form
  static constexpr Type frozen_occupied{0b0000001};
  /// space of sp states that are fully occupied (i.e., non-correlated) in the reference (vacuum) state but
  /// can be rotated by mixing with the rest of non-frozen orbitals
  static constexpr Type inactive_occupied{0b0000010};
  /// space of sp states that are fully occupied (i.e., non-correlated) in the reference (vacuum) state but
  /// can be correlated and rotated by mixing with the rest of non-frozen orbitals
  static constexpr Type active_occupied{0b0000100};
  /// space of sp states that are partially occupied (i.e., correlated, or open shells in spin-free single-determinant reference) in the reference (vacuum) state;
  static constexpr Type active{0b0001000};
  /// space of sp states that are fully occupied in the reference (vacuum) state
  /// @note this is the union of IndexSpace::frozen_occupied , IndexSpace::inactive_occupied , IndexSpace::active_occupied
  static constexpr Type occupied = IndexSpace::frozen_occupied.unIon(IndexSpace::inactive_occupied).unIon(IndexSpace::active_occupied);
  /// space of sp states that are fully or partially occupied in the reference (vacuum) state
  /// @note this is the union of IndexSpace::occupied and IndexSpace::active
  static constexpr Type maybe_occupied = IndexSpace::occupied.unIon(IndexSpace::active);
  /// space of sp states that are fully or partially occupied in the reference (vacuum) state and can be correlated
  /// @note this is the union of IndexSpace::active_occupied and IndexSpace::active
  static constexpr Type active_maybe_occupied = IndexSpace::active_occupied.unIon(IndexSpace::active);
  /// space of sp states that are not used to define the reference (vacuum) state (i.e., they are unoccupied) but
  /// can be correlated and rotated by mixing with the rest of non-frozen orbitals
  /// @note unlike IndexSpace::other_unoccupied, these states are supported by a finite computational basis
  static constexpr Type active_unoccupied{0b0010000};
  /// space of sp states that are not used to define the reference (vacuum) state (i.e., they are unoccupied) but
  /// can be rotated by mixing with the rest of non-frozen orbitals
  /// @note unlike IndexSpace::other_unoccupied, these states are supported by a finite computational basis
  static constexpr Type inactive_unoccupied{0b0100000};
  /// space of sp states that are fully unoccupied in the reference (vacuum) state
  /// @note this is the union of IndexSpace::inactive_unoccupied and IndexSpace::active_unoccupied
  /// @note unlike IndexSpace::other_unoccupied, these states are supported by a finite computational basis
  static constexpr Type unoccupied = IndexSpace::active_unoccupied.unIon(IndexSpace::inactive_unoccupied);
  /// space of sp states that are fully or partially unoccupied in the reference (vacuum) state
  /// @note this is the union of IndexSpace::unoccupied  and IndexSpace::active
  /// @note unlike IndexSpace::other_unoccupied, these states are supported by a finite computational basis
  static constexpr Type maybe_unoccupied = IndexSpace::unoccupied.unIon(IndexSpace::active);
  /// space of sp states that are fully or partially unoccupied in the reference (vacuum) state and can be correlated
  /// @note this is the union of IndexSpace::active_unoccupied and IndexSpace::active
  static constexpr Type active_maybe_unoccupied = IndexSpace::active_unoccupied.unIon(IndexSpace::active);
  /// space of sp states that can be correlated
  /// @note this is the union of IndexSpace::active_occupied , IndexSpace::active_unoccupied and IndexSpace::active
  static constexpr Type all_active = IndexSpace::active_occupied.unIon(IndexSpace::active_unoccupied).unIon(IndexSpace::active);
  /// space of sp states represented in computational basis
  /// @note this is the union of IndexSpace::maybe_occupied and IndexSpace::maybe_unoccupied
  static constexpr Type all = IndexSpace::maybe_occupied.unIon(IndexSpace::maybe_unoccupied);
  /// space of sp states that are not used to define the reference (vacuum) state (i.e., they are unoccupied) and not supported
  /// by a supported by a finite computational basis; i.e., these states are the rest of the sp Hilbert space
  static constexpr Type other_unoccupied{0b1000000};
  /// a union of the IndexSpace::other_unoccupied space and IndexSpace::inactive_unoccupied space
  /// @note useful when treating active and inactive unoccupieds differently in e.g. valence-correlated F12 theory
  static constexpr Type complete_inactive_unoccupied = IndexSpace::other_unoccupied.unIon(IndexSpace::inactive_unoccupied);
  /// set of all fully unoccupied states
  /// @note this is a union of IndexSpace::unoccupied and IndexSpace::other_unoccupied
  static constexpr Type complete_unoccupied = IndexSpace::unoccupied.unIon(IndexSpace::other_unoccupied);
  /// set of arbitrary fully or partially unoccupied states
/// @note this is a union of IndexSpace::complete_unoccupied and IndexSpace::active
  static constexpr Type complete_maybe_unoccupied = IndexSpace::complete_unoccupied.unIon(IndexSpace::active);
  /// union of all previous spaces
  /// @note this is a union of IndexSpace::all and IndexSpace::other_unoccupied
  static constexpr Type complete = IndexSpace::all.unIon(IndexSpace::other_unoccupied);
  // clang-format on

  /// list of all standard types
  static constexpr Type standard_types[] = {frozen_occupied,
                                            inactive_occupied,
                                            active_occupied,
                                            occupied,
                                            active,
                                            maybe_occupied,
                                            active_maybe_occupied,
                                            active_unoccupied,
                                            inactive_unoccupied,
                                            unoccupied,
                                            maybe_unoccupied,
                                            active_maybe_unoccupied,
                                            all_active,
                                            all,
                                            other_unoccupied,
                                            complete_inactive_unoccupied,
                                            complete_unoccupied,
                                            complete_maybe_unoccupied,
                                            complete};

  template <int32_t typeint>
  static constexpr bool is_standard_type() {
    return ranges::any_of(standard_types,
                          [](const auto t) { return t == Type{typeint}; });
  }
  /// @}

  /// \name standard quantum numbers tags
  /// \note spin quantum number takes 2 rightmost bits
  /// @{
  /// no quantum numbers
  constexpr static QuantumNumbers nullqns{0b000000};
  /// spin-up
  constexpr static QuantumNumbers alpha{0b000001};
  /// spin-down
  constexpr static QuantumNumbers beta{0b000010};
  /// spin mask
  constexpr static QuantumNumbers spinmask{0b000011};

  /// list of all standard quantum numbers
  static constexpr QuantumNumbers standard_qns[] = {nullqns, alpha, beta};

  template <int32_t qnsint>
  static const constexpr bool is_standard_qns() {
    return ranges::any_of(
        standard_qns, [](const auto t) { return t == QuantumNumbers{qnsint}; });
  }
  /// @}

  struct bad_key : std::invalid_argument {
    bad_key() : std::invalid_argument("bad key") {}
  };
  struct bad_attr : std::invalid_argument {
    bad_attr() : std::invalid_argument("bad attribute") {}
    using std::invalid_argument::invalid_argument;
  };

  struct KeyCompare {
    using is_transparent = void;
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

  /// IndexSpace needs null IndexSpace
  static const IndexSpace &null_instance() { return null_instance_; }
  /// the null IndexSpace is keyed by this key
  static std::wstring null_key() { return L""; }

  /// @brief returns the instance of an IndexSpace object
  /// @param attr the space attribute
  /// @throw bad_attr if \p attr has not been registered
  static const IndexSpace &instance(Attr attr) {
    assert(attr.is_valid());
    if (attr == Attr::null()) return null_instance();
    if (!instance_exists(attr)) throw bad_attr();
    return instances_.find(attr)->second;
  }

  /// @brief returns the instance of an IndexSpace object
  /// @param type the type attribute
  /// @param qns the quantum numbers attribute
  /// @throw bad_key if key not found
  static const IndexSpace &instance(Type type, QuantumNumbers qns = nullqns) {
    const auto attr = Attr(type, qns);
    assert(attr.is_valid());
    if (attr == Attr::null()) return null_instance();
    if (!instance_exists(attr)) {
      std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
      throw bad_attr(converter.to_bytes(L"Request to non-registered space: " +
                                        sequant::to_wstring(type) + L" " +
                                        sequant::to_wstring(qns)));
    }
    return instances_.find(attr)->second;
  }

  /// @brief returns the instance of an IndexSpace object associated
  ///        with the given key
  /// @param key the key associated with this space; this can be either
  ///            the base key used to invoke `IndexSpace::register_instance()`
  ///            or a key used to invoke `IndexSpace::register_key()`
  /// @throw bad_key if key not found
  static const IndexSpace &instance(const std::wstring_view key) {
    if (key == null_key()) return null_instance();
    const auto attr = to_attr(reduce_key(key));
    assert(attr.is_valid());
    if (!instance_exists(attr)) throw bad_key();
    return instances_.find(attr)->second;
  }

  /// @brief constructs a registered instance of an IndexSpace object,
  ///        associates it with a base key
  /// @param base_key string key that will be used as the "base key" for this
  ///        particular space, i.e. the default used for example for
  ///        constructing temporary indices for this space
  /// @param type the IndexSpace::Type attribute to associate \p base_key with
  /// @param qns the IndexSpace::QuantumNumbers attribute to associate \p
  /// base_key with
  /// @param throw_if_already_registered if true, throws an exception if \p
  /// base_key has already been registered
  /// @throw bad_key if \p base_key has already been registered
  /// and \p throw_if_already_registered is true
  static void register_instance(const std::wstring_view base_key, Type type,
                                QuantumNumbers qn = nullqns,
                                bool throw_if_already_registered = true) {
    const auto attr = Attr(type, qn);
    assert(attr.is_valid());
    if (instance_exists(attr) && throw_if_already_registered) throw bad_key();
    const auto irreducible_basekey = reduce_key(base_key);
    const auto irreducible_basekey_str = to_wstring(irreducible_basekey);
    attr2basekey_[attr] = irreducible_basekey_str;
    key2attr_.emplace(irreducible_basekey_str, attr);
    instances_.emplace(std::make_pair(attr, IndexSpace(attr)));
  }

  static bool instance_exists(std::wstring_view key) noexcept {
    return instance_exists(to_attr(reduce_key(key)));
  }

  /// @brief associate a given key with the IndexSpace
  /// @note every IndexSpace constructed via
  ///       `register_instance(base_key,...)` is associated
  ///       with `base_key`; this allows to associated additional
  ///       keys to map to the same IndexSpace
  /// @param key string key that will map to this particular space
  static void register_key(const std::wstring_view key, Type type,
                           QuantumNumbers qn = nullqns,
                           bool throw_if_already_registered = true) {
    const auto attr = Attr(type, qn);
    assert(attr.is_valid());
    const auto irreducible_key = reduce_key(key);
    const auto irreducible_key_str = to_wstring(irreducible_key);
    if (key2attr_.find(irreducible_key_str) != key2attr_.end() &&
        throw_if_already_registered)
      throw bad_key();
    key2attr_.emplace(irreducible_key_str, attr);
  }

  Attr attr() const noexcept {
    assert(attr_.is_valid());
    return attr_;
  }
  Type type() const noexcept { return attr().type(); }
  QuantumNumbers qns() const noexcept { return attr().qns(); }

  /// @brief returns the base key for IndexSpace objects
  /// @param space an IndexSpace object
  /// @throw bad_key if this space has not been registered
  static std::wstring base_key(const IndexSpace &space) {
    return base_key(space.attr());
  }

  /// @brief returns the base key for IndexSpace objects of the given attribute
  /// @param attr the space attribute
  /// @throw bad_attr if \p attr has not been registered
  static std::wstring base_key(Attr attr) {
    assert(attr.is_valid());
    if (attr == Attr::null()) return L"";
    if (!instance_exists(attr)) throw bad_attr();
    return attr2basekey_.find(attr)->second;
  }

  /// Default ctor creates space with nonnull type and null quantum numbers
  IndexSpace() noexcept : attr_(nonnulltype, nullqns) {}

  IndexSpace(const IndexSpace &other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument(
          "IndexSpace copy ctor received invalid argument");
    attr_ = other.attr_;
  }
  IndexSpace(IndexSpace &&other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument(
          "IndexSpace move ctor received invalid argument");
    attr_ = other.attr_;
  }
  IndexSpace &operator=(const IndexSpace &other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument(
          "IndexSpace copy assignment operator received invalid argument");
    attr_ = other.attr_;
    return *this;
  }
  IndexSpace &operator=(IndexSpace &&other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument(
          "IndexSpace move assignment operator received invalid argument");
    attr_ = other.attr_;
    return *this;
  }

  /// @brief constructs an instance of an IndexSpace object
  IndexSpace(Type type, QuantumNumbers qns) noexcept : attr_(type, qns) {
    assert(attr_.is_valid());
  }

 private:
  Attr attr_ = Attr::invalid();
  /// @brief constructs an instance of an IndexSpace object
  explicit IndexSpace(Attr attr) noexcept : attr_(attr) {
    assert(attr.is_valid());
  }

  static container::map<Attr, std::wstring> attr2basekey_;
  static container::map<std::wstring, Attr, KeyCompare> key2attr_;
  static container::map<Attr, IndexSpace> instances_;
  static IndexSpace null_instance_;

  static std::wstring_view reduce_key(std::wstring_view key) {
    const auto underscore_position = key.rfind(L'_');
    return key.substr(0, underscore_position);
  }

  /// @param key the key associated with a registered IndexSpace; this can be
  /// either
  ///            the base key used to invoke `IndexSpace::register_instance()`
  ///            or a key used to invoke `IndexSpace::register_key()`
  /// @return the attribute of the IndexSpace object corresponding to @p key
  /// @throw bad_key if \p key has not been registered
  static Attr to_attr(std::wstring_view key) {
    const auto found_it = key2attr_.find(key);
    if (found_it != key2attr_.end()) return found_it->second;
    throw bad_key();
  }

  static std::wstring to_wstring(std::wstring_view key) {
    return std::wstring(key.begin(), key.end());
  }

  static bool instance_exists(Attr attr) {
    return instances_.find(attr) != instances_.end();
  }
};

inline bool operator==(const IndexSpace &space, IndexSpace::Type t) {
  return space.type() == t;
}
inline bool operator==(IndexSpace::Type t, const IndexSpace &space) {
  return space.type() == t;
}
inline bool operator!=(const IndexSpace &space, IndexSpace::Type t) {
  return !(space == t);
}
inline bool operator!=(IndexSpace::Type t, const IndexSpace &space) {
  return !(t == space);
}
inline bool operator==(const IndexSpace &space,
                       IndexSpace::QuantumNumbers qns) {
  return space.qns() == qns;
}
inline bool operator==(IndexSpace::QuantumNumbers qns,
                       const IndexSpace &space) {
  return space.qns() == qns;
}
inline bool operator!=(const IndexSpace &space,
                       IndexSpace::QuantumNumbers qns) {
  return !(space == qns);
}
inline bool operator!=(IndexSpace::QuantumNumbers qns,
                       const IndexSpace &space) {
  return !(qns == space);
}
inline bool operator==(const IndexSpace &space1, const IndexSpace &space2) {
  return space1.type() == space2.type() && space1.qns() == space2.qns();
}
inline bool operator!=(const IndexSpace &space1, const IndexSpace &space2) {
  return !(space1 == space2);
}
inline IndexSpace::Type intersection(IndexSpace::Type type1,
                                     IndexSpace::Type type2) {
  return type1 == type2 ? type1 : type1.intersection(type2);
}
inline IndexSpace::QuantumNumbers intersection(IndexSpace::QuantumNumbers v1,
                                               IndexSpace::QuantumNumbers v2) {
  return v1 == v2 ? v1 : v1.intersection(v2);
}
inline const IndexSpace &intersection(const IndexSpace &space1,
                                      const IndexSpace &space2) {
  return space1 == space2
             ? space1
             : IndexSpace::instance(space1.attr().intersection(space2.attr()));
}
inline const IndexSpace &intersection(const IndexSpace &space1,
                                      const IndexSpace &space2,
                                      const IndexSpace &space3) {
  return space1 == space2 && space1 == space3
             ? space1
             : IndexSpace::instance(space1.attr().intersection(
                   space2.attr().intersection(space3.attr())));
}
inline IndexSpace::Type unIon(IndexSpace::Type type1, IndexSpace::Type type2) {
  return type1 == type2 ? type1 : type1.unIon(type2);
}
inline IndexSpace::QuantumNumbers unIon(IndexSpace::QuantumNumbers qns1,
                                        IndexSpace::QuantumNumbers qns2) {
  return qns1 == qns2 ? qns1 : qns1.unIon(qns2);
}
inline const IndexSpace &unIon(const IndexSpace &space1,
                               const IndexSpace &space2) {
  return space1 == space2
             ? space1
             : IndexSpace::instance(space1.attr().unIon(space2.attr()));
}
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
inline bool operator<(const IndexSpace &space1, const IndexSpace &space2) {
  return space1.attr() < space2.attr();
}

/// @return -1 if @c space only include orbitals with complete occupancy, +1 if
/// it includes no orbitals with complete occupancy,
///         and 0 of it includes some orbitals with with complete occupancy.
inline int occupancy_class(const IndexSpace &space) {
  const auto included_in_occupied =
      includes(IndexSpace::occupied, space.type());
  const auto included_in_unoccupied =
      includes(IndexSpace::complete_maybe_unoccupied, space.type());
  assert(!(included_in_occupied && included_in_unoccupied));
  if (included_in_occupied && !included_in_unoccupied)
    return -1;
  else if (!included_in_occupied && !included_in_unoccupied)
    return 0;
  else if (!included_in_occupied && included_in_unoccupied)
    return 1;
  abort();  // unreachable
}

std::wstring to_wolfram(const IndexSpace &space);

}  // namespace sequant

#endif  // SEQUANT_SPACE_H
