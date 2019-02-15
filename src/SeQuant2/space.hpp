//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT2_SPACE_H
#define SEQUANT2_SPACE_H

#include <cassert>
#include <bitset>

#include "attr.hpp"
#include "container.hpp"

namespace sequant2 {

/// @brief space of Index objects
///
/// IndexSpace is a set of attributes associated 1-to-1 with keys
class IndexSpace {
 public:

  /// @brief TypeAttr is the type of index space.
  ///
  /// The type is described as a set of (orthogonal) attributes; for simplicity it is encoded as a bitset for ease of
  /// computing.
  struct TypeAttr : std::bitset<32> {
    constexpr explicit TypeAttr(int32_t value) noexcept : std::bitset<32>(static_cast<unsigned long long>(value)) {}
    operator int64_t() const { return static_cast<int64_t>(this->to_ulong()); }
    int32_t to_int32() const { return static_cast<int32_t>(this->to_ulong()); }
    TypeAttr intersection(TypeAttr other) const { return TypeAttr(this->to_int32() & other.to_int32()); }
    TypeAttr unIon(TypeAttr other) const { return TypeAttr(this->to_int32() | other.to_int32()); }
    /// @return true if \c other is included in this object
    bool includes(TypeAttr other) const {
      return intersection(other) == other;
    }
    /// @return true if in canonical order this object preceeds \c other
    bool operator<(TypeAttr other) const {
      return this->to_int32() < other.to_int32();
    }

    /// @return an invalid TypeAttr
    static constexpr TypeAttr invalid() noexcept { return TypeAttr(0xffff); }
  };
  /// denotes other quantum numbers (particle type, spin, etc.)
  struct QuantumNumbersAttr : std::bitset<32> {
    constexpr explicit QuantumNumbersAttr(int32_t value) noexcept : std::bitset<32>(static_cast<unsigned long long>(value)) {}
    explicit operator int64_t() const { return static_cast<int64_t>(this->to_ulong()); }
    int32_t to_int32() const { return static_cast<int32_t>(this->to_ulong()); }
    QuantumNumbersAttr intersection(QuantumNumbersAttr other) const {
      return QuantumNumbersAttr(this->to_int32() & other.to_int32());
    }
    QuantumNumbersAttr unIon(QuantumNumbersAttr other) const {
      return QuantumNumbersAttr(this->to_int32() | other.to_int32());
    }
    /// @return true if \c other is included in this object
    bool includes(QuantumNumbersAttr other) const {
      return intersection(other) == other;
    }
    /// @return true if in canonical order this object preceeds \c other
    bool operator<(QuantumNumbersAttr other) const {
      return this->to_int32() < other.to_int32();
    }

    /// @return an invalid TypeAttr
    static constexpr QuantumNumbersAttr invalid() noexcept { return QuantumNumbersAttr(0xffff); }
  };

  /// @brief Attr describes all attibutes of a space (occupancy + quantum numbers)
  struct Attr : TypeAttr, QuantumNumbersAttr {
    Attr(TypeAttr type, QuantumNumbersAttr qns) noexcept : TypeAttr(type), QuantumNumbersAttr(qns) {};
    Attr(int32_t type, int32_t qns) noexcept : TypeAttr(type), QuantumNumbersAttr(qns) {};
//    explicit Attr(int64_t value) : TypeAttr((value & 0xffffffff00000000) >> 32), QuantumNumbersAttr(value & 0x00000000ffffffff) {}
    Attr(const Attr &) = default;
    Attr(Attr &&) = default;
    Attr &operator=(const Attr &) = default;
    Attr &operator=(Attr &&) = default;

    const TypeAttr &type() const { return static_cast<const TypeAttr &>(*this); }
    TypeAttr &type() { return static_cast<TypeAttr &>(*this); }
    const QuantumNumbersAttr &qns() const { return static_cast<const QuantumNumbersAttr &>(*this); }
    QuantumNumbersAttr &qns() { return static_cast<QuantumNumbersAttr &>(*this); }

    explicit operator int64_t() const {
      return (static_cast<int64_t>(this->type()) << 32) + static_cast<int64_t>(this->qns());
    }

    Attr intersection(Attr other) const {
      return Attr(this->type().intersection(other.type()),
                  this->qns().intersection(other.qns()));
    }
    Attr unIon(Attr other) const { return Attr(this->type().unIon(other.type()), this->qns().unIon(other.qns())); }

    /// @return true if \c other is included in this object
    bool includes(Attr other) const {
      return this->type().includes(other.type()) && this->qns().includes(other.qns());
    }

    bool operator==(Attr other) const { return this->type() == other.type() && this->qns() == other.qns(); }
    bool operator!=(Attr other) const { return !(*this == other); }

    static Attr null() noexcept { return Attr{TypeAttr{0}, QuantumNumbersAttr{0}}; }
    static Attr invalid() noexcept { return Attr{TypeAttr::invalid(), QuantumNumbersAttr::invalid()}; }
    bool is_valid() const noexcept { return *this != Attr::invalid(); }

    /// Attr objects are ordered by quantum numbers, then by type
    bool operator<(Attr other) const {
      if (this->qns() < other.qns()) {
        return true;
      } else {
        return !(this->qns() == other.qns()) || this->type() < other.type();
      }
    }
  };
  using Type = TypeAttr;
  using QuantumNumbers = QuantumNumbersAttr;

  /// standard space tags are predefined that helps implement set theory of standard spaces as binary ops on bitsets
  static Type frozen_occupied;
  static Type inactive_occupied;
  static Type active_occupied;
  static Type occupied;
  static Type active_unoccupied;
  static Type inactive_unoccupied;
  static Type unoccupied;
  static Type all_active;
  static Type all;
  static Type other_unoccupied;
  static Type complete_unoccupied;
  static Type complete;
  template <int32_t typeint> static const constexpr bool is_standard_type() {
    const Type type{typeint};
    return (type == frozen_occupied || type == inactive_occupied || type == active_occupied ||
        type == occupied || type == active_unoccupied || type == inactive_unoccupied ||
        type == unoccupied || type == all_active || type == all || type == other_unoccupied ||
        type == complete_unoccupied || type == complete);
  }

  /// standard space tags are predefined that helps implement set theory of standard spaces as binary ops on bitsets
  static QuantumNumbers nullqns;  //!< no quantum numbers
  static QuantumNumbers alpha;  //!< spin-up
  static QuantumNumbers beta;  //!< spin-down
  template <int32_t qnsint> static const constexpr bool is_standard_qns() {
    const QuantumNumbers qns{qnsint};
    return (qns == nullqns || qns == alpha || qns == beta);
  }

  struct bad_key : std::invalid_argument {
    bad_key() : std::invalid_argument("bad key") {}
  };
  struct bad_attr : std::invalid_argument {
    bad_attr() : std::invalid_argument("bad attribute") {}
  };

  /// IndexSpace needs null IndexSpace
  static const IndexSpace &null_instance() {
    return null_instance_;
  }
  /// the null IndexSpace is keyed by this key
  static std::wstring null_key() {
    return L"";
  }

  /// @brief returns the instance of an IndexSpace object
  /// @param attr the space attribute
  /// @throw bad_key if key not found
  static const IndexSpace &instance(Attr attr) {
    assert(attr.is_valid());
    if (attr == Attr::null())
      return null_instance();
    if (!instance_exists(attr))
      throw bad_attr();
    return instances_.find(attr)->second;
  }

  /// @brief returns the instance of an IndexSpace object
  /// @param type the type attribute
  /// @param qns the quantum numbers attribute
  /// @throw bad_key if key not found
  static const IndexSpace &instance(Type type, QuantumNumbers qns = nullqns) {
    const auto attr = Attr(type, qns);
    assert(attr.is_valid());
    if (attr == Attr::null())
      return null_instance();
    if (!instance_exists(attr))
      throw bad_attr();
    return instances_.find(attr)->second;
  }

  /// @brief returns the instance of an IndexSpace object
  /// @param key a string key describing a particular space that has been registered before
  /// @throw bad_key if key not found
  static const IndexSpace &instance(const std::wstring_view key) {
    if (key == null_key())
      return null_instance();
    const auto attr = to_attr(reduce_key(key));
    assert(attr.is_valid());
    if (!instance_exists(attr))
      throw bad_key();
    return instances_.find(attr)->second;
  }

  /// @brief returns the instance of an IndexSpace object
  /// @param key string key describing a particular space
  static void register_instance(const std::wstring_view key,
                                Type type,
                                QuantumNumbers qn = nullqns,
                                bool throw_if_already_registered = true) {
    const auto attr = Attr(type, qn);
    assert(attr.is_valid());
    if (instance_exists(attr) && throw_if_already_registered)
      throw bad_key();
    const auto irreducible_key = reduce_key(key);
    keys_[attr] = to_wstring(irreducible_key);
    instances_.emplace(std::make_pair(attr, IndexSpace(attr)));
  }

  /// @brief returns the instance of an IndexSpace object
  /// @param key string key describing a particular space
  static void register_standard_instances() {
    const bool do_not_throw = false;
    register_instance(L"i", active_occupied, nullqns, do_not_throw);
    register_instance(L"m", occupied, nullqns, do_not_throw);
    register_instance(L"a", active_unoccupied, nullqns, do_not_throw);
    register_instance(L"e", unoccupied, nullqns, do_not_throw);
    register_instance(L"p", all, nullqns, do_not_throw);
    register_instance(L"⍺'", other_unoccupied, nullqns, do_not_throw);
    register_instance(L"⍺", complete_unoccupied, nullqns, do_not_throw);
    register_instance(L"κ", complete, nullqns, do_not_throw);
  }

  static bool instance_exists(std::wstring_view key) noexcept {
    return instance_exists(to_attr(reduce_key(key)));
  }

  Attr attr() const noexcept {
    assert(attr_.is_valid());
    return attr_;
  }
  Type type() const noexcept {
    return attr().type();
  }
  QuantumNumbers qns() const noexcept {
    return attr().qns();
  }

  /// @brief returns the base key for IndexSpace objects
  /// @param space an IndexSpace object
  /// @throw bad_key if this space has not beed registered
  static std::wstring base_key(const IndexSpace& space) {
    return base_key(space.attr());
  }

  /// @brief returns the base key for IndexSpace objects of the given attribute
  /// @param attr the space attribute
  /// @throw bad_key if this object has not beed registered
  static std::wstring base_key(Attr attr) {
    assert(attr.is_valid());
    if (attr == Attr::null())
      return L"";
    if (!instance_exists(attr))
      throw bad_attr();
    return keys_.find(attr)->second;
  }

  /// Default ctor creates an invalid space
  IndexSpace() : attr_(Attr::invalid()) {}

  IndexSpace(const IndexSpace &other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument("IndexSpace copy ctor received invalid argument");
    attr_ = other.attr_;
  }
  IndexSpace(IndexSpace &&other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument("IndexSpace move ctor received invalid argument");
    attr_ = other.attr_;
  }
  IndexSpace &operator=(const IndexSpace &other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument("IndexSpace copy assignment operator received invalid argument");
    attr_ = other.attr_;
    return *this;
  }
  IndexSpace &operator=(IndexSpace &&other) {
    if (!other.attr().is_valid())
      throw std::invalid_argument("IndexSpace move assignment operator received invalid argument");
    attr_ = other.attr_;
    return *this;
  }

 private:

  Attr attr_ = Attr::invalid();
  /// @brief constructs an instance of an IndexSpace object
  explicit IndexSpace(Attr attr) noexcept : attr_(attr) {
    assert(attr.is_valid());
  }
  /// @brief constructs an instance of an IndexSpace object
  IndexSpace(Type type, QuantumNumbers qns) noexcept : attr_(type, qns) {
    assert(attr_.is_valid());
  }

  static container::map<Attr, std::wstring> keys_;
  static container::map<Attr, IndexSpace> instances_;
  static IndexSpace null_instance_;

  static std::wstring_view reduce_key(std::wstring_view key) {
    const auto underscore_position = key.find(L'_');
    return key.substr(0, underscore_position);
  }

  static Attr to_attr(std::wstring_view key) {
    for (const auto &attr_key: keys_) {
      if (attr_key.second == key)
        return attr_key.first;
    }
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
inline bool operator==(const IndexSpace &space, IndexSpace::QuantumNumbers qns) {
  return space.qns() == qns;
}
inline bool operator==(IndexSpace::QuantumNumbers qns, const IndexSpace &space) {
  return space.qns() == qns;
}
inline bool operator!=(const IndexSpace &space, IndexSpace::QuantumNumbers qns) {
  return !(space == qns);
}
inline bool operator!=(IndexSpace::QuantumNumbers qns, const IndexSpace &space) {
  return !(qns == space);
}
inline bool operator==(const IndexSpace &space1, const IndexSpace &space2) {
  return space1.type() == space2.type() && space1.qns() == space2.qns();
}
inline bool operator!=(const IndexSpace &space1, const IndexSpace &space2) {
  return !(space1 == space2);
}
inline IndexSpace::Type intersection(IndexSpace::Type type1, IndexSpace::Type type2) {
  return type1.intersection(type2);
}
inline IndexSpace::QuantumNumbers intersection(IndexSpace::QuantumNumbers v1, IndexSpace::QuantumNumbers v2) {
  return v1.intersection(v2);
}
inline const IndexSpace &intersection(const IndexSpace &space1, const IndexSpace &space2) {
  return IndexSpace::instance(space1.attr().intersection(space2.attr()));
}
inline const IndexSpace &intersection(const IndexSpace &space1, const IndexSpace &space2, const IndexSpace &space3) {
  return IndexSpace::instance(space1.attr().intersection(space2.attr().intersection(space3.attr())));
}
inline IndexSpace::Type unIon(IndexSpace::Type type1, IndexSpace::Type type2) {
  return type1.unIon(type2);
}
inline IndexSpace::QuantumNumbers unIon(IndexSpace::QuantumNumbers qns1, IndexSpace::QuantumNumbers qns2) {
  return qns1.unIon(qns2);
}
inline const IndexSpace &unIon(const IndexSpace &space1, const IndexSpace &space2) {
  return IndexSpace::instance(space1.attr().unIon(space2.attr()));
}
/// @return true if type2 is included in type1, i.e. intersection(type1, type2) == type2
inline bool includes(IndexSpace::Type type1, IndexSpace::Type type2) {
  return type1.includes(type2);
}
/// @return true if qns1 is included in qns2, i.e. \code intersection(qns1, qns2) == qns2 \endcode is true
inline bool includes(IndexSpace::QuantumNumbers qns1, IndexSpace::QuantumNumbers qns2) {
  return qns1.includes(qns2);
}
/// @return true if space1 is included in space2, i.e. intersection(space1, space2) == space2
inline bool includes(const IndexSpace &space1, const IndexSpace &space2) {
  return space1.attr().includes(space2.attr());
}

/// IndexSpace are ordered by their attributes (i.e. labels do not matter one bit)
inline bool operator<(const IndexSpace &space1, const IndexSpace &space2) {
  return space1.attr() < space2.attr();
}


/// @return -1 if @c space includes no orbitals with zero occupancy, +1 if it includes only orbitals with zero occupancy,
///         and 0 of it includes some orbitals with zero occupancy.
inline int occupancy_class(const IndexSpace& space) {
  const auto included_in_occupied = includes(IndexSpace::occupied, space.type());
  const auto included_in_unoccupied = includes(IndexSpace::complete_unoccupied, space.type());
  assert(!(included_in_occupied && included_in_unoccupied));
  if (included_in_occupied && !included_in_unoccupied)
    return -1;
  else if (!included_in_occupied && !included_in_unoccupied)
    return 0;
  else if (!included_in_occupied && included_in_unoccupied)
    return 1;
  abort(); // unreachable
}

std::wstring to_wolfram(const IndexSpace &space);

}  // namespace sequant2

#endif //SEQUANT2_SPACE_H
