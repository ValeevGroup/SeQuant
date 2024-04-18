//
// Created by Eduard Valeyev on 3/20/18.
//

#ifndef SEQUANT_SPACE_H
#define SEQUANT_SPACE_H

#include <bitset>
#include <cassert>
#include <cmath>

#include "attr.hpp"
#include "container.hpp"

#include <range/v3/algorithm/any_of.hpp>

namespace sequant {

/// @brief TypeAttr denotes the type of index space.
///
/// This class models a host (complete) space partitioned into disjoint
/// subspaces. To simplify implementation of set operations
/// (intersection, union, etc.) it is encoded as a fixed-width (32) bitset.
struct TypeAttr {
  int32_t bitset = 0;

  TypeAttr(int32_t value) : bitset(value) {}

  constexpr explicit operator int64_t() const {
    return static_cast<int64_t>(bitset);
  }
  const int32_t to_int32() const { return bitset; }
  const TypeAttr intersection(TypeAttr other) const {
    return TypeAttr(this->to_int32() & other.to_int32());
  }
  const TypeAttr unIon(TypeAttr other) const {
    return TypeAttr(this->to_int32() | other.to_int32());
  }
  const TypeAttr exclusionary_or(TypeAttr other) const {
    return TypeAttr(this->to_int32() xor other.to_int32());
  }

  bool const operator==(const TypeAttr other) const {
    return this->to_int32() == other.to_int32();
  }
  bool const operator!=(TypeAttr other) const { return !(*this == other); }

  /// @return true if \c other is included in this object
  const bool includes(TypeAttr other) const {
    return intersection(other) == other;
  }
  /// @return true if in canonical order this object preceeds \c other
  bool const operator<(TypeAttr other) const {
    return this->to_int32() < other.to_int32();
  }

  /// @return an invalid TypeAttr
  static TypeAttr invalid() noexcept { return TypeAttr(0xffff); }
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

  friend constexpr bool operator==(QuantumNumbersAttr, QuantumNumbersAttr);
  friend constexpr bool operator!=(QuantumNumbersAttr, QuantumNumbersAttr);

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

constexpr bool operator==(QuantumNumbersAttr lhs, QuantumNumbersAttr rhs) {
  return lhs.to_int32() == rhs.to_int32();
}

constexpr bool operator!=(QuantumNumbersAttr lhs, QuantumNumbersAttr rhs) {
  return !(lhs == rhs);
}

/// @brief a collection of attributes which define a space or set
/// IndexSpaces have a base_label, TypeAttr or interpretable bitset in the
/// context of other spaces. IndexSpaces may also have QuantumNumberAttr which
/// are differentiate spaces with different quanta such as fermionic spin
/// IndesSpaces may additionally have an approximate size which is useful in
/// expression optimization in the context of multiple spaces.
class IndexSpace {
 public:
  using TypeAttr = sequant::TypeAttr;
  using QuantumNumbersAttr = sequant::QuantumNumbersAttr;

  // IndexSpace(IndexSpace& idxspace);

  IndexSpace(std::wstring_view base_label, TypeAttr typeattr_,
             QuantumNumbersAttr qnattr_ = QuantumNumbersAttr{0},
             unsigned long approximate_size = 10) {
    attr_ = Attr(typeattr_, qnattr_);
    base_key_ = base_label;
    approximate_size_ = approximate_size;
  }

  IndexSpace(std::wstring_view base_label, int32_t space_bitset,
             int32_t qn_bitset = 0) {
    attr_ = Attr(space_bitset, qn_bitset);
    base_key_ = base_label;
  }

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

    std::vector<Attr> irreducible_reps() const {
      std::vector<Attr> result;
      std::bitset<32> bit32(this->sequant::TypeAttr::to_int32());
      for (int i = 0; i < bit32.size(); i++) {
        if (bit32[i]) {
          Attr temp(static_cast<int>(std::pow(2, i)), this->qns().to_int32());
          result.push_back(temp);
        }
      }
      return result;
    }

    // make a list of all excluded spaces between two index spaces that at least
    // one of them contains.

    std::vector<Attr> excluded_spaces(Attr other) const {
      std::vector<Attr> result;

      // if the excluded space simply forms a union of the two spaces, return
      // the two spaces unchanged order smallest to largest
      // TODO try and break this in unit tests
      if (this->type().unIon(other.type()).to_int32() ==
          this->exclusionary_or(other).to_int32()) {
        if (this->type().to_int32() < other.type().to_int32()) {
          result.push_back(*this);
          result.push_back(other);
        } else {
          result.push_back(other);
          result.push_back(*this);
        }
        return result;
      }
      std::bitset<32> xor_bitset(this->exclusionary_or(other).to_int32());
      std::vector<std::pair<int, int>> start_stop_ranges;
      /// TODO need to make a cleaner implementation here.
      // std::bitset does not have an iterator
      int temp_start = 33;
      int temp_end = -1;
      for (int i = 0; i < 32; i++) {
        if (xor_bitset[i]) {
          if (i > temp_end && i < temp_start) {
            temp_start = i;
          }
          temp_end = i;
        } else {
          if (temp_end >= 0 && temp_start != 33) {
            start_stop_ranges.push_back({temp_start, temp_end});
            temp_start = 33;
          }
        }
      }
      for (int i = 0; i < start_stop_ranges.size(); i++) {
        std::bitset<32> new_bitspace;
        for (int j = start_stop_ranges[i].first;
             j <= start_stop_ranges[i].second; j++) {
          new_bitspace.set(j, true);
        }
        Attr new_attr(new_bitspace.to_ulong(), this->qns().to_int32());
        result.push_back(new_attr);
      }
      if (result.empty()) {  // allow null result
        result.push_back({0, 0});
      }
      return result;
    }

    /// @return true if \p other is included in this object
    bool includes(Attr other) const {
      return this->type().includes(other.type()) &&
             this->qns().includes(other.qns());
    }

    bool operator==(Attr other) const {
      return this->type() == other.type() && this->qns() == other.qns();
    }
    bool operator!=(Attr other) const { return !(*this == other); }

    static Attr null() noexcept { return Attr{0, nullqns}; }
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
  /// @}

  /// @{
  // clang-format off



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
  bool operator==(IndexSpace IS)const{return this->type() == IS.type() && this->get_base_key() == IS.get_base_key() && this->qns() == IS.qns() ? true : false;}

  bool operator!=(IndexSpace IS)const{return !(*this == IS);}


  Attr attr() const noexcept {
    assert(attr_.is_valid());
    return attr_;
  }
  Type type() const noexcept { return attr().type(); }
  QuantumNumbers qns() const noexcept { return attr().qns(); }



  /// Default ctor creates space with nonnull type and null quantum numbers
  IndexSpace() {
    attr_ = {{0x7fffffff}, nullqns};
    base_key_ =L"";
    approximate_size_ = 0;
  }

  IndexSpace(const IndexSpace &other) {
  attr_ = other.get_attr();
  base_key_ = other.get_base_key();
  approximate_size_ = other.get_approximate_size();
  }
  IndexSpace(IndexSpace &&other) {
     attr_ = other.get_attr();
     base_key_ = other.get_base_key();
     approximate_size_ = other.get_approximate_size();
  }
  IndexSpace &operator=(const IndexSpace &other) {
    attr_ = other.get_attr();
    base_key_ = other.get_base_key();
    approximate_size_ = other.get_approximate_size();
    return *this;
  }
  IndexSpace &operator=(IndexSpace &&other) {
    attr_ = other.get_attr();
    base_key_ = other.get_base_key();
    approximate_size_ = other.get_approximate_size();
    return *this;
  }


  const Attr get_attr()const {
    return attr_;
  }
  std::wstring get_base_key() const{
    return base_key_;
  }
  static std::wstring_view reduce_key(std::wstring_view key) {
    const auto underscore_position = key.rfind(L'_');
    if(underscore_position != -1){ // key can be reduced
      return key.substr(0, underscore_position);
    }
    else{
      return key;
    }
  }

  // having approximate sizes of spaces can be helpful when optimizing evaluation order.
  unsigned long get_approximate_size() const {
    return approximate_size_;
  }

 private:
  Attr attr_ = Attr::invalid();
  std::wstring base_key_;
  std::size_t approximate_size_;
  /// @brief constructs an instance of an IndexSpace object
  //explicit IndexSpace(Attr attr) noexcept : attr_(attr) {
  //  assert(attr.is_valid());
  //}



  static std::wstring to_wstring(std::wstring_view key) {
    return std::wstring(key.begin(), key.end());
  }

};

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


/*std::wstring to_wolfram(const IndexSpace space){
 throw std::logic_error("not implemented");
};*/

}  // namespace sequant

#endif  // SEQUANT_SPACE_H
