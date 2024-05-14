//
// Created by Conner Masteran on 4/16/24.
//

#ifndef SEQUANT_INDEX_SPACE_REGISTRY_HPP
#define SEQUANT_INDEX_SPACE_REGISTRY_HPP

#include "space.hpp"

#include <boost/hana.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>

namespace sequant {

inline namespace space_tags {
struct IsVacuumOccupied {};
struct IsReferenceOccupied {};
struct IsComplete {};
struct IsHole {};
struct IsParticle {};

constexpr auto is_vacuum_occupied = IsVacuumOccupied{};
constexpr auto is_reference_occupied = IsReferenceOccupied{};
constexpr auto is_complete = IsComplete{};
constexpr auto is_hole = IsHole{};
constexpr auto is_particle = IsParticle{};

}  // namespace space_tags

/// @brief set of known IndexSpace objects

/// Each IndexSpace object has hardwired base key (label) that gives
/// indexed expressions appropriate semantics; e.g., spaces referred to by
/// indices in \f$ t_{p_1}^{i_1} \f$ are defined if IndexSpace objects with
/// base keys \f$ p \f$ and \f$ i \f$ are registered.
/// Since index spaces have set-theoretic semantics, the user must
/// provide complete set of unions/intersects of the base spaces to
/// cover all possible IndexSpace objects that can be generated in their
/// program.
/// @note  @c IndexSpaces with @c typeattr corresponding to occupied indices
/// will always be lower than typeattr corresponding to unoccupied orbitals.
class IndexSpaceRegistry {
 public:
  IndexSpaceRegistry() {
    // register nullspace
    this->add(nullspace);
  }

  IndexSpaceRegistry(const IndexSpaceRegistry& other) = default;
  IndexSpaceRegistry(IndexSpaceRegistry&& other) = default;
  IndexSpaceRegistry& operator=(const IndexSpaceRegistry& other) = default;
  IndexSpaceRegistry& operator=(IndexSpaceRegistry&& other) = default;

  decltype(auto) begin() const { return spaces.begin(); }
  decltype(auto) end() const { return spaces.end(); }

  /// @brief retrieve an IndexSpace from the registry by the label
  /// @param label can be numbered or the @c base_key
  /// @return IndexSpace associated with that key
  /// @throw IndexSpace::bad_key if matching space is not found
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  const IndexSpace& retrieve(S&& label) const {
    auto it = spaces.find(IndexSpace::reduce_key(to_basic_string_view(label)));
    if (it == spaces.end()) {
      throw IndexSpace::bad_key(label);
    }
    return *it;
  }

  /// @brief retrieve an IndexSpace from the registry by its type and quantum
  /// numbers
  /// @param type IndexSpace::Type
  /// @param qns IndexSpace::QuantumNumbers
  /// @return IndexSpace associated with that key.
  /// @throw std::invalid_argument if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Type& type,
                             const IndexSpace::QuantumNumbers& qns) const {
    auto it = std::find_if(spaces.begin(), spaces.end(), [&](const auto& is) {
      return is.type() == type && is.qns() == qns;
    });
    if (it == spaces.end()) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::retrieve(type,qn): missing { IndexSpace::Type=" +
          std::to_string(type.to_int32()) + " , IndexSpace::QuantumNumbers=" +
          std::to_string(qns.to_int32()) + " } combination");
    }
    return *it;
  }

  /// @name adding IndexSpace objects to the registry
  /// @{

  /// @brief add an IndexSpace to this registry.
  /// @param IS an IndexSpace
  /// @return reference to `this`
  /// @throw std::invalid_argument if trying to add an IndexSpace with that
  /// existing base key or attr.
  IndexSpaceRegistry& add(const IndexSpace& IS) {
    auto it = spaces.find(IS.base_key());
    if (it != spaces.end()) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add(is): already have an IndexSpace associated "
          "with is.base_key(); if you are trying to replace the IndexSpace use "
          "IndexSpaceRegistry::replace(is)");
    } else {
#ifndef NDEBUG
      // make sure there are no duplicate IndexSpaces whose attribute is
      // IS.attr()
      if (std::find_if(spaces.begin(), spaces.end(), [&IS](auto&& is) {
            return IS.attr() == is.attr();
          }) != spaces.end()) {
        throw std::invalid_argument(
            "IndexSpaceRegistry::add(is): already have an IndexSpace "
            "associated with is.attr(); if you are trying to replace the "
            "IndexSpace use IndexSpaceRegistry::replace(is)");
      }
#endif
      spaces.emplace(IS);
    }
    return *this;
  }

  template <typename S, typename... OptionalArgs,
            typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& add(S&& type_label, IndexSpace::Type type,
                          OptionalArgs&&... args) {
    auto h_args = boost::hana::make_tuple(args...);
    auto h_types = boost::hana::to<boost::hana::tuple_tag>(
        boost::hana::tuple_t<std::remove_reference_t<OptionalArgs>...>);

    // process IndexSpace::QuantumNumbers, set to default is not given
    auto h_qns = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::type_c<decltype(arg)> ==
             boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    constexpr auto nqns = boost::hana::size(h_qns);
    static_assert(
        nqns == boost::hana::size_c<0> || nqns == boost::hana::size_c<1>,
        "IndexSpaceRegistry::add: only one IndexSpace::QuantumNumbers argument "
        "is allowed");
    constexpr auto have_qns = nqns == boost::hana::size_c<1>;
    IndexSpace::QuantumNumbersAttr qns;
    if constexpr (have_qns) {
      qns = boost::hana::at_c<0>(h_qns);
    }

    // process approximate_size, set to default is not given
    auto h_ints = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::traits::is_integral(boost::hana::decltype_(arg));
    });
    constexpr auto nints = boost::hana::size(h_ints);
    static_assert(
        nints == boost::hana::size_c<0> || nints == boost::hana::size_c<1>,
        "IndexSpaceRegistry::add: only one integral argument is allowed");
    constexpr auto have_approximate_size = nints == boost::hana::size_c<1>;
    unsigned long approximate_size = 10;
    if constexpr (have_approximate_size) {
      approximate_size = boost::hana::at_c<0>(h_ints);
    }

    // make space
    IndexSpace space(std::forward<S>(type_label), type, qns, approximate_size);
    this->add(space);

    // process non-integer tags
    auto h_nonints = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
          boost::hana::type_c<decltype(arg)>);
    });
    boost::hana::for_each(h_nonints, [this, &type](auto arg) {
      if constexpr (boost::hana::type_c<decltype(arg)> ==
                    boost::hana::type_c<space_tags::IsVacuumOccupied>) {
        this->vacuum_occupied_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<
                               space_tags::IsReferenceOccupied>) {
        this->reference_occupied_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsComplete>) {
        this->complete_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsHole>) {
        this->active_hole_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsParticle>) {
        this->active_particle_space(type);
      } else {
        static_assert(meta::always_false<decltype(arg)>::value,
                      "IndexSpaceRegistry::add: unknown tag");
      }
    });

    return *this;
  }

  /// @}

  /// @brief removes an IndexSpace associated with `IS.base_key()` from this
  /// @param IS an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& remove(const IndexSpace& IS) {
    auto it = spaces.find(IS.base_key());
    if (it != spaces.end()) {
      spaces.erase(IS);
    }
    return *this;
  }

  /// @brief equivalent to `remove(this->retrieve(label))`
  /// @param label space label
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& remove(S&& label) {
    auto&& IS = this->retrieve(std::forward<S>(label));
    return this->remove(IS);
  }

  /// @brief replaces an IndexSpace registered in the registry under
  /// IS.base_key()
  ///        with @p IS
  /// @param IS an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& replace(const IndexSpace& IS) {
    this->remove(IS);
    return this->add(IS);
  }

  /// @brief returns the list of base IndexSpace objects in the order of
  /// increasing @c typeattr

  /// An IndexSpace object is base if it has a single bit set in its @c
  /// typeattr.
  std::vector<IndexSpace> base_spaces() const {
    std::vector<IndexSpace> result;
    // if the registered instance has a single bit true, it is a base space.
    //  std::map should order this new map via IndexSpace.
    for (auto&& space : spaces) {
      if (has_single_bit(space.type().to_int32())) {
        result.push_back(space);
      }
    }
    std::sort(result.begin(), result.end(),
              [](auto&& a, auto&& b) -> bool { return a.type() < b.type(); });
    return result;
  }

  /// @brief does the space have only one bit set
  /// @param IS IndexSpace
  /// @note list of all base spaces is computed each time. Otherwise, a separate
  /// list would need to be maintained for every mutation of the registry. This
  /// seems safer.
  static bool is_base_space(const IndexSpace& IS) {
    return has_single_bit(IS.type().to_int32());
  }

  /// @brief clear the IndexSpaceRegistry map
  /// @return reference to `this`
  IndexSpaceRegistry& clear_registry() {
    spaces.clear();
    return this->add(nullspace);
  }

  ///@brief Is the result of a binary operation null or not registered return
  /// false.
  /// a user may wish to know if an operation returns a space they have
  /// registered.
  /// @param i1 IndexSpace
  /// @param i2 IndexSpace
  /// @param op a function which takes two int32 as arguments and returns an
  /// int32
  /// @note IndexSpaces must have the same @c QuantumNumberAttr to be a valid
  /// bitop
  bool valid_bitop(const IndexSpace& i1, const IndexSpace& i2,
                   const std::function<int32_t(int32_t, int32_t)>& op) {
    auto bitop_int = op(i1.type().to_int32(), i2.type().to_int32());
    bool same_qn = i1.qns() == i2.qns();
    if (!same_qn) return false;
    auto& temp_space = find_by_attr({bitop_int, i1.qns()});
    return temp_space == nullspace ? false : true;
  }

  ///@brief return the resulting space corresponding to a bitwise intersection
  /// between two spaces.
  /// @param space1
  /// @param space2
  /// @return the resulting space after intesection
  /// @note can return nullspace
  /// @note throw invalid_argument if the bitwise result is not registered
  const IndexSpace& intersection(const IndexSpace& space1,
                                 const IndexSpace& space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {  // spaces with different quantum numbers do not
                        // intersect.
        return nullspace;
      }
      auto intersection_attr = space1.type().intersection(space2.type());
      const IndexSpace& intersection_space =
          find_by_attr({intersection_attr, space1.qns()});
      // the nullspace is a reasonable return value for intersection
      if (intersection_space == nullspace && intersection_attr != 0) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return intersection_space;
      }
    }
  }

  ///@brief return the resulting space corresponding to a bitwise intersection
  /// between two spaces.
  /// @param space1
  /// @param space2
  /// @param space3
  /// @return the resulting space after intesection
  /// @note can return nullspace
  /// @note throw invalid_argument if the bitwise result is not registered
  const IndexSpace& intersection(const IndexSpace& space1,
                                 const IndexSpace& space2,
                                 const IndexSpace& space3) const {
    if (space1 == space2 && space1 == space3) {
      return space1;
    } else {
      bool same_qns =
          ((space1.qns() == space2.qns()) && (space1.qns() == space3.qns()));
      if (!same_qns) {  // spaces with different quantum numbers do not
                        // intersect.
        return nullspace;
      }
      auto intersection_attr =
          space1.type().intersection(space2.type()).intersection(space3.type());
      const IndexSpace& intersection_space =
          find_by_attr({intersection_attr, space1.qns()});
      if (intersection_space == nullspace && intersection_attr != 0) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return intersection_space;
      }
    }
  }

  /// @param space1
  /// @param space2
  /// @return the union of two spaces.
  /// @note can only return registered spaces
  /// @note never returns nullspace
  const IndexSpace& unIon(const IndexSpace& space1,
                          const IndexSpace& space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {
        throw std::invalid_argument(
            "asking for the intersection of spaces with incompatible quantum "
            "number attributes.");
      }
      auto unIontype = space1.type().unIon(space2.type());
      const IndexSpace& unIonSpace = find_by_attr({unIontype, space1.qns()});
      if (unIonSpace == nullspace) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return unIonSpace;
      }
    }
  }
  /// @brief which spaces result from XOR bitwise operation of two spaces. Only
  /// keep connected spaces.
  /// @param space1
  /// @param space2
  /// @return a list of spaces
  /// @note nullspace is not a valid return
  /// @note finding unregistered @c typeattr will throw
  std::vector<IndexSpace> non_overlapping_spaces(
      const IndexSpace& space1, const IndexSpace& space2) const {
    auto attributes = space1.attr().excluded_spaces(space2.attr());
    std::vector<IndexSpace> result;
    for (int i = 0; i < attributes.size(); i++) {
      auto excluded_space = find_by_attr(attributes[i]);
      if (excluded_space == nullspace) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      }
      result.push_back(excluded_space);
    }
    return result;
  }

  ///@brief do two spaces have non-overlaping bitsets.
  /// @note does not probe the registry for these spaces
  bool has_non_overlapping_spaces(const IndexSpace& space1,
                                  const IndexSpace& space2) const {
    return space1.type().exclusionary_or(space2.type()).to_int32() != 0;
  }

  ///@brief an @c IndexSpace is occupied with respect to the fermi vacuum or a
  /// subset of that space
  /// @note only makes sense to ask this if in a SingleProduct vacuum context.
  bool is_pure_occupied(const IndexSpace& IS) const {
    if (IS == nullspace) {
      return false;
    }
    if (IS.type().to_int32() <=
        vacuum_occupied_space(IS.qns()).type().to_int32()) {
      return true;
    } else {
      return false;
    }
  }

  ///@brief all states are unoccupied in the fermi vacuum
  ///@note again, this only makes sense to ask if in a SingleProduct vacuum
  /// context.
  bool is_pure_unoccupied(const IndexSpace& IS) const {
    if (IS == nullspace) {
      return false;
    } else {
      return IS.type().intersection(vacuum_occupied_space(IS.qns()).type()) ==
             nulltype;
    }
  }

  ///@brief some states are fermi vacuum occupied
  bool contains_occupied(const IndexSpace& IS) const {
    return IS.type().intersection(vacuum_occupied_space(IS.qns()).type()) !=
           nulltype;
  }

  ///@brief some states are fermi vacuum unoccupied
  bool contains_unoccupied(const IndexSpace& IS) const {
    if (IS == nullspace) {
      return false;
    } else {
      return vacuum_occupied_space(IS.qns()).type() < IS.type();
    }
  }

  /// @name  specifies which spaces have nonzero occupancy in the vacuum wave
  ///        function
  /// @note needed for applying Wick theorem with Fermi vacuum
  /// @{

  /// @param t an IndexSpace::Type specifying which base spaces have nonzero
  ///          occupancy in
  ///          the vacuum wave function by default (i.e. for any quantum number
  ///          choice); to specify occupied space per specific QN set use the
  ///          other overload
  /// @return reference to `this`
  IndexSpaceRegistry& vacuum_occupied_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "vacuum_occupied_space");
    std::get<0>(vacocc_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  ///                nonzero occupancy in the reference wave function
  /// @return reference to `this`
  IndexSpaceRegistry& vacuum_occupied_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "vacuum_occupied_space");
    std::get<1>(vacocc_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `vacuum_occupied_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& vacuum_occupied_space(const IndexSpace& s) {
    return vacuum_occupied_space(s.type());
  }

  /// equivalent to `vacuum_occupied_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& vacuum_occupied_space(S&& l) {
    return vacuum_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in vacuum state for any set of quantum numbers
  const IndexSpace::Type& vacuum_occupied_space() const {
    if (std::get<0>(vacocc_) == nulltype) {
      throw std::invalid_argument(
          "vacuum occupied space has not been specified, invoke "
          "vacuum_occupied_space(IndexSpace::Type) or "
          "vacuum_occupied_space(std::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(vacocc_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space occupied in vacuum state for the given set of quantum
  /// numbers
  const IndexSpace& vacuum_occupied_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(vacocc_).find(qn);
    if (it != std::get<1>(vacocc_).end()) {
      return retrieve(it->second, qn);
    } else {
      return retrieve(this->vacuum_occupied_space(), qn);
    }
  }

  /// @}

  /// @name  assign which spaces have nonzero occupancy in the reference wave
  ///        function (i.e., the wave function uses to compute reference
  ///        expectation value)
  /// @note needed for computing expectation values when the vacuum state does
  /// not match the wave function of interest.
  /// @{

  /// @param t an IndexSpace::Type specifying which base spaces have nonzero
  /// occupancy in
  ///          the reference wave function by default (i.e., for any choice of
  ///          quantum numbers); to specify occupied space per specific QN set
  ///          use the other overload
  /// @return reference to `this`
  IndexSpaceRegistry& reference_occupied_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "reference_occupied_space");
    std::get<0>(refocc_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  /// nonzero occupancy in
  ///          the reference wave function
  /// @return reference to `this`
  IndexSpaceRegistry& reference_occupied_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "reference_occupied_space");
    std::get<1>(refocc_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `reference_occupied_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& reference_occupied_space(const IndexSpace& s) {
    return reference_occupied_space(s.type());
  }

  /// equivalent to `reference_occupied_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& reference_occupied_space(S&& l) {
    return reference_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in reference state for any set of quantum
  /// numbers
  const IndexSpace::Type& reference_occupied_space() const {
    if (std::get<0>(refocc_) == nulltype) {
      throw std::invalid_argument(
          "reference occupied space has not been specified, invoke "
          "reference_occupied_space(IndexSpace::Type) or "
          "reference_occupied_space(std::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(refocc_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space occupied in vacuum state for the given set of quantum
  /// numbers
  const IndexSpace& reference_occupied_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(refocc_).find(qn);
    if (it != std::get<1>(refocc_).end()) {
      return retrieve(it->second, qn);
    } else {
      return retrieve(this->vacuum_occupied_space(), qn);
    }
  }

  /// @}

  /// @name  specifies which spaces comprise the entirety of Hilbert space
  /// @note needed for creating general operators in mbpt/op
  /// @{

  /// @param t an IndexSpace::Type specifying the complete Hilbert space;
  ///          to specify occupied space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& complete_space(const IndexSpace::Type& s) {
    throw_if_missing(s, "complete_space");
    std::get<0>(complete_) = s;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  /// nonzero occupancy in
  ///          the reference wave function
  IndexSpaceRegistry& complete_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "complete_space");
    std::get<1>(complete_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `complete_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& complete_space(const IndexSpace& s) {
    return complete_space(s.type());
  }

  /// equivalent to `complete_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& complete_space(S&& l) {
    return complete_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the complete Hilbert space for any set of quantum numbers
  const IndexSpace::Type& complete_space() const {
    if (std::get<0>(complete_) == nulltype) {
      throw std::invalid_argument(
          "complete space has not been specified, call "
          "complete_space(IndexSpace::Type)");
    } else
      return std::get<0>(complete_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the complete Hilbert space for the given set of quantum numbers
  const IndexSpace& complete_space(const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(complete_).find(qn);
    if (it != std::get<1>(complete_).end()) {
      return retrieve(it->second, qn);
    } else {
      return retrieve(this->complete_space(), qn);
    }
  }

  /// @}

  /// @name specifies in which space holes can be created successfully from the
  /// reference wave function
  /// @note convenience for making operators
  /// @{

  /// @param t an IndexSpace::Type specifying where holes can be created;
  ///          to specify hole space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& active_hole_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "active_hole_space");
    std::get<0>(active_hole_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which holes
  /// can be created
  IndexSpaceRegistry& active_hole_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "active_hole_space");
    std::get<1>(active_hole_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `active_hole_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& active_hole_space(const IndexSpace& s) {
    return active_hole_space(s.type());
  }

  /// equivalent to `active_hole_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& active_hole_space(S&& l) {
    return active_hole_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which holes can be created
  const IndexSpace::Type& active_hole_space() const {
    if (std::get<0>(active_hole_space_) == nulltype) {
      throw std::invalid_argument(
          "active hole space has not been specified, invoke "
          "active_hole_space(IndexSpace::Type) or "
          "active_hole_space(std::map<IndexSpace::QuantumNumbers,IndexSpace::"
          "Type>)");
    } else
      return std::get<0>(active_hole_space_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space in which holes can be created for the given set of
  /// quantum numbers
  const IndexSpace& active_hole_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(active_hole_space_).find(qn);
    if (it != std::get<1>(active_hole_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->active_hole_space(), qn);
    }
  }

  /// @}

  /// @name specifies in which space particles can be created successfully from
  /// the reference wave function
  /// @note convenience for making operators
  /// @{

  /// @param t an IndexSpace::Type specifying where particles can be created;
  ///          to specify particle space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& active_particle_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "active_particle_space");
    std::get<0>(active_particle_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which
  /// particles can be created
  IndexSpaceRegistry& active_particle_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "active_particle_space");
    std::get<1>(active_particle_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `active_hole_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& active_particle_space(const IndexSpace& s) {
    return active_particle_space(s.type());
  }

  /// equivalent to `active_hole_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& active_particle_space(S&& l) {
    return active_particle_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which particles can be created
  const IndexSpace::Type& active_particle_space() const {
    if (std::get<0>(active_particle_space_) == nulltype) {
      throw std::invalid_argument(
          "active particle space has not been specified, invoke "
          "active_particle_space(IndexSpace::Type) or "
          "active_particle_space(std::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(active_particle_space_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space in which particles can be created for the given set of
  /// quantum numbers
  const IndexSpace& active_particle_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(active_particle_space_).find(qn);
    if (it != std::get<1>(active_particle_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->active_particle_space(), qn);
    }
  }

  /// @}

  const static IndexSpace nullspace;

 private:
  // N.B. need transparent comparator, see https://stackoverflow.com/a/35525806
  std::set<IndexSpace, IndexSpace::KeyCompare> spaces;

  const static IndexSpace::Type nulltype;

  /// TODO use c++20 std::has_single_bit() when we update to this version
  static bool has_single_bit(std::uint32_t bits) {
    return bits & (((bool)(bits & (bits - 1))) - 1);
  }

  ///@brief find an IndexSpace from its type. return nullspace if not present.
  ///@param find the IndexSpace via it's @c attr
  const IndexSpace& find_by_attr(const IndexSpace::Attr& attr) const {
    for (auto&& space : spaces) {
      if (space.attr() == attr) {
        return space;
      }
    }
    return nullspace;
  }

  void throw_if_missing(const IndexSpace::Type& t,
                        const IndexSpace::QuantumNumbers& qn,
                        std::string call_context = "") {
    for (auto&& space : spaces) {
      if (space.type() == t && space.qns() == qn) {
        return;
      }
    }
    throw std::invalid_argument(
        call_context +
        ": missing { IndexSpace::Type=" + std::to_string(t.to_int32()) +
        " , IndexSpace::QuantumNumbers=" + std::to_string(qn.to_int32()) +
        " } combination");
  }

  // same as above, but ignoring qn
  void throw_if_missing(const IndexSpace::Type& t,
                        std::string call_context = "") {
    for (auto&& space : spaces) {
      if (space.type() == t) {
        return;
      }
    }
    throw std::invalid_argument(call_context + ": missing { IndexSpace::Type=" +
                                std::to_string(t.to_int32()) +
                                " , any IndexSpace::QuantumNumbers } space");
  }

  void throw_if_missing_any(
      const std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>& qn2type,
      std::string call_context = "") {
    std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type_found;
    for (auto&& space : spaces) {
      for (auto&& [qn, t] : qn2type) {
        if (space.type() == t && space.qns() == qn) {
          auto [it, found] = qn2type_found.try_emplace(qn, t);
          assert(!found);
          // found all? return
          if (qn2type_found.size() == qn2type.size()) {
            return;
          }
        }
      }
    }

    std::string errmsg;
    for (auto&& [qn, t] : qn2type) {
      if (!qn2type_found.contains(qn)) {
        errmsg +=
            call_context +
            ": missing { IndexSpace::Type=" + std::to_string(t.to_int32()) +
            " , IndexSpace::QuantumNumbers=" + std::to_string(qn.to_int32()) +
            " } combination\n";
      }
    }
    throw std::invalid_argument(errmsg);
  }

  // Need to define defaults for various traits, like which spaces are occupied
  // in vacuum, etc. Makes sense to make these part of the registry to avoid
  // having to pass these around in every call N.B. default and QN-specific
  // space selections merged into single tuple

  // defines active bits in TypeAttr; used by general operators in mbpt/op
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      complete_ = {nulltype, {}};

  // used for fermi vacuum wick application
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      vacocc_ = {nulltype, {}};

  // used for MR MBPT to take average over multiconfiguration reference
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      refocc_ = {nulltype, {}};

  // both needed to make excitation and de-excitation operators. not
  // necessarily equivalent in the case of multi-reference context.
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      active_particle_space_ = {nulltype, {}};
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      active_hole_space_ = {nulltype, {}};

  friend bool operator==(const IndexSpaceRegistry& isr1,
                         const IndexSpaceRegistry& isr2) {
    return isr1.spaces == isr2.spaces;
  }
};

inline const IndexSpace IndexSpaceRegistry::nullspace{L"", 0, 0};
inline const IndexSpace::Type IndexSpaceRegistry::nulltype{0};

}  // namespace sequant
#endif  // SEQUANT_INDEX_SPACE_REGISTRY_HPP
