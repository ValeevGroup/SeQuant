//
// Created by Conner Masteran on 4/16/24.
//

#ifndef SEQUANT_INDEX_SPACE_REGISTRY_HPP
#define SEQUANT_INDEX_SPACE_REGISTRY_HPP

#include "space.hpp"

namespace sequant {

/// @brief a map of labels and @c IndexSpaces interpretable by the user.
/// An IndexSpaceRegistry is a context specific, modifiable, list of registered
/// IndexSpaces. It is the task of the USER to properly construct a complete set
/// for their needs. Thus, the IndexSpaceRegistry is the way a user can retrieve
/// IndexSpaces specific to their context. because set operations(overlap,
/// union) are context dependent, these functions must also live here.
/// @note unified spaces have connected @c typeattr.
/// @note  @c IndexSpaces with @c typeattr corresponding to occupied indices
/// will always be lower than typeattr corresponding to unoccupied orbitals.
class IndexSpaceRegistry {
 public:
  IndexSpaceRegistry() {
    // register nulltype
    this->add(nulltype);
  }

  /// @brief retrieve an IndexSpace from the registry by the label
  /// @param label can be numbered or the @c base_key
  /// @return IndexSpace associated with that key.
  /// @warning must provide spin decoration in label to access space with
  /// non-null @c QuantumNumbersAttr
  IndexSpace retrieve(std::wstring_view label) const {
    for (auto it = label_space.begin(); it != label_space.end(); ++it) {
      if (IndexSpace::reduce_key(label) == it->first) {
        return it->second;
      }
    }
    throw IndexSpace::bad_key(label);
  }

  /// @brief add an IndexSpace to this registry.
  /// @param IS a non-null IndexSpace
  /// @warning duplicate type bitsets forbidden
  void add(IndexSpace IS) {
    if (!(find_indexspace(IS) == nulltype)) {
      throw std::invalid_argument(
          "The registry already has a TypeAttr(bitset) corresponding to the "
          "IndexSpace you are trying to add! "
          "If you are trying to replace the index use the replace(IndexSpace) "
          "function");
    } else {
      label_space.emplace(IS.get_base_key(), IS);
    }
  }

  /// @brief all IndexSpace label pairs with one @c typeattr bit set ordered by
  /// increasing @c typeattr
  std::vector<std::pair<IndexSpace, std::wstring_view>> base_spaces_label() {
    std::vector<std::pair<IndexSpace, std::wstring_view>> result;
    // if the registered instance has a single bit true, it is a base space.
    //  std::map should order this new map via IndexSpace.
    for (auto begin = label_space.begin(); begin != label_space.end();
         begin++) {
      if (has_single_bit(begin->second.type().to_int32())) {
        result.push_back({begin->second, begin->first});
      }
    }
    std::sort(result.begin(), result.end(),
              [](std::pair<IndexSpace, std::wstring_view> a,
                 std::pair<IndexSpace, std::wstring_view> b) -> bool {
                return a.first.type() < b.first.type();
              });
    return result;
  }

  /// @brief does the space have only one bit set
  /// @param IS IndexSpace
  /// @note list of all base spaces is computed each time. Otherwise, a separate
  /// list would need to be maintained for every mutation of the registry. This
  /// seems safer.
  bool is_base_space(IndexSpace IS) {
    auto spaces_label = base_spaces_label();
    bool result = false;
    for (std::size_t i = 0; i < spaces_label.size(); i++) {
      if (IS == spaces_label[i].first) {
        result = true;
      }
    }
    return result;
  }

  /// @brief clear the IndexSpaceRegistry map
  void clear_registry() {
    label_space.clear();
    IndexSpace nulltype(L"", 0, 0);
    this->add(nulltype);
  }

  /// @brief remove an IndexSpace by its string label
  void remove(std::wstring_view label) {
    auto it = label_space.begin();
    bool found = false;
    // iterating through a container whose size changes requires special
    // handling
    while (it != label_space.end() && !found) {
      if (label == it->first) {
        label_space.erase(it);
        found = true;
      }
      if (it != label_space.end() && !found) {
        it++;
      }
    }
  }

  /// @brief replace the label of an existing IndexSpace in the registry
  /// @param original_label the original base_key of a space
  /// @param new_label the new chosen @c base_key
  void relabel(std::wstring original_label, std::wstring new_label) {
    bool found = false;
    auto it = label_space.begin();
    while (!found && it != label_space.end()) {
      if (original_label == it->first) {
        auto original_attr = it->second.attr();
        label_space.erase(it);
        IndexSpace new_space(new_label, original_attr.type(),
                             original_attr.qns());
        label_space.emplace(new_space.get_base_key(), new_space);
        found = true;
      }
      if (it != label_space.end() && !found) {
        it++;
      }
    }
    if (!found) {
      throw "orginal label not found!";
    }
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
  bool valid_bitop(const IndexSpace &i1, const IndexSpace &i2,
                   const std::function<int32_t(int32_t, int32_t)> &op) {
    auto bitop_int = op(i1.type().to_int32(), i2.type().to_int32());
    bool same_qn = i1.qns() == i2.qns();
    if (!same_qn) return false;
    auto temp_space = find_indexspace({bitop_int, i1.qns()});
    return temp_space == nulltype ? false : true;
  }

  ///@brief return the resulting space corresponding to a bitwise intersection
  /// between two spaces.
  /// @param space1
  /// @param space2
  /// @return the resulting space after intesection
  /// @note nulltype_() is a valid return for intersection
  /// @note throw invalid_argument if the bitwise result is not registered
  const IndexSpace intersection(const IndexSpace &space1,
                                const IndexSpace &space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {  // spaces with different quantum numbers do not
                        // intersect.
        return nulltype_();
      }
      auto intersection_attr = space1.type().intersection(space2.type());
      IndexSpace intersection_space =
          find_indexspace({intersection_attr, space1.qns()});
      // the nullspace is a reasonable return value for intersection
      if (intersection_space == nulltype && intersection_attr != 0) {
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
  /// @note nulltype_() is a valid return for intersection
  /// @note throw invalid_argument if the bitwise result is not registered
  IndexSpace intersection(const IndexSpace &space1, const IndexSpace &space2,
                          const IndexSpace &space3) {
    if (space1 == space2 && space1 == space3) {
      return space1;
    } else {
      bool same_qns =
          ((space1.qns() == space2.qns()) && (space1.qns() == space3.qns()));
      if (!same_qns) {  // spaces with different quantum numbers do not
                        // intersect.
        return nulltype_();
      }
      auto intersection_attr =
          space1.type().intersection(space2.type()).intersection(space3.type());
      IndexSpace intersection_space =
          find_indexspace({intersection_attr, space1.qns()});
      if (intersection_space == nulltype && intersection_attr != 0) {
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
  /// @note nulltype_() is not a valid return for this operation.
  IndexSpace unIon(const IndexSpace &space1, const IndexSpace &space2) {
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
      IndexSpace unIonSpace = find_indexspace({unIontype, space1.qns()});
      if (unIonSpace == nulltype) {
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
  /// @note nulltype is not a valid return
  /// @note finding unregistered @c typeattr will throw
  std::vector<IndexSpace> non_overlapping_spaces(
      const IndexSpace &space1, const IndexSpace &space2) const {
    auto attributes = space1.attr().excluded_spaces(space2.attr());
    std::vector<IndexSpace> result;
    for (int i = 0; i < attributes.size(); i++) {
      auto excluded_space = find_indexspace(attributes[i]);
      if (excluded_space == nulltype) {
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
  bool has_non_overlapping_spaces(const IndexSpace &space1,
                                  const IndexSpace &space2) const {
    return space1.type().exclusionary_or(space2.type()).to_int32() != 0;
  }
  IndexSpace nulltype_() const { return nulltype; }

  ///@brief an @c IndexSpace is occupied with respect to the fermi vacuum or a
  /// subset of that space
  /// @note only makes sense to ask this if in a SingleProduct vacuum context.
  bool is_pure_occupied(IndexSpace IS) const {
    if (IS == nulltype_()) {
      return false;
    }
    if (IS.type().to_int32() <= vacuum_occupied().type().to_int32()) {
      return true;
    } else {
      return false;
    }
  }

  ///@brief all states are unoccupied in the fermi vacuum
  ///@note again, this only makes sense to ask if in a SingleProduct vacuum
  /// context.
  bool is_pure_unoccupied(IndexSpace IS) const {
    if (IS == nulltype_()) {
      return false;
    } else {
      return intersection(IS, vacuum_occupied()) == nulltype_();
    }
  }

  ///@brief some states are fermi vacuum occupied
  bool contains_occupied(IndexSpace IS) {
    return this->intersection(IS, vacuum_occupied()) != nulltype_();
  }

  ///@brief some states are fermi vacuum unoccupied
  bool contains_unoccupied(IndexSpace IS) const {
    if (IS == nulltype_()) {
      return false;
    } else {
      return vacuum_occupied().type() < IS.type();
    }
  }

  /// @brief complete space is the unIon of all base spaces
  /// @note useful for operator construction
  /// @note complete space must be assigned
  IndexSpace complete() {
    if (complete_ == nulltype) {
      throw std::invalid_argument(
          "complete has not been registered. please assign the complete"
          "space label by calling assign_complete(label)");
    } else
      return complete_;
  }

  /// @brief union of all fermi occupied base spaces
  /// @note essential for SingleProduct vacuum wick application.
  /// @note this space must be assigned by the user
  IndexSpace vacuum_occupied() const {
    if (vacuum_occupied_ == nulltype) {
      throw std::invalid_argument(
          "occupied has not been registered. please assign the occupied"
          "space label by calling assign_vacuum_occupied(label)");
    } else
      return vacuum_occupied_;
  }

  /// @brief space where active particles can reside
  /// @note necessary for operator construction
  /// @note allows user to explicitly decide what an excitation or de-excitation
  /// looks like in their context
  IndexSpace active_particle_space() const {
    if (active_particle_space_ == nulltype) {
      throw std::invalid_argument(
          "active_particle_space has not been registered. please assign the"
          "space label by calling assign_active_particle_space(label)");
    } else
      return active_particle_space_;
  }

  /// @brief space where active holes can reside
  /// @note necessary for operator construction
  /// @note allows user to explicitly decide what an excitation or de-excitation
  /// looks like in their context
  IndexSpace active_hole_space() const {
    if (active_hole_space_ == nulltype) {
      throw std::invalid_argument(
          "active_hole_space has not been registered. please assign the"
          "space label by calling assign_hole_particle_space(label)");
    } else
      return active_hole_space_;
  }

  /// @brief which spaces have density in the wave function
  /// @note needed for computing expectation values when the vacuum state does
  /// not match the wave function of interest.
  IndexSpace density_occupied() const {
    if (density_occuiped_ == nulltype_()) {
      throw std::invalid_argument(
          "density_occupied has not been registered. please assign the"
          "space label by calling assign_density_occupied(label)");
    } else {
      return density_occuiped_;
    }
  }
  /// @brief assign a space to set the fermi vacuum occupancy
  /// @param label the base_key of an IndexSpace in the registry.
  /// @warning The user must be smart about this! proper wick application
  /// requires this to be set properly!
  void assign_vacuum_occupied(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    vacuum_occupied_ = label_space[label];
  }

  /// @brief assign the space used to represent generic mbpt::Operators.
  /// @param label the base_key of an IndexSpace in the registry.
  void assign_complete(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    complete_ = label_space[label];
  }

  ///@brief assign a space where particles can exist and be displaced from
  /// @param label the base_key of an IndexSpace in the registry.
  ///@note this concept is useful when constructing fermionic excitation and
  /// de-excitation operators.
  ///@note the choice here will effect the behavior of SeQuant's provided
  /// operators in mbpt::op and mbpt::TensorOp.
  void assign_active_particle_space(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    active_particle_space_ = label_space[label];
  }

  ///@brief assign a space where holes can exist and be displaced from
  /// @param label the base_key of an IndexSpace in the registry.
  ///@note this concept is useful when constructing fermionic excitation and
  /// de-excitation operators.
  ///@note the choice here will effect the behavior of SeQuant's provided
  /// operators in mbpt::op and mbpt::TensorOp.
  void assign_active_hole_space(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    active_hole_space_ = label_space[label];
  }

  /// @brief assign which spaces have density in the wave function
  /// @param label the base_key of an IndexSpace in the registry.
  /// @note needed for computing expectation values when the vacuum state does
  /// not match the wave function of interest.
  void assign_density_occupied(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    density_occuiped_ = label_space[label];
  }

  /// @return a copy of the label_space map which defines this registry.
  std::map<std::wstring, IndexSpace> get_map() const { return label_space; }

 private:
  std::map<std::wstring, IndexSpace> label_space;
  const IndexSpace nulltype = {L"", 0, 0};

  /// TODO use c++20 std::has_single_bit() when we update to this version
  bool has_single_bit(std::uint32_t bits) {
    return bits & (((bool)(bits & (bits - 1))) - 1);
  }
  ///@brief find an IndexSpace from its type. return nullspace if not present.
  ///@param IS find via the IndexSpace
  const IndexSpace find_indexspace(IndexSpace IS) const {
    for (auto it = label_space.begin(); it != label_space.end(); it++) {
      if (it->second.attr() == IS.get_attr()) {
        return it->second;
      }
    }
    return nulltype;
  }

  ///@brief find an IndexSpace from its type. return nullspace if not present.
  ///@param find the IndexSpace via it's @c attr
  const IndexSpace find_indexspace(IndexSpace::Attr attr) const {
    for (auto it = label_space.begin(); it != label_space.end(); it++) {
      if (it->second.attr() == attr) {
        return it->second;
      }
    }
    return nulltype;
  }
  IndexSpace complete_ = {L"", 0,
                          0};  // need to define to use generic operator class
  IndexSpace vacuum_occupied_ = {
      L"", 0, 0};  // needed for fermi vacuum wick application
  IndexSpace density_occuiped_ = {L"", 0, 0};
  // both needed to make excitation and de-excitation operators. not
  // necessarily equivalent in the case of multi-reference context.
  IndexSpace active_particle_space_ = {L"", 0, 0};
  IndexSpace active_hole_space_ = {L"", 0, 0};
};

inline bool operator==(const IndexSpaceRegistry &isr1,
                       const IndexSpaceRegistry &isr2) {
  return isr1.get_map() == isr2.get_map();
}
}  // namespace sequant
#endif  // SEQUANT_INDEX_SPACE_REGISTRY_HPP
