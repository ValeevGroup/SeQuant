//
// Created by Conner Masteran on 4/16/24.
//

#ifndef SEQUANT_INDEX_SPACE_REGISTRY_H
#define SEQUANT_INDEX_SPACE_REGISTRY_H

#include "space.hpp"

namespace sequant {

// An IndexSpaceRegistry is a context specific, modifiable, list of registered
// IndexSpaces. It is the task of the USER to properly construct a complete set
// for their needs. Thus, the IndexSpaceRegistry is the way a user can retrieve
// IndexSpaces specific to their context. because set operations(overlap, union)
// are context dependent, these functions must also live here. Currently, there
// are two important restrictions, 1. unified spaces must be connected (1010
// invalid typeattr).
// 2. IndexSpaces with typeattr corresponding to occupied indices will always be
// lower than typeattr corresponding to unoccupied orbitals.  In the future,
// these restrictions may be lifted if a need arises.
//
class IndexSpaceRegistry {
 public:
  IndexSpaceRegistry() {
    // register nulltype
    this->add(nulltype);
  }

  // retrieve an IndexSpace from the registry by the label
  IndexSpace retrieve(std::wstring_view label) const {
    for (auto it = label_space.begin(); it != label_space.end(); ++it) {
      if (IndexSpace::reduce_key(label) == it->first) {
        return it->second;
      }
    }
    throw IndexSpace::bad_key(label);
  }

  // add an IndexSpace to this registry. duplicate type bitsets forbidden
  void add(IndexSpace IS) {
    if (!(find_IndexSpace(IS) == nulltype)) {
      throw std::invalid_argument(
          "The registry already has a TypeAttr(bitset) corresponding to the "
          "IndexSpace you are trying to add! "
          "If you are trying to replace the index use the replace(IndexSpace) "
          "function");
    } else {
      label_space.emplace(IS.get_base_key(), IS);
    }
  }

  // return a vector of pairs since we will be accessing by index not key.
  std::vector<std::pair<IndexSpace, std::wstring_view>> base_spaces_label() {
    std::vector<std::pair<IndexSpace, std::wstring_view>> result;
    // if the registered instance has a single bit true, it is a base space.
    //  std::map should order this new map via IndexSpace.
    /// TODO check that this is ordered correctly
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

  // clear the label_space map essentially clearing the registry
  void clear_registry() {
    label_space.clear();
    IndexSpace nulltype(L"", 0, 0);
    this->add(nulltype);
  }

  // remove an IndexSpace by its string label
  void remove(std::wstring_view label) {
    auto it = label_space.begin();
    bool found = false;
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

  // replace a member of the registry by passing the original label and the new
  // space to replace it
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

  // pass a function which computes a logical bit operation between two
  // IndexSpace.type()
  bool valid_bitop(const IndexSpace &i1, const IndexSpace &i2,
                   const std::function<int32_t(int32_t, int32_t)> &op) {
    auto bitop_int = op(i1.type().to_int32(), i2.type().to_int32());
    bool same_qn = i1.qns() == i2.qns();
    if (!same_qn) return false;
    auto temp_space = find_IndexSpace({bitop_int, i1.qns()});
    return temp_space == nulltype ? false : true;
  }

  // return the resulting space corresponding to a bitwise intersection between
  // two spaces.
  //  check to see if the resulting space is actually registered.
  const IndexSpace intersection(const IndexSpace &space1,
                                const IndexSpace &space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {
        throw std::invalid_argument(
            "asking for the intersection of spaces with incompatible quantum "
            "number attributes.");
      }
      auto intersection_attr = space1.type().intersection(space2.type());
      IndexSpace intersection_space =
          find_IndexSpace({intersection_attr, space1.qns()});
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

  // return the resulting space corresponding to a bitwise intersection between
  // three spaces.
  //  check to see if the resulting space is actually registered.
  IndexSpace intersection(const IndexSpace &space1, const IndexSpace &space2,
                          const IndexSpace &space3) {
    if (space1 == space2 && space1 == space3) {
      return space1;
    } else {
      bool same_qns =
          ((space1.qns() == space2.qns()) && (space1.qns() == space3.qns()));
      if (!same_qns) {
        throw std::invalid_argument(
            "asking for the intersection of spaces with incompatible quantum "
            "number attributes.");
      }
      auto intersection_attr =
          space1.type().intersection(space2.type()).intersection(space3.type());
      IndexSpace intersection_space =
          find_IndexSpace({intersection_attr, space1.qns()});
      if (intersection_space == nulltype && intersection_attr != 0) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return intersection_space;
      }
    }
  }

  // return the resulting space corresponding to a bitwise unIon between two
  // spaces.
  //  check to see if the resulting space is actually registered.
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
      IndexSpace unIonSpace = find_IndexSpace({unIontype, space1.qns()});
      if (unIonSpace == nulltype) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return unIonSpace;
      }
    }
  }

  // return the resulting spaces corresponding to a bitwise excluded or between
  // two spaces.
  //  check to see if the resulting space is actually registered.
  std::vector<IndexSpace> non_overlapping_spaces(
      const IndexSpace &space1, const IndexSpace &space2) const {
    auto attributes = space1.attr().excluded_spaces(space2.attr());
    std::vector<IndexSpace> result;
    for (int i = 0; i < attributes.size(); i++) {
      auto excluded_space = find_IndexSpace(attributes[i]);
      if (excluded_space == nulltype) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      }
      result.push_back(excluded_space);
    }
    return result;
  }

  // do two spaces have non_overlapping spaces
  bool has_non_overlapping_spaces(const IndexSpace &space1,
                                  const IndexSpace &space2) const {
    return space1.type().exclusionary_or(space2.type()).to_int32() != 0;
  }
  IndexSpace nulltype_() const { return nulltype; }

  // pure occupied means all states are occupied in a given reference
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

  // all states are unoccupied with respect to a given reference
  bool is_pure_unoccupied(IndexSpace IS) const {
    if (IS == nulltype_()) {
      return false;
    } else {
      return intersection(IS, vacuum_occupied()) == nulltype_();
    }
  }

  // some states are occupied in a given reference
  bool contains_occupied(IndexSpace IS) {
    return this->intersection(IS, vacuum_occupied()) != nulltype_();
  }

  // some states are unoccupied in a given reference
  bool contains_unoccupied(IndexSpace IS) const {
    if (IS == nulltype_()) {
      return false;
    } else {
      return vacuum_occupied().type() < IS.type();
    }
  }

  // complete space with all possible IndexSpaces included
  IndexSpace complete() {
    if (complete_ == nulltype) {
      throw std::invalid_argument(
          "complete has not been registered. please assign the complete"
          "space label by calling assign_complete(label)");
    } else
      return complete_;
  }

  // occupied with respect to vacuum reference!! defining the occupied orbitals
  // in the Fermi vacuum is absolutely essential.
  IndexSpace vacuum_occupied() const {
    if (vacuum_occupied_ == nulltype) {
      throw std::invalid_argument(
          "occupied has not been registered. please assign the occupied"
          "space label by calling assign_vacuum_occupied(label)");
    } else
      return vacuum_occupied_;
  }

  // space where active particles can reside, in Single Reference case this is
  // usually active_occupied. Multireference will also include an active space.
  IndexSpace active_particle_space() const {
    if (active_particle_space_ == nulltype) {
      throw std::invalid_argument(
          "active_particle_space has not been registered. please assign the"
          "space label by calling assign_active_particle_space(label)");
    } else
      return active_particle_space_;
  }

  // space where all active holes can reside, in Single Reference case this is
  // usually active_unoccupied. Multireference will also include an active
  // space.
  IndexSpace active_hole_space() const {
    if (active_hole_space_ == nulltype) {
      throw std::invalid_argument(
          "active_hole_space has not been registered. please assign the"
          "space label by calling assign_hole_particle_space(label)");
    } else
      return active_hole_space_;
  }

  // needed to compute densites in physical vacuum.
  IndexSpace density_occupied() const {
    if (density_occuiped_ == nulltype_()) {
      throw std::invalid_argument(
          "density_occupied has not been registered. please assign the"
          "space label by calling assign_density_occupied(label)");
    } else {
      return density_occuiped_;
    }
  }
  void assign_vacuum_occupied(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    vacuum_occupied_ = label_space[label];
  }

  void assign_complete(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    complete_ = label_space[label];
  }

  void assign_active_particle_space(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    active_particle_space_ = label_space[label];
  }

  void assign_active_hole_space(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    active_hole_space_ = label_space[label];
  }

  // which spaces could contain particles in your representation
  void assign_density_occupied(std::wstring label) {
    if (label_space.find(label) == label_space.end()) {
      throw std::invalid_argument("label not added to registry");
    }
    density_occuiped_ = label_space[label];
  }
  std::map<std::wstring, IndexSpace> get_map() const { return label_space; }

 private:
  std::map<std::wstring, IndexSpace> label_space;
  const IndexSpace nulltype = {L"", 0, 0};

  /// TODO use c++20 std::has_single_bit() when we update to this version
  bool has_single_bit(std::uint32_t bits) {
    return bits & (((bool)(bits & (bits - 1))) - 1);
  }
  // find an indexspace from its type. return nullspace if not present.
  // a bit strange, but prevents an additional map that needs to be maintained
  const IndexSpace find_IndexSpace(IndexSpace IS) const {
    for (auto it = label_space.begin(); it != label_space.end(); it++) {
      if (it->second.attr() == IS.get_attr()) {
        return it->second;
      }
    }
    return nulltype;
  }

  // find an IndexSpace from attribute. return nullspace if not present.
  // sometimes we wish to check whether if an antribute is in the registry
  const IndexSpace find_IndexSpace(IndexSpace::Attr attr) const {
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
  // neccessarily exclusionary in the case of multi-reference context.
  IndexSpace active_particle_space_ = {L"", 0, 0};
  IndexSpace active_hole_space_ = {L"", 0, 0};
};

inline bool operator==(const IndexSpaceRegistry &isr1,
                       const IndexSpaceRegistry &isr2) {
  return isr1.get_map() == isr2.get_map();
}
}  // namespace sequant
#endif  // SEQUANT_INDEX_SPACE_REGISTRY_H
