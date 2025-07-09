//
// Created by Conner Masteran on 4/16/24.
//

#ifndef SEQUANT_INDEX_SPACE_REGISTRY_HPP
#define SEQUANT_INDEX_SPACE_REGISTRY_HPP

#include <SeQuant/core/space.hpp>

#include <range/v3/algorithm/sort.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/unique.hpp>

#include <boost/hana.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>

#include <mutex>

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

// clang-format off
/// @brief set of known IndexSpace objects

/// Each IndexSpace object has hardwired base key (label) that gives
/// indexed expressions appropriate semantics; e.g., spaces referred to by
/// indices in \f$ t_{p_1}^{i_1} \f$ are defined if IndexSpace objects with
/// base keys \f$ p \f$ and \f$ i \f$ are registered.
/// Since index spaces have set-theoretic semantics, the user must
/// provide complete set of unions/intersects of the base spaces to
/// cover all possible IndexSpace objects that can be generated in their
/// program.
///
/// Registry contains 2 parts: set of IndexSpace objects (managed by a
/// `std::shared_ptr`, see IndexSpaceRegistry::spaces()) and specification of
/// various spaces (vacuum, reference, complete, etc.). Copy semantics is thus
/// partially shallow, with spaces shared between copies. This allows to have
/// multiple registries share same set of spaces but have different
/// specifications of vacuum, reference, etc.; this is useful for providing
/// different contexts for fermions and bosons, for example.
///
/// Spaces that can be occupied by physical particles need to be
/// introspected for their structure, occupancy, etc. The registry provides the
/// API needed for dealing with such states.
/// - IndexSpaceRegistry::physical_particle_attribute_mask specify which states are occupied by
/// physical particles
/// - every space occupied by physical particles is a union of basis
/// ("base") spaces. IndexSpaceRegistry::is_base detects such spaces and
/// IndexSpaceRegistry::base_space_types/IndexSpaceRegistry::base_spaces report the list of base spaces
/// - in SingleProduct vacuum IndexSpaceRegistry::is_pure_occupied/IndexSpaceRegistry::is_pure_unoccupied report
/// whether a space is occupied or unoccupied (always the case for base
/// spaces, but neither may be true for composite  spaces)
/// - IndexSpaceRegistry::vacuum_occupied_space report whether a space has nonzero occupancy
/// in the vacuum state (that defines the normal order); this is needed to
/// apply Wick theorem with Fermi vacuum
/// - IndexSpaceRegistry::reference_occupied_space reports whether a space has nonzero occupancy
/// in the reference state used to compute reference
/// expectation value; only needed for computing expectation values when the
/// vacuum state does not match the expectation value state.
/// - IndexSpaceRegistry::complete_space specifies which spaces comprise the entirety of Hilbert
/// space; needed for creating general operators in mbpt/op
/// - IndexSpaceRegistry::particle_space and IndexSpaceRegistry::hole_space specify in which space particles/holes
/// can be created successfully from the reference state; this is a
/// convenience for making operators
// clang-format on
class IndexSpaceRegistry {
 public:
  /// default constructor creates a registry containing only IndexSpace::null
  /// @note null space is registered so we don't have to handle it as a corner
  /// case in retrieve() and other methods
  IndexSpaceRegistry()
      : spaces_(std::make_shared<
                container::set<IndexSpace, IndexSpace::KeyCompare>>()) {
    // register nullspace
    this->add(IndexSpace::null);
  }

  /// constructs an IndexSpaceRegistry from an existing set of IndexSpace
  /// objects
  IndexSpaceRegistry(
      std::shared_ptr<container::set<IndexSpace, IndexSpace::KeyCompare>>
          spaces)
      : spaces_(std::move(spaces)) {}

  /// copy constructor
  IndexSpaceRegistry(const IndexSpaceRegistry& other)
      : spaces_(other.spaces_),
        physical_particle_attribute_mask_(
            other.physical_particle_attribute_mask_),
        vacocc_(other.vacocc_),
        refocc_(other.refocc_),
        complete_(other.complete_),
        hole_space_(other.hole_space_),
        particle_space_(other.particle_space_) {}

  /// move constructor
  IndexSpaceRegistry(IndexSpaceRegistry&& other)
      : spaces_(std::move(other.spaces_)),
        physical_particle_attribute_mask_(
            std::move(other.physical_particle_attribute_mask_)),
        vacocc_(std::move(other.vacocc_)),
        refocc_(std::move(other.refocc_)),
        complete_(std::move(other.complete_)),
        hole_space_(std::move(other.hole_space_)),
        particle_space_(std::move(other.particle_space_)) {}

  /// copy assignment operator
  IndexSpaceRegistry& operator=(const IndexSpaceRegistry& other) {
    spaces_ = other.spaces_;
    physical_particle_attribute_mask_ = other.physical_particle_attribute_mask_;
    vacocc_ = other.vacocc_;
    refocc_ = other.refocc_;
    complete_ = other.complete_;
    hole_space_ = other.hole_space_;
    particle_space_ = other.particle_space_;
    return *this;
  }

  /// move assignment operator
  IndexSpaceRegistry& operator=(IndexSpaceRegistry&& other) {
    spaces_ = std::move(other.spaces_);
    physical_particle_attribute_mask_ =
        std::move(other.physical_particle_attribute_mask_);
    vacocc_ = std::move(other.vacocc_);
    refocc_ = std::move(other.refocc_);
    complete_ = std::move(other.complete_);
    hole_space_ = std::move(other.hole_space_);
    particle_space_ = std::move(other.particle_space_);
    return *this;
  }

  /// deep copy of this object, creates a copy of its spaces
  IndexSpaceRegistry clone() const;

  const auto& spaces() const { return spaces_; }

  decltype(auto) begin() const { return spaces_->cbegin(); }
  decltype(auto) end() const { return spaces_->cend(); }

  /// @brief retrieve a pointer to IndexSpace from the registry by the label
  /// @param label a @c base_key of an IndexSpace, or a label of an Index (see
  /// Index::label() )
  /// @return pointer to IndexSpace associated with that key, or nullptr if not
  /// found
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  const IndexSpace* retrieve_ptr(S&& label) const {
    auto it =
        spaces_->find(IndexSpace::reduce_key(to_basic_string_view(label)));
    return it != spaces_->end() ? &(*it) : nullptr;
  }

  /// @brief retrieve a pointer to IndexSpace from the registry by the label
  /// @param label a @c base_key of an IndexSpace, or a label of an Index (see
  /// Index::label() )
  /// @return pointer to IndexSpace associated with that key, or nullptr if not
  /// found
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpace* retrieve_ptr(S&& label) {
    auto it =
        spaces_->find(IndexSpace::reduce_key(to_basic_string_view(label)));
    return it != spaces_->end() ? &(*it) : nullptr;
  }

  /// @brief retrieve an IndexSpace from the registry by the label
  /// @param label a @c base_key of an IndexSpace, or a label of an Index (see
  /// Index::label() )
  /// @return IndexSpace associated with that key
  /// @throw IndexSpace::bad_key if matching space is not found
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  const IndexSpace& retrieve(S&& label) const {
    if (const auto* ptr = retrieve_ptr(std::forward<S>(label))) {
      return *ptr;
    } else
      throw IndexSpace::bad_key(label);
  }

  /// @brief retrieve a pointer to IndexSpace from the registry by its type and
  /// quantum numbers
  /// @param type IndexSpace::Type
  /// @param qns IndexSpace::QuantumNumbers
  /// @return pointer to the IndexSpace associated with that key, or nullptr if
  /// not found
  const IndexSpace* retrieve_ptr(const IndexSpace::Type& type,
                                 const IndexSpace::QuantumNumbers& qns) const {
    auto it = std::find_if(
        spaces_->begin(), spaces_->end(),
        [&](const auto& is) { return is.type() == type && is.qns() == qns; });
    return it != spaces_->end() ? &(*it) : nullptr;
  }

  /// @brief retrieve an IndexSpace from the registry by its type and quantum
  /// numbers
  /// @param type IndexSpace::Type
  /// @param qns IndexSpace::QuantumNumbers
  /// @return IndexSpace associated with that key.
  /// @throw std::invalid_argument if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Type& type,
                             const IndexSpace::QuantumNumbers& qns) const {
    if (const auto* ptr = retrieve_ptr(type, qns)) {
      return *ptr;
    } else
      throw std::invalid_argument(
          "IndexSpaceRegistry::retrieve(type,qn): missing { IndexSpace::Type=" +
          std::to_string(type.to_int32()) + " , IndexSpace::QuantumNumbers=" +
          std::to_string(qns.to_int32()) + " } combination");
  }

  /// @brief retrieve pointer to the IndexSpace from the registry by the
  /// IndexSpace::Attr
  /// @param space_attr an IndexSpace::Attr
  /// @return pointer to the IndexSpace associated with that key, or nullptr if
  /// not found
  const IndexSpace* retrieve_ptr(const IndexSpace::Attr& space_attr) const {
    auto it = std::find_if(
        spaces_->begin(), spaces_->end(),
        [&space_attr](const IndexSpace& s) { return s.attr() == space_attr; });
    return it != spaces_->end() ? &(*it) : nullptr;
  }

  /// @brief retrieve an IndexSpace from the registry by the IndexSpace::Attr
  /// @param space_attr an IndexSpace::Attr
  /// @return IndexSpace associated with that key.
  /// @throw std::invalid_argument if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Attr& space_attr) const {
    if (const auto* ptr = retrieve_ptr(space_attr)) {
      return *ptr;
    } else
      throw std::invalid_argument(
          "IndexSpaceRegistry::retrieve(attr): missing { IndexSpace::Type=" +
          std::to_string(space_attr.type().to_int32()) +
          " , IndexSpace::QuantumNumbers=" +
          std::to_string(space_attr.qns().to_int32()) + " } combination");
  }

  /// queries presence of a registered IndexSpace
  /// @param label a @c base_key of an IndexSpace, or a label of an Index (see
  /// Index::label() )
  /// @return true, if an IndexSpace with key @p label is registered
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  bool contains(S&& label) const {
    return this->retrieve_ptr(std::forward<S>(label));
  }

  /// queries presence of a registered IndexSpace
  /// @param space an IndexSpace object
  /// @return true, if an IndexSpace with key `{type,qns}` is registered
  bool contains(const IndexSpace& space) const {
    return this->retrieve_ptr(space.type(), space.qns());
  }

  /// queries presence of a registered IndexSpace
  /// @param type an IndexSpace::Type object
  /// @param qns an IndexSpace::QuantumNumbers object
  /// @return true, if an IndexSpace with key `{type,qns}` is registered
  bool contains(const IndexSpace::Type& type,
                const IndexSpace::QuantumNumbers& qns) const {
    return this->retrieve_ptr(type, qns);
  }

  /// queries presence of a registered IndexSpace
  /// @param space_attr an IndexSpace::Attr object
  /// @return true, if an IndexSpace with key @p space_attr is registered
  bool contains(const IndexSpace::Attr& space_attr) const {
    return this->retrieve_ptr(space_attr);
  }

  /// @name adding IndexSpace objects to the registry
  /// @{

  /// @brief add an IndexSpace to this registry.
  /// @param IS an IndexSpace
  /// @return reference to `this`
  /// @throw std::invalid_argument if `IS.base_key()` or `IS.attr()` matches
  /// an already registered IndexSpace
  IndexSpaceRegistry& add(const IndexSpace& IS) {
    auto it = spaces_->find(IS.base_key());
    if (it != spaces_->end()) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add(is): already have an IndexSpace associated "
          "with is.base_key(); if you are trying to replace the IndexSpace use "
          "IndexSpaceRegistry::replace(is)");
    } else {
      // make sure there are no duplicate IndexSpaces whose attribute is
      // IS.attr()
      if (ranges::any_of(*spaces_,
                         [&IS](auto&& is) { return IS.attr() == is.attr(); })) {
        throw std::invalid_argument(
            "IndexSpaceRegistry::add(is): already have an IndexSpace "
            "associated with is.attr(); if you are trying to replace the "
            "IndexSpace use IndexSpaceRegistry::replace(is)");
      }
      spaces_->emplace(IS);
    }

    return clear_memoized_data_and_return_this();
  }

  /// @brief add an IndexSpace to this registry.
  /// @param type_label a label that will denote the space type,
  ///                   must be convertible to a std::string
  /// @param type an IndexSpace::Type
  /// @param args optional arguments consisting of a mix of zero or more of
  /// the following:
  ///   - IndexSpace::QuantumNumbers
  ///   - approximate size of the space (unsigned long)
  ///   - any of { is_vacuum_occupied , is_reference_occupied , is_complete ,
  ///   is_hole , is_particle }
  /// @return reference to `this`
  /// @throw std::invalid_argument if `type_label` or `type` matches
  /// an already registered IndexSpace
  template <typename S, typename... OptionalArgs,
            typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpaceRegistry& add(S&& type_label, IndexSpace::Type type,
                          OptionalArgs&&... args) {
    auto h_args = boost::hana::make_tuple(args...);

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
    IndexSpace::QuantumNumbers qns;
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

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// @brief add a union of IndexSpace objects to this registry.
  /// @param type_label a label that will denote the space type,
  ///                   must be convertible to a std::string
  /// @param components sequence of IndexSpace objects or labels (known to this)
  /// whose union will be known by @p type_label
  /// @param args optional arguments consisting of a mix of zero or more of
  /// { is_vacuum_occupied , is_reference_occupied , is_complete , is_hole ,
  /// is_particle }
  /// @return reference to `this`
  template <typename S, typename IndexSpaceOrLabel, typename... OptionalArgs,
            typename = meta::EnableIfAllBasicStringConvertible<S>,
            typename = std::enable_if_t<
                (std::is_same_v<std::decay_t<IndexSpaceOrLabel>, IndexSpace> ||
                 meta::is_basic_string_convertible_v<
                     std::decay_t<IndexSpaceOrLabel>>)>>
  IndexSpaceRegistry& add_unIon(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    assert(components.size() > 1);

    auto h_args = boost::hana::make_tuple(args...);

    // make space
    IndexSpace::Attr space_attr;
    long count = 0;
    if (components.size() <= 1) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add_union: must have at least two components");
    }
    for (auto&& component : components) {
      const IndexSpace* component_ptr;
      if constexpr (std::is_same_v<std::decay_t<IndexSpaceOrLabel>,
                                   IndexSpace>) {
        component_ptr = &component;
      } else {
        component_ptr = &(this->retrieve(component));
      }
      if (count == 0)
        space_attr = component_ptr->attr();
      else
        space_attr = space_attr.unIon(component_ptr->attr());
      ++count;
    }
    const auto approximate_size = compute_approximate_size(space_attr);

    IndexSpace space(std::forward<S>(type_label), space_attr.type(),
                     space_attr.qns(), approximate_size);
    this->add(space);
    auto type = space.type();

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// alias to add_unIon
  template <typename S, typename IndexSpaceOrLabel, typename... OptionalArgs,
            typename = meta::EnableIfAllBasicStringConvertible<S>,
            typename = std::enable_if_t<
                (std::is_same_v<std::decay_t<IndexSpaceOrLabel>, IndexSpace> ||
                 meta::is_basic_string_convertible_v<
                     std::decay_t<IndexSpaceOrLabel>>)>>
  IndexSpaceRegistry& add_union(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    return this->add_unIon(std::forward<S>(type_label), components,
                           std::forward<OptionalArgs>(args)...);
  }

  /// @brief add a union of IndexSpace objects to this registry.
  /// @param type_label a label that will denote the space type,
  ///                   must be convertible to a std::string
  /// @param components sequence of IndexSpace objects or labels (known to this)
  /// whose intersection will be known by @p type_label
  /// @param args optional arguments consisting of a mix of zero or more of
  /// { is_vacuum_occupied , is_reference_occupied , is_complete , is_hole ,
  /// is_particle }
  /// @return reference to `this`
  template <typename S, typename IndexSpaceOrLabel, typename... OptionalArgs,
            typename = meta::EnableIfAllBasicStringConvertible<S>,
            typename = std::enable_if_t<
                (std::is_same_v<std::decay_t<IndexSpaceOrLabel>, IndexSpace> ||
                 meta::is_basic_string_convertible_v<
                     std::decay_t<IndexSpaceOrLabel>>)>>
  IndexSpaceRegistry& add_intersection(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    assert(components.size() > 1);

    auto h_args = boost::hana::make_tuple(args...);

    // make space
    IndexSpace::Attr space_attr;
    long count = 0;
    if (components.size() <= 1) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add_intersection: must have at least two "
          "components");
    }
    for (auto&& component : components) {
      const IndexSpace* component_ptr;
      if constexpr (std::is_same_v<std::decay_t<IndexSpaceOrLabel>,
                                   IndexSpace>) {
        component_ptr = &component;
      } else {
        component_ptr = &(this->retrieve(component));
      }
      if (count == 0)
        space_attr = component_ptr->attr();
      else
        space_attr = space_attr.intersection(component_ptr->attr());
      ++count;
    }
    const auto approximate_size = compute_approximate_size(space_attr);

    IndexSpace space(std::forward<S>(type_label), space_attr.type(),
                     space_attr.qns(), approximate_size);
    this->add(space);
    auto type = space.type();

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// @}

  /// @brief removes an IndexSpace associated with `IS.base_key()` from this
  /// @param IS an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& remove(const IndexSpace& IS) {
    auto it = spaces_->find(IS.base_key());
    if (it != spaces_->end()) {
      spaces_->erase(IS);
    }
    return clear_memoized_data_and_return_this();
  }

  /// @brief equivalent to `remove(this->retrieve(label))`
  /// @param label space label
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
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

  /// @brief clear the contents of *this
  /// @return reference to `this`
  IndexSpaceRegistry& clear() {
    *this = IndexSpaceRegistry{};
    return *this;
  }

  /// @brief queries if the intersection space is registered
  /// @param space1
  /// @param space2
  /// @return true if `space1.intersection(space2)` is registered
  bool valid_intersection(const IndexSpace& space1,
                          const IndexSpace& space2) const {
    auto result_attr = space1.attr().intersection(space2.attr());
    return retrieve_ptr(result_attr);
  }

  /// @brief return the resulting space corresponding to a bitwise intersection
  /// between two spaces.
  /// @param space1 a registered IndexSpace
  /// @param space2 a registered IndexSpace
  /// @return the intersection of @p space1 and @p space2
  /// @note can return nullspace
  /// @note throw invalid_argument if the nonnull intersection is not registered
  const IndexSpace& intersection(const IndexSpace& space1,
                                 const IndexSpace& space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      const auto target_qns = space1.qns().intersection(space2.qns());
      bool same_qns = space1.qns() == space2.qns();
      if (!target_qns && !same_qns) {  // spaces with different quantum numbers
                                       // do not intersect.
        return IndexSpace::null;
      }

      // check the registry
      auto intersection_attr = space1.type().intersection(space2.type());
      const IndexSpace& intersection_space =
          find_by_attr({intersection_attr, space1.qns()});
      // the nullspace is a reasonable return value for intersection
      if (intersection_space == IndexSpace::null && intersection_attr) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return intersection_space;
      }
    }
  }

  /// @param space1_key base key of a registered IndexSpace
  /// @param space2_key base key of a registered IndexSpace
  /// @return the intersection of @p space1 and @p space2
  /// @note can return nullspace
  /// @note throw invalid_argument if the nonnull intersection is not registered
  template <typename S1, typename S2,
            typename = meta::EnableIfAllBasicStringConvertible<S1, S2>>
  const IndexSpace& intersection(S1&& space1_key, S2&& space2_key) const {
    if (!contains(space1_key) || !contains(space2_key))
      throw std::invalid_argument(
          "IndexSpaceRegistry::intersection(s1,s2): s1 and s2 must both be "
          "registered");
    return this->intersection(this->retrieve(std::forward<S1>(space1_key)),
                              this->retrieve(std::forward<S2>(space2_key)));
  }

  /// @brief is a union between spaces eligible and registered
  /// @param space1
  /// @param space2
  /// @return true if space is constructable and registered
  bool valid_unIon(const IndexSpace& space1, const IndexSpace& space2) const {
    // check typeattr
    if (!space1.type().includes(space2.type()) &&
        space1.qns() == space2.qns()) {
      // union possible
      auto union_type = space1.type().unIon(space2.type());
      IndexSpace::Attr union_attr{union_type, space1.qns()};
      if (!find_by_attr(union_attr)) {  // possible but not registered
        return false;
      } else
        return true;
    }
    // check qn
    else if (!space1.qns().includes(space2.qns()) &&
             space1.type() == space2.type()) {
      // union possible
      auto union_qn = space1.qns().unIon(space2.qns());
      IndexSpace::Attr union_attr{space1.type(), union_qn};
      if (!find_by_attr(union_attr)) {  // possible but not registered
        return false;
      } else
        return true;
    } else {  // union not mathematically allowed.
      return false;
    }
  }

  /// @param space1
  /// @param space2
  /// @return the union of two spaces.
  /// @note can only return registered spaces
  /// @note never returns nullspace
  const IndexSpace& unIon(const IndexSpace& space1,
                          const IndexSpace& space2) const {
    if (!contains(space1) || !contains(space2))
      throw std::invalid_argument(
          "IndexSpaceRegistry::unIon(s1,s2): s1 and s2 must both be "
          "registered");

    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {
        throw std::invalid_argument(
            "IndexSpaceRegistry::unIon(s1,s2): s1 and s2 must have identical "
            "quantum number attributes.");
      }
      auto unIontype = space1.type().unIon(space2.type());
      const IndexSpace& unIonSpace = find_by_attr({unIontype, space1.qns()});
      if (unIonSpace == IndexSpace::null) {
        throw std::invalid_argument(
            "IndexSpaceRegistry::unIon(s1,s2): the result is not registered, "
            "must register first.");
      } else {
        return unIonSpace;
      }
    }
  }

  /// @param space1_key base key of a registered IndexSpace
  /// @param space2_key base key of a registered IndexSpace
  /// @return the union of two spaces.
  /// @note can only return registered spaces
  /// @note never returns nullspace
  template <typename S1, typename S2,
            typename = meta::EnableIfAllBasicStringConvertible<S1, S2>>
  const IndexSpace& unIon(S1&& space1_key, S2&& space2_key) const {
    if (!contains(space1_key) || !contains(space2_key))
      throw std::invalid_argument(
          "IndexSpaceRegistry::unIon(s1,s2): s1 and s2 must both be "
          "registered");
    return this->unIon(this->retrieve(std::forward<S1>(space1_key)),
                       this->retrieve(std::forward<S2>(space2_key)));
  }

  /// @name physical particle space structure introspection
  /// @{

  /// @brief sets the mask of attributes of `IndexSpace`s  that can be occupied
  /// by physical particles.

  /// Some states do not correspond to physical particles, but may be present
  /// in the registry. Such spaces are not considered when e.g. determining
  /// the base spaces.
  /// @param m the mask of attributes of `IndexSpace`s  that can be occupied by
  /// physical particles.
  void physical_particle_attribute_mask(bitset_t m);

  /// @brief accesses the mask of attributes of `IndexSpace`s  that can be
  /// occupied by physical particles.

  /// Some states do not correspond to physical particles, but may be present
  /// in the registry. Such spaces are not considered when e.g. determining
  /// the base spaces.
  /// @return the mask of attributes of `IndexSpace`s  that can be occupied by
  /// physical particles.
  bitset_t physical_particle_attribute_mask() const;

  /// @brief returns the list of _basis_ IndexSpace::Type objects

  /// A base IndexSpace::Type object has 1 bit in its bitstring.
  /// @sa IndexSpaceRegistry::is_base
  /// @return (memoized) set of base IndexSpace::Type objects, sorted in
  /// increasing order
  const std::vector<IndexSpace::Type>& base_space_types() const {
    if (!base_space_types_) {
      auto types = *spaces_ | ranges::views::transform([](const auto& s) {
        return s.type();
      }) | ranges::views::filter([](const auto& t) { return is_base(t); }) |
                   ranges::views::unique | ranges::to_vector;
      ranges::sort(types, [](auto t1, auto t2) { return t1 < t2; });
      std::scoped_lock guard{mtx_memoized_};
      if (!base_space_types_) {
        base_space_types_ =
            std::make_shared<std::vector<IndexSpace::Type>>(std::move(types));
      }
    }
    return *base_space_types_;
  }

  /// @brief returns the list of _basis_ IndexSpace objects

  /// A base IndexSpace object has 1 bit in its type() bitstring.
  /// @sa IndexSpaceRegistry::is_base
  /// @return (memoized) set of base IndexSpace objects, sorted in the order of
  /// increasing type()
  const std::vector<IndexSpace>& base_spaces() const {
    if (!base_spaces_) {
      auto spaces =
          *spaces_ |
          ranges::views::filter([this](const auto& s) { return is_base(s); }) |
          ranges::views::unique | ranges::to_vector;
      ranges::sort(spaces,
                   [](auto s1, auto s2) { return s1.type() < s2.type(); });
      std::scoped_lock guard{mtx_memoized_};
      if (!base_spaces_) {
        base_spaces_ =
            std::make_shared<std::vector<IndexSpace>>(std::move(spaces));
      }
    }
    return *base_spaces_;
  }

  /// @brief checks if an IndexSpace is in the basis
  /// @param IS IndexSpace
  /// @return true if @p IS is in the basis
  /// @sa base_spaces
  bool is_base(const IndexSpace& IS) const {
    // is base if has base type and has no bits outsize of the physical particle
    // attribute mask
    return is_base(IS.type()) &&
           ((bitset_t(IS.qns()) & (~(physical_particle_attribute_mask()))) ==
            0);
  }

  /// @brief checks if an IndexSpace::Type is in the basis
  /// @param t IndexSpace::Type
  /// @return true if @p t is in the basis
  /// @sa space_type_basis
  static bool is_base(const IndexSpace::Type& t) {
    return has_single_bit(t.to_int32());
  }

  /// @}

  /// @brief an @c IndexSpace is occupied with respect to the fermi vacuum or a
  /// subset of that space
  /// @note only makes sense to ask this if in a SingleProduct vacuum context.
  bool is_pure_occupied(const IndexSpace& IS) const {
    if (!IS) {
      return false;
    }
    if (IS.type().to_int32() <=
        vacuum_occupied_space(IS.qns()).type().to_int32()) {
      return true;
    } else {
      return false;
    }
  }

  /// @brief all states are unoccupied in the fermi vacuum
  /// @note again, this only makes sense to ask if in a SingleProduct vacuum
  /// context.
  bool is_pure_unoccupied(const IndexSpace& IS) const {
    if (!IS) {
      return false;
    } else {
      return !IS.type().intersection(vacuum_occupied_space(IS.qns()).type());
    }
  }

  /// @brief some states are fermi vacuum occupied
  bool contains_occupied(const IndexSpace& IS) const {
    return IS.type().intersection(vacuum_occupied_space(IS.qns()).type()) !=
           IndexSpace::Type::null;
  }

  /// @brief some states are fermi vacuum unoccupied
  bool contains_unoccupied(const IndexSpace& IS) const {
    return IS.type().intersection(vacuum_unoccupied_space(IS.qns()).type()) !=
           IndexSpace::Type::null;
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
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
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
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpaceRegistry& vacuum_occupied_space(S&& l) {
    return vacuum_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in vacuum state for any set of quantum numbers
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// vacuum_occupied_space had not been specified
  const IndexSpace::Type& vacuum_occupied_space(
      bool nulltype_ok = false) const {
    if (!std::get<0>(vacocc_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "vacuum occupied space has not been specified, invoke "
          "vacuum_occupied_space(IndexSpace::Type) or "
          "vacuum_occupied_space(container::map<IndexSpace::QuantumNumbers,"
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
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
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
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpaceRegistry& reference_occupied_space(S&& l) {
    return reference_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in reference state for any set of quantum
  /// numbers
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// reference_occupied_space had not been specified
  const IndexSpace::Type& reference_occupied_space(
      bool nulltype_ok = false) const {
    if (!std::get<0>(refocc_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "reference occupied space has not been specified, invoke "
          "reference_occupied_space(IndexSpace::Type) or "
          "reference_occupied_space(container::map<IndexSpace::QuantumNumbers,"
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
      return retrieve(this->reference_occupied_space(), qn);
    }
  }

  /// @}

  /// @name  specifies which spaces comprise the entirety of Hilbert space
  /// @note needed for creating general operators in mbpt/op
  /// @{

  /// @param s an IndexSpace::Type specifying the complete Hilbert space;
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
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
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
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpaceRegistry& complete_space(S&& l) {
    return complete_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the complete Hilbert space for any set of quantum numbers
  /// @throw std::invalid_argument if @p nulltype_ok is false and complete_space
  /// had not been specified
  const IndexSpace::Type& complete_space(bool nulltype_ok = false) const {
    if (!std::get<0>(complete_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
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

  /// @return the space that is unoccupied in the vacuum state
  const IndexSpace& vacuum_unoccupied_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto complete_type = this->complete_space(qn).type();
    auto vacocc_type = this->vacuum_occupied_space(qn).type();
    auto vacuocc_type =
        complete_type.xOr(vacocc_type).intersection(complete_type);
    return this->retrieve(vacuocc_type, qn);
  }

  /// @name specifies in which space holes can be created successfully from the
  /// reference wave function
  /// @note convenience for making operators
  /// @{

  /// @param t an IndexSpace::Type specifying where holes can be created;
  ///          to specify hole space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& hole_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "hole_space");
    std::get<0>(hole_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which holes
  /// can be created
  IndexSpaceRegistry& hole_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "hole_space");
    std::get<1>(hole_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `hole_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& hole_space(const IndexSpace& s) {
    return hole_space(s.type());
  }

  /// equivalent to `hole_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpaceRegistry& hole_space(S&& l) {
    return hole_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which holes can be created
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// hole_space had not been specified
  const IndexSpace::Type& hole_space(bool nulltype_ok = false) const {
    if (!std::get<0>(hole_space_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "active hole space has not been specified, invoke "
          "hole_space(IndexSpace::Type) or "
          "hole_space(container::map<IndexSpace::QuantumNumbers,IndexSpace::"
          "Type>)");
    } else
      return std::get<0>(hole_space_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space in which holes can be created for the given set of
  /// quantum numbers
  const IndexSpace& hole_space(const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(hole_space_).find(qn);
    if (it != std::get<1>(hole_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->hole_space(), qn);
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
  IndexSpaceRegistry& particle_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "particle_space");
    std::get<0>(particle_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which
  /// particles can be created
  IndexSpaceRegistry& particle_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "particle_space");
    std::get<1>(particle_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `particle_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& particle_space(const IndexSpace& s) {
    return particle_space(s.type());
  }

  /// equivalent to `particle_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAllBasicStringConvertible<S>>
  IndexSpaceRegistry& particle_space(S&& l) {
    return particle_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which particles can be created
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// particle_space had not been specified
  const IndexSpace::Type& particle_space(bool nulltype_ok = false) const {
    if (!std::get<0>(particle_space_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "active particle space has not been specified, invoke "
          "particle_space(IndexSpace::Type) or "
          "particle_space(container::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(particle_space_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space in which particles can be created for the given set of
  /// quantum numbers
  const IndexSpace& particle_space(const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(particle_space_).find(qn);
    if (it != std::get<1>(particle_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->particle_space(), qn);
    }
  }

  /// @}

  /// @}

 private:
  // N.B. need transparent comparator, see https://stackoverflow.com/a/35525806
  std::shared_ptr<container::set<IndexSpace, IndexSpace::KeyCompare>> spaces_;

  bitset_t physical_particle_attribute_mask_ = bitset::null;

  // memoized data
  mutable std::shared_ptr<std::vector<IndexSpace::Type>> base_space_types_;
  mutable std::shared_ptr<std::vector<IndexSpace>> base_spaces_;
  mutable std::recursive_mutex
      mtx_memoized_;  // used to update the memoized data
  IndexSpaceRegistry& clear_memoized_data_and_return_this() {
    std::scoped_lock guard{mtx_memoized_};
    base_space_types_.reset();
    base_spaces_.reset();
    return *this;
  }

  ///@brief true if has one and only one bit set.
  static bool has_single_bit(std::uint32_t bits) {
    return bits && !(bits & (bits - 1));
  }

  /// @brief find an IndexSpace from its attr. return nullspace if not present.
  /// @param attr the attribute of the IndexSpace
  const IndexSpace& find_by_attr(const IndexSpace::Attr& attr) const {
    for (auto&& space : *spaces_) {
      if (space.attr() == attr) {
        return space;
      }
    }
    return IndexSpace::null;
  }

  void throw_if_missing(const IndexSpace::Type& t,
                        const IndexSpace::QuantumNumbers& qn,
                        std::string call_context = "") {
    for (auto&& space : *spaces_) {
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
    for (auto&& space : *spaces_) {
      if (space.type() == t) {
        return;
      }
    }
    throw std::invalid_argument(call_context + ": missing { IndexSpace::Type=" +
                                std::to_string(t.to_int32()) +
                                " , any IndexSpace::QuantumNumbers } space");
  }

  void throw_if_missing_any(const container::map<IndexSpace::QuantumNumbers,
                                                 IndexSpace::Type>& qn2type,
                            std::string call_context = "") {
    container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type_found;
    for (auto&& space : *spaces_) {
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
#if __cplusplus < 202002L
      if (qn2type_found.find(qn) == qn2type_found.end()) {
#else
      if (!qn2type_found.contains(qn)) {
#endif
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

  // used for fermi vacuum wick application
  std::tuple<IndexSpace::Type,
             container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      vacocc_ = {{}, {}};

  // used for MR MBPT to take average over multiconfiguration reference
  std::tuple<IndexSpace::Type,
             container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      refocc_ = {{}, {}};

  // defines active bits in TypeAttr; used by general operators in mbpt/op
  std::tuple<IndexSpace::Type,
             container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      complete_ = {{}, {}};

  // both needed to make excitation and de-excitation operators. not
  // necessarily equivalent in the case of multi-reference context.
  std::tuple<IndexSpace::Type,
             container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      hole_space_ = {{}, {}};
  std::tuple<IndexSpace::Type,
             container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      particle_space_ = {{}, {}};

  // Boost.Hana snippet to process attribute tag arguments
  template <typename ArgsHanaTuple>
  void process_attribute_tags(ArgsHanaTuple h_tuple,
                              const IndexSpace::Type& type) {
    boost::hana::for_each(h_tuple, [this, &type](auto arg) {
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
        this->hole_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsParticle>) {
        this->particle_space(type);
      } else {
        static_assert(meta::always_false<decltype(arg)>::value,
                      "IndexSpaceRegistry::add{,_union,_intersect}: unknown "
                      "attribute tag");
      }
    });
  }

  /// @brief computes the approximate size of the space

  /// for a base space return its extent, for a composite space compute as a sum
  /// of extents of base subspaces
  /// @param space_attr the IndexSpace attribute
  /// @return the approximate size of the space
  unsigned long compute_approximate_size(
      const IndexSpace::Attr& space_attr) const {
    if (is_base(space_attr.type())) {
      return this->retrieve(space_attr).approximate_size();
    } else {
      // compute_approximate_size is used when populating the registry
      // so don't use base_spaces() here
      unsigned long size = ranges::accumulate(
          *spaces_ | ranges::views::filter([this, &space_attr](auto& s) {
            return s.qns() == space_attr.qns() && is_base(s.type()) &&
                   space_attr.type().intersection(s.type());
          }),
          0ul, [](unsigned long size, const IndexSpace& s) {
            return size + s.approximate_size();
          });
      return size;
    }
  }

  friend bool operator==(const IndexSpaceRegistry& isr1,
                         const IndexSpaceRegistry& isr2) {
    return *isr1.spaces_ == *isr2.spaces_;
  }
};  // class IndexSpaceRegistry

}  // namespace sequant
#endif  // SEQUANT_INDEX_SPACE_REGISTRY_HPP
