//
// Created by Conner Masteran on 4/16/24.
//

#ifndef SEQUANT_INDEX_BASIS_REGISTRY_HPP
#define SEQUANT_INDEX_BASIS_REGISTRY_HPP

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/bitset.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/unique.hpp>

#include <boost/hana.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <mutex>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>

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
/// @brief set of known IndexSpace objects and named basis instances

/// Each IndexSpace object has hardwired base key (label) that gives
/// indexed expressions appropriate semantics; e.g., spaces referred to by
/// indices in \f$ t_{p_1}^{i_1} \f$ are defined if IndexSpace objects with
/// base keys \f$ p \f$ and \f$ i \f$ are registered.
/// Since index spaces have set-theoretic semantics, the user must
/// provide complete set of unions/intersects of the base spaces to
/// cover all possible IndexSpace objects that can be generated in their
/// program.
///
/// Registry contains 2 parts: a table of IndexBasis objects keyed by label
/// (IndexBasisRegistry::bases(); a registered IndexSpace is the entry under its
/// base_key() with a null basis instance, see IndexBasisRegistry::spaces(); a
/// *named basis instance* is an entry under its own label, printed and parsed
/// by that label), and the specification of various spaces (vacuum, reference,
/// complete, etc.). Registries can be constructed from the same table but have
/// different specifications of vacuum, reference, etc.; this is useful for
/// providing different contexts for fermions and bosons, for example.
///
/// Spaces that can be occupied by physical particles need to be
/// introspected for their structure, occupancy, etc. The registry provides the
/// API needed for dealing with such states.
/// - IndexBasisRegistry::physical_particle_attribute_mask specify which states are occupied by
/// physical particles
/// - every space occupied by physical particles is a union of basis
/// ("base") spaces. IndexBasisRegistry::is_base detects such spaces and
/// IndexBasisRegistry::base_space_types/IndexBasisRegistry::base_spaces report the list of base spaces
/// - in SingleProduct vacuum IndexBasisRegistry::is_pure_occupied/IndexBasisRegistry::is_pure_unoccupied report
/// whether a space is occupied or unoccupied (always the case for base
/// spaces, but neither may be true for composite  spaces)
/// - IndexBasisRegistry::vacuum_occupied_space report whether a space has nonzero occupancy
/// in the vacuum state (that defines the normal order); this is needed to
/// apply Wick theorem with Fermi vacuum
/// - IndexBasisRegistry::reference_occupied_space reports whether a space has nonzero occupancy
/// in the reference state used to compute reference
/// expectation value; only needed for computing expectation values when the
/// vacuum state does not match the expectation value state.
/// - IndexBasisRegistry::complete_space specifies which spaces comprise the entirety of Hilbert
/// space; needed for creating general operators in mbpt/op
/// - IndexBasisRegistry::particle_space and IndexBasisRegistry::hole_space specify in which space particles/holes
/// can be created successfully from the reference state; this is a
/// convenience for making operators
// clang-format on
class IndexBasisRegistry {
 public:
  /// label -> basis: a space is the entry under its base_key() with a null
  /// instance, a named basis instance any other entry
  using table_type = container::map<std::wstring, IndexBasis, std::less<>>;

  /// exception type thrown when a label names a basis instance where a space
  /// is required
  struct not_a_space : IndexSpace::bad_key {
    using bad_key::bad_key;
  };

  /// the space entries of a registry, as IndexSpace objects, in label order;
  /// valid as long as the registry is neither modified nor destroyed
  class SpacesView {
   public:
    class iterator {
     public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = IndexSpace;
      using difference_type = std::ptrdiff_t;
      using pointer = const IndexSpace*;
      using reference = const IndexSpace&;

      iterator() = default;
      iterator(table_type::const_iterator it, table_type::const_iterator end)
          : it_(it), end_(end) {
        skip();
      }

      reference operator*() const { return it_->second.space(); }
      pointer operator->() const { return &it_->second.space(); }
      iterator& operator++() {
        ++it_;
        skip();
        return *this;
      }
      iterator operator++(int) {
        auto old = *this;
        ++*this;
        return old;
      }
      friend bool operator==(const iterator&, const iterator&) = default;

     private:
      void skip() {
        while (it_ != end_ && it_->second.has_basis_instance()) ++it_;
      }
      table_type::const_iterator it_{}, end_{};
    };

    explicit SpacesView(const table_type& table) : table_(&table) {}
    iterator begin() const { return {table_->cbegin(), table_->cend()}; }
    iterator end() const { return {table_->cend(), table_->cend()}; }

   private:
    const table_type* table_;
  };

  /// default constructor creates a registry containing only IndexSpace::null
  /// @note null space is registered so we don't have to handle it as a corner
  /// case in retrieve() and other methods
  IndexBasisRegistry() {
    // register nullspace
    this->add(IndexSpace::null);
  }

  /// constructs an IndexBasisRegistry from an existing table (e.g. the
  /// bases() of another registry), spaces and named basis instances alike
  /// @note the table is taken as given, it is not validated against the
  /// invariants that add() enforces, except that each label is a valid base
  /// key (the null space's empty key aside) and each named basis instance
  /// carries the label it is registered under
  /// @throw Exception if a label is not a valid base key (see
  /// io::serialization::v1::is_base_key())
  explicit IndexBasisRegistry(table_type bases)
      : bases_(std::move(bases)),
        named_count_(ranges::count_if(
            bases_, [](const auto& e) { return !is_space(e); })) {
    for (auto& [label, basis] : bases_) {
      if (basis.space()) validate_label(label, "IndexBasisRegistry(table)");
      if (basis.has_basis_instance() && basis.name() != label)
        basis = IndexBasis(basis.space(), basis.basis_instance(), label,
                           basis.extent_, basis.metric_, basis.field_);
    }
  }

  /// copy constructor
  IndexBasisRegistry(const IndexBasisRegistry& other)
      : bases_(other.bases_),
        named_count_(other.named_count_),
        physical_particle_attribute_mask_(
            other.physical_particle_attribute_mask_),
        vacocc_(other.vacocc_),
        refocc_(other.refocc_),
        complete_(other.complete_),
        hole_space_(other.hole_space_),
        particle_space_(other.particle_space_) {}

  /// move constructor
  IndexBasisRegistry(IndexBasisRegistry&& other)
      : bases_(std::move(other.bases_)),
        named_count_(other.named_count_),
        physical_particle_attribute_mask_(
            std::move(other.physical_particle_attribute_mask_)),
        vacocc_(std::move(other.vacocc_)),
        refocc_(std::move(other.refocc_)),
        complete_(std::move(other.complete_)),
        hole_space_(std::move(other.hole_space_)),
        particle_space_(std::move(other.particle_space_)) {
    other.named_count_ = 0;
    // what other has memoized describes the spaces it gave up
    other.clear_memoized_data_and_return_this();
  }

  /// copy assignment operator
  IndexBasisRegistry& operator=(const IndexBasisRegistry& other) {
    bases_ = other.bases_;
    named_count_ = other.named_count_;
    physical_particle_attribute_mask_ = other.physical_particle_attribute_mask_;
    vacocc_ = other.vacocc_;
    refocc_ = other.refocc_;
    complete_ = other.complete_;
    hole_space_ = other.hole_space_;
    particle_space_ = other.particle_space_;
    return clear_memoized_data_and_return_this();
  }

  /// move assignment operator
  IndexBasisRegistry& operator=(IndexBasisRegistry&& other) {
    if (this == &other) return *this;
    bases_ = std::move(other.bases_);
    named_count_ = other.named_count_;
    physical_particle_attribute_mask_ =
        std::move(other.physical_particle_attribute_mask_);
    vacocc_ = std::move(other.vacocc_);
    refocc_ = std::move(other.refocc_);
    complete_ = std::move(other.complete_);
    hole_space_ = std::move(other.hole_space_);
    particle_space_ = std::move(other.particle_space_);
    other.named_count_ = 0;
    // what other has memoized describes the spaces it gave up
    other.clear_memoized_data_and_return_this();
    return clear_memoized_data_and_return_this();
  }

  /// @return view of the registered spaces (the entries with a null basis
  /// instance)
  SpacesView spaces() const { return SpacesView{bases_}; }

  SpacesView::iterator begin() const { return spaces().begin(); }
  SpacesView::iterator end() const { return spaces().end(); }

  /// @return the table of all entries, spaces and named basis instances
  const table_type& bases() const { return bases_; }

  /// @brief retrieve a pointer to the IndexBasis registered under a label
  /// @param label a label of a space or of a named basis instance, or a label
  /// of an Index (see Index::label() )
  /// @return pointer to the IndexBasis registered under that key, or nullptr
  /// if not found
  template <basic_string_convertible S>
  const IndexBasis* retrieve_basis_ptr(S&& label) const {
    auto it = bases_.find(IndexSpace::reduce_key(to_basic_string_view(label)));
    return it != bases_.end() ? &it->second : nullptr;
  }

  /// @brief retrieve the IndexBasis registered under a label
  /// @param label a label of a space or of a named basis instance, or a label
  /// of an Index (see Index::label() )
  /// @return the IndexBasis registered under that key
  /// @throw IndexSpace::bad_key if no entry is registered under that key
  template <basic_string_convertible S>
  const IndexBasis& retrieve_basis(S&& label) const {
    if (const auto* b = retrieve_basis_ptr(label)) return *b;
    throw IndexSpace::bad_key(label);
  }

  /// @return the label under which the space and instance of @p b are
  /// registered (see same_instance()), or std::nullopt if @p b has no basis
  /// instance or they are not named
  std::optional<std::wstring_view> basis_label(const IndexBasis& b) const {
    if (!b.has_basis_instance() || named_count_ == 0) return std::nullopt;
    for (auto const& [label, basis] : bases_)
      if (basis.has_basis_instance() && same_instance(basis, b))
        return std::wstring_view(label);
    return std::nullopt;
  }

  /// @return the named entry of the space and instance of @p b (see
  /// same_instance()), which carries the entry's name, extent, metric and
  /// field, if they are registered under a name; otherwise @p b unchanged
  IndexBasis resolve(const IndexBasis& b) const {
    if (!b.has_basis_instance() || named_count_ == 0) return b;
    for (auto const& [label, basis] : bases_)
      if (same_instance(basis, b)) return basis;
    return b;
  }

  /// @brief sets the extent of the entry registered under a label: the
  /// dimension of a space (the extent of its own basis, and of its named
  /// basis instances registered without an extent of their own), or the
  /// extent of a named basis instance
  /// @param label a label of a space or of a named basis instance
  /// @param n the extent
  /// @return reference to `this`
  /// @throw IndexSpace::bad_key if no entry is registered under @p label
  template <basic_string_convertible S>
  IndexBasisRegistry& extent(S&& label, std::size_t n) {
    auto& entry = entry_or_throw(label);
    if (entry.has_basis_instance())
      entry.extent(n);
    else
      for (IndexSpace& space : space_copies_of(entry)) space.dimension(n);
    return clear_memoized_data_and_return_this();
  }

  /// @deprecated use extent(label, n)
  template <basic_string_convertible S>
  [[deprecated("use extent(label, n)")]] IndexBasisRegistry& approximate_size(
      S&& label, std::size_t n) {
    return extent(std::forward<S>(label), n);
  }

  /// @brief sets the Field of the entry registered under a label: of a
  /// space's own basis (and of its named basis instances registered without
  /// a field of their own), or of a named basis instance
  /// @param label a label of a space or of a named basis instance
  /// @param f the Field
  /// @return reference to `this`
  /// @throw IndexSpace::bad_key if no entry is registered under @p label
  template <basic_string_convertible S>
  IndexBasisRegistry& field(S&& label, Field f) {
    auto& entry = entry_or_throw(label);
    if (entry.has_basis_instance())
      entry.field(f);
    else
      for (IndexSpace& space : space_copies_of(entry)) space.field(f);
    return clear_memoized_data_and_return_this();
  }

  /// @brief sets the metric of the named basis instance registered under a
  /// label
  /// @param label a label of a named basis instance
  /// @param m the metric
  /// @return reference to `this`
  /// @throw IndexSpace::bad_key if no entry is registered under @p label
  /// @throw Exception if @p label is that of a space: the own basis of a
  /// space is orthonormal
  template <basic_string_convertible S>
  IndexBasisRegistry& metric(S&& label, IndexSpaceMetric m) {
    auto& entry = entry_or_throw(label);
    if (!entry.has_basis_instance())
      throw Exception("IndexBasisRegistry::metric: '" + toUtf8(label) +
                      "' is a space, whose own basis is orthonormal; register "
                      "a non-orthonormal basis of it under a name");
    entry.metric(m);
    return clear_memoized_data_and_return_this();
  }

  /// @brief retrieve a pointer to IndexSpace from the registry by the label
  /// @param label a @c base_key of an IndexSpace, or a label of an Index (see
  /// Index::label() )
  /// @return pointer to IndexSpace associated with that key, or nullptr if not
  /// found or if the key names a basis instance
  /// @note a space's metadata is written with extent(label, n) and
  /// field(label, f), which also reach its named basis instances
  template <basic_string_convertible S>
  const IndexSpace* retrieve_ptr(S&& label) const {
    const auto* b = retrieve_basis_ptr(std::forward<S>(label));
    return b && !b->has_basis_instance() ? &b->space() : nullptr;
  }

  /// @brief retrieve an IndexSpace from the registry by the label
  /// @param label a @c base_key of an IndexSpace, or a label of an Index (see
  /// Index::label() )
  /// @return IndexSpace associated with that key
  /// @throw not_a_space if the key names a basis instance
  /// @throw IndexSpace::bad_key if matching space is not found
  template <basic_string_convertible S>
  const IndexSpace& retrieve(S&& label) const {
    if (const auto* ptr = retrieve_ptr(label)) {
      return *ptr;
    } else if (retrieve_basis_ptr(label)) {
      throw not_a_space(label);
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
    auto it = std::find_if(bases_.begin(), bases_.end(), [&](const auto& e) {
      return is_space(e) && e.second.space().type() == type &&
             e.second.space().qns() == qns;
    });
    return it != bases_.end() ? &it->second.space() : nullptr;
  }

  /// @brief retrieve an IndexSpace from the registry by its type and quantum
  /// numbers
  /// @param type IndexSpace::Type
  /// @param qns IndexSpace::QuantumNumbers
  /// @return IndexSpace associated with that key.
  /// @throw Exception if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Type& type,
                             const IndexSpace::QuantumNumbers& qns) const {
    if (const auto* ptr = retrieve_ptr(type, qns)) {
      return *ptr;
    } else
      throw Exception(
          "IndexBasisRegistry::retrieve(type,qn): missing { IndexSpace::Type=" +
          to_string(type) + " , IndexSpace::QuantumNumbers=" + to_string(qns) +
          " } combination");
  }

  /// @brief retrieve pointer to the IndexSpace from the registry by the
  /// IndexSpace::Attr
  /// @param space_attr an IndexSpace::Attr
  /// @return pointer to the IndexSpace associated with that key, or nullptr if
  /// not found
  const IndexSpace* retrieve_ptr(const IndexSpace::Attr& space_attr) const {
    auto it = std::find_if(bases_.begin(), bases_.end(), [&](const auto& e) {
      return is_space(e) && e.second.space().attr() == space_attr;
    });
    return it != bases_.end() ? &it->second.space() : nullptr;
  }

  /// @brief retrieve an IndexSpace from the registry by the IndexSpace::Attr
  /// @param space_attr an IndexSpace::Attr
  /// @return IndexSpace associated with that key.
  /// @throw Exception if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Attr& space_attr) const {
    if (const auto* ptr = retrieve_ptr(space_attr)) {
      return *ptr;
    } else
      throw Exception(std::string("IndexBasisRegistry::retrieve(attr): attr=") +
                      to_string(space_attr) + " is missing");
  }

  /// queries presence of a registered IndexSpace or named basis instance
  /// @param label a @c base_key of an IndexSpace, a label of a named basis
  /// instance, or a label of an Index (see Index::label() )
  /// @return true, if an entry with key @p label is registered
  template <basic_string_convertible S>
  bool contains(S&& label) const {
    return retrieve_basis_ptr(std::forward<S>(label)) != nullptr;
  }

  /// queries presence of a registered basis
  /// @param b an IndexBasis object
  /// @return true, if @p b has a basis instance and is named, or has none and
  /// its space is registered
  bool contains(const IndexBasis& b) const {
    return b.has_basis_instance() ? basis_label(b).has_value()
                                  : contains(b.space());
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
  /// @throw Exception if `IS.base_key()` is not a valid base key
  /// (see io::serialization::v1::is_base_key()), or if it or
  /// `IS.attr()` matches an already registered IndexSpace
  IndexBasisRegistry& add(const IndexSpace& IS) {
    if (IS)
      validate_label(IS.base_key(), "IndexBasisRegistry::add(index_space)");
    auto it = bases_.find(IS.base_key());
    if (it != bases_.end()) {
      throw Exception(
          (std::string("IndexBasisRegistry::add(index_space): an entry with "
                       "index_space.base_key()=") +
           toUtf8(IS.base_key()) +
           " already in the registry; if you are trying to replace the "
           "IndexSpace use "
           "IndexBasisRegistry::replace(is)"));
    } else {
      // make sure there are no duplicate IndexSpaces whose attribute is
      // IS.attr()
      if (ranges::any_of(bases_, [&IS](const auto& e) {
            return is_space(e) && e.second.space().attr() == IS.attr();
          })) {
        throw Exception(
            (std::string("IndexBasisRegistry::add(index_space): space with "
                         "index_space.attr()=") +
             to_string(IS.attr()) +
             " already in the registry; if you are trying to replace the "
             "IndexSpace use "
             "IndexBasisRegistry::replace(is)"));
      }
      bases_.emplace(IS.base_key(), IndexBasis{IS});
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
  ///   - dimension of the space (unsigned long)
  ///   - any of { is_vacuum_occupied , is_reference_occupied , is_complete ,
  ///   is_hole , is_particle }
  /// @return reference to `this`
  /// @throw Exception if `type_label` or `type` matches
  /// an already registered IndexSpace
  template <basic_string_convertible S, typename... OptionalArgs>
  IndexBasisRegistry& add(S&& type_label, IndexSpace::Type type,
                          OptionalArgs&&... args) {
    auto h_args = boost::hana::make_tuple(args...);

    // process IndexSpace::QuantumNumbers, set to default if not given
    auto h_qns = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::type_c<decltype(arg)> ==
             boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    constexpr auto nqns = boost::hana::size(h_qns);
    static_assert(
        nqns == boost::hana::size_c<0> || nqns == boost::hana::size_c<1>,
        "IndexBasisRegistry::add: only one IndexSpace::QuantumNumbers argument "
        "is allowed");
    constexpr auto have_qns = nqns == boost::hana::size_c<1>;
    IndexSpace::QuantumNumbers qns;
    if constexpr (have_qns) {
      qns = boost::hana::at_c<0>(h_qns);
    }

    // process dimension and Field, set to defaults if not given
    const auto [size, field] = parse_size_and_field(args...);

    // make space
    IndexSpace space(std::forward<S>(type_label), type, qns, size.value_or(10),
                     field.value_or(Field::Complex));
    this->add(space);

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers> &&
             boost::hana::type_c<decltype(arg)> != boost::hana::type_c<Field>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// @brief registers a basis instance of a registered space under its own
  /// label
  /// @param label the label of the basis instance; must be a valid base
  /// key (see io::serialization::v1::is_base_key()): letters, ⁺, ⁻,
  /// combining diacritics, arrows and primes, so no digits or `_`, which an
  /// index label reserves for the ordinal
  /// @param basis an IndexBasis with a basis instance whose space is
  /// registered; the entry is named @p label whatever @p basis is named
  /// @param args optional arguments consisting of a mix of zero or one of
  /// each of the following, each defaulting to what @p basis carries (for a
  /// basis built from a space and an instance: the dimension of the space,
  /// IndexSpaceMetric::Unit and the field of the space's own basis):
  ///   - extent of the basis (unsigned long)
  ///   - IndexSpaceMetric
  ///   - Field
  /// @return reference to `this`
  /// @throw Exception if @p label is not a valid label or is already
  /// registered, if @p basis has no basis instance, if its space is not
  /// registered, or if @p basis is already named
  template <basic_string_convertible S, typename... Args>
  IndexBasisRegistry& add(S&& label, const IndexBasis& basis, Args&&... args) {
    auto h_args = boost::hana::make_tuple(args...);
    auto h_tags = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> != boost::hana::type_c<Field> &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpaceMetric>;
    });
    static_assert(boost::hana::size(h_tags) == boost::hana::size_c<0>,
                  "IndexBasisRegistry::add(label, basis): attribute tags are "
                  "per space; only an integral extent, an IndexSpaceMetric "
                  "and a Field may be given for a basis instance");
    const auto [size, field] = parse_size_and_field(args...);
    auto h_metric = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::type_c<decltype(arg)> ==
             boost::hana::type_c<IndexSpaceMetric>;
    });
    constexpr auto nmetrics = boost::hana::size(h_metric);
    static_assert(
        nmetrics == boost::hana::size_c<0> ||
            nmetrics == boost::hana::size_c<1>,
        "IndexBasisRegistry::add(label, basis): only one IndexSpaceMetric "
        "argument is allowed");
    IndexSpaceMetric metric = basis.metric_;
    if constexpr (nmetrics == boost::hana::size_c<1>) {
      metric = boost::hana::at_c<0>(h_metric);
    }
    std::wstring key = toUtf16(std::forward<S>(label));
    if (!basis.has_basis_instance())
      throw Exception("IndexBasisRegistry::add(label, basis): '" + toUtf8(key) +
                      "' needs a basis instance; register a space with "
                      "add(IndexSpace)");
    validate_label(key, "IndexBasisRegistry::add(label, basis)");
    if (bases_.contains(key))
      throw Exception("IndexBasisRegistry::add(label, basis): label '" +
                      toUtf8(key) + "' is already registered");
    const auto* space = retrieve_ptr(basis.space().type(), basis.space().qns());
    if (!space || *space != basis.space())
      throw Exception("IndexBasisRegistry::add(label, basis): the space of '" +
                      toUtf8(key) + "' is not registered");
    if (auto taken = basis_label(basis))
      throw Exception("IndexBasisRegistry::add(label, basis): the basis of '" +
                      toUtf8(key) + "' is already named '" + toUtf8(*taken) +
                      "'");
    IndexBasis named{
        *space, basis.basis_instance(),
        key,    size ? std::optional<std::size_t>(*size) : basis.extent_,
        metric, field ? field : basis.field_};
    bases_.emplace(std::move(key), std::move(named));
    ++named_count_;
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
  template <basic_string_convertible S, index_space_or_label IndexSpaceOrLabel,
            typename... OptionalArgs>
  IndexBasisRegistry& add_unIon(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    SEQUANT_ASSERT(components.size() > 1);

    auto h_args = boost::hana::make_tuple(args...);

    // make space
    IndexSpace::Attr space_attr;
    long count = 0;
    if (components.size() <= 1) {
      throw Exception(
          "IndexBasisRegistry::add_unIon: must have at least two components");
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
    const auto dimension = compute_dimension(space_attr);
    const Field field = compute_field(space_attr);

    IndexSpace space(std::forward<S>(type_label), space_attr.type(),
                     space_attr.qns(), dimension, field);
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
  template <basic_string_convertible S, index_space_or_label IndexSpaceOrLabel,
            typename... OptionalArgs>
  IndexBasisRegistry& add_union(
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
  template <basic_string_convertible S, index_space_or_label IndexSpaceOrLabel,
            typename... OptionalArgs>
  IndexBasisRegistry& add_intersection(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    SEQUANT_ASSERT(components.size() > 1);

    auto h_args = boost::hana::make_tuple(args...);

    // make space
    IndexSpace::Attr space_attr;
    long count = 0;
    if (components.size() <= 1) {
      throw Exception(
          "IndexBasisRegistry::add_intersection: must have at least two "
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
    const auto dimension = compute_dimension(space_attr);

    IndexSpace space(std::forward<S>(type_label), space_attr.type(),
                     space_attr.qns(), dimension);
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
  /// @throw Exception if basis instances of the space are named
  IndexBasisRegistry& remove(const IndexSpace& IS) {
    auto it = bases_.find(IS.base_key());
    if (it != bases_.end() && is_space(*it)) {
      std::string named;
      for (auto const& [label, basis] : bases_)
        if (basis.has_basis_instance() && basis.space() == it->second.space())
          named += (named.empty() ? "" : ", ") + toUtf8(label);
      if (!named.empty())
        throw Exception("IndexBasisRegistry::remove: space " +
                        toUtf8(IS.base_key()) +
                        " still has named basis instances: " + named);
      bases_.erase(it);
    }
    return clear_memoized_data_and_return_this();
  }

  /// @brief removes the named basis instance registered under @p label, else
  /// equivalent to `remove(this->retrieve(label))`
  /// @param label label of a space or of a named basis instance
  /// @return reference to `this`
  template <basic_string_convertible S>
  IndexBasisRegistry& remove(S&& label) {
    auto it = bases_.find(IndexSpace::reduce_key(to_basic_string_view(label)));
    if (it != bases_.end() && !is_space(*it)) {
      bases_.erase(it);
      --named_count_;
      return clear_memoized_data_and_return_this();
    }
    auto&& IS = this->retrieve(std::forward<S>(label));
    return this->remove(IS);
  }

  /// @brief replaces an IndexSpace registered in the registry under
  /// IS.base_key()
  ///        with @p IS
  /// @param IS an IndexSpace
  /// @return reference to `this`
  /// @throw Exception if the replaced space has named basis instances
  IndexBasisRegistry& replace(const IndexSpace& IS) {
    this->remove(IS);
    return this->add(IS);
  }

  /// @brief clear the contents of *this
  /// @return reference to `this`
  IndexBasisRegistry& clear() {
    *this = IndexBasisRegistry{};
    return *this;
  }

  /// @brief queries if the intersection space is nonnull and registered
  /// @param space1
  /// @param space2
  /// @return true if `space1.intersection(space2)` is nonnull and registered
  bool valid_intersection(const IndexSpace& space1,
                          const IndexSpace& space2) const {
    auto result_attr = space1.attr().intersection(space2.attr());
    return result_attr != IndexSpace::Type::null && retrieve_ptr(result_attr);
  }

  /// @brief queries if the intersection space is nonnull and registered
  /// @param space1_key base key of a registered IndexSpace
  /// @param space2_key base key of a registered IndexSpace
  /// @return true if
  /// `valid_intersection(retrieve(space1_key),retrieve(space2_key))` is
  /// nonnull and registered
  template <basic_string_convertible S1, basic_string_convertible S2>
  bool valid_intersection(S1&& space1_key, S2&& space2_key) const {
    const auto& space1 = retrieve_or_throw(
        space1_key,
        "IndexBasisRegistry::valid_intersection(s1,s2): space with key s1=");
    const auto& space2 = retrieve_or_throw(
        space2_key,
        "IndexBasisRegistry::valid_intersection(s1,s2): space key s2=");
    return this->valid_intersection(space1, space2);
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
      auto intersection_type = space1.type().intersection(space2.type());
      const IndexSpace& intersection_space =
          find_by_attr({intersection_type, space1.qns()});
      // the nullspace is a reasonable return value for intersection
      if (intersection_space == IndexSpace::null && intersection_type) {
        throw Exception(
            std::string("intersection(s1=") + to_string(space1) +
            ",s2=" + to_string(space2) +
            ": no space with resulting type=" + to_string(intersection_type) +
            " is found in the registry. Add a "
            "space with this type to the registry.");
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
  template <basic_string_convertible S1, basic_string_convertible S2>
  const IndexSpace& intersection(S1&& space1_key, S2&& space2_key) const {
    const auto& space1 = retrieve_or_throw(
        space1_key,
        "IndexBasisRegistry::intersection(s1,s2): space with key s1=");
    const auto& space2 = retrieve_or_throw(
        space2_key, "IndexBasisRegistry::intersection(s1,s2): space key s2=");
    return this->intersection(space1, space2);
  }

  /// @brief is a union between spaces is registered
  /// @param space1
  /// @param space2
  /// @return true if space is registered
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

  /// @brief is a union between spaces is registered
  /// @param space1_key base key of a registered IndexSpace
  /// @param space2_key base key of a registered IndexSpace
  /// @return true if `valid_unIon(retrieve(space1_key),retrieve(space2_key))`
  /// is registered
  template <basic_string_convertible S1, basic_string_convertible S2>
  bool valid_unIon(S1&& space1_key, S2&& space2_key) const {
    const auto& space1 = retrieve_or_throw(
        space1_key,
        "IndexBasisRegistry::valid_unIon(s1,s2): space with key s1=");
    const auto& space2 = retrieve_or_throw(
        space2_key, "IndexBasisRegistry::valid_unIon(s1,s2): space key s2=");
    return this->valid_unIon(space1, space2);
  }

  /// alias for valid_unIon
  template <index_space_or_label T1, index_space_or_label T2>
  bool valid_union(T1&& t1, T2&& t2) const {
    return valid_unIon(std::forward<T1>(t1), std::forward<T2>(t2));
  }

  /// @param space1
  /// @param space2
  /// @return the union of two spaces.
  /// @note can only return registered spaces
  /// @note never returns nullspace
  const IndexSpace& unIon(const IndexSpace& space1,
                          const IndexSpace& space2) const {
    if (!contains(space1))
      throw Exception(
          std::string("IndexBasisRegistry::valid_unIon(s1,s2): space s1=") +
          to_string(space1) + " must be added to the registry first");
    if (!contains(space2))
      throw Exception(
          std::string("IndexBasisRegistry::valid_unIon(s1,s2): space s2=") +
          to_string(space2) + " must be added to the registry first");

    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {
        throw Exception(
            std::string("IndexBasisRegistry::valid_unIon(s1,s2): spaces s1=") +
            to_string(space1) + " and s2=" + to_string(space2) +
            " must have identical "
            "quantum number attributes.");
      }
      auto unIontype = space1.type().unIon(space2.type());
      const IndexSpace& unIonSpace = find_by_attr({unIontype, space1.qns()});
      if (unIonSpace == IndexSpace::null) {
        throw Exception(
            std::string("IndexBasisRegistry::valid_unIon(s1,s2), s1=") +
            to_string(space1) + " and s2=" + to_string(space2) +
            ": the result is not registered, must register first.");
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
  template <basic_string_convertible S1, basic_string_convertible S2>
  const IndexSpace& unIon(S1&& space1_key, S2&& space2_key) const {
    const auto& space1 = retrieve_or_throw(
        space1_key, "IndexBasisRegistry::unIon(s1,s2): space with key s1=");
    const auto& space2 = retrieve_or_throw(
        space2_key, "IndexBasisRegistry::unIon(s1,s2): space key s2=");
    return this->unIon(space1, space2);
  }

  /// @name physical particle space structure introspection
  /// @{

  /// @brief sets the mask of QN attributes of `IndexSpace`s  that can be
  /// occupied by physical particles.

  /// Some states do not correspond to physical particles, but may be present
  /// in the registry. Such spaces are not considered when e.g. determining
  /// the base spaces.
  /// @param m the mask of QN attributes of `IndexSpace`s  that can be occupied
  ///        by physical particles.
  void physical_particle_attribute_mask(bitset_t m);

  /// @brief sets the mask of QN attributes of `IndexSpace`s
  /// that can be occupied by physical particles.

  /// Some states do not correspond to physical particles, but may be present
  /// in the registry. Such spaces are not considered when e.g. determining
  /// the base spaces.
  /// @param m the mask of QN attributes of `IndexSpace`s  that can be occupied
  ///        by physical particles.
  template <convertible_to_bitset T>
  void physical_particle_attribute_mask(T m) {
    this->physical_particle_attribute_mask(to_bitset(m));
  }

  /// @brief accesses the mask of QN attributes of `IndexSpace`s  that can be
  /// occupied by physical particles.

  /// Some states do not correspond to physical particles, but may be present
  /// in the registry. Such spaces are not considered when e.g. determining
  /// the base spaces.
  /// @return the mask of QN attributes of `IndexSpace`s  that can be occupied
  ///         by physical particles.
  bitset_t physical_particle_attribute_mask() const;

  /// @return the bits of \p qn used to specify states of physical particles.
  IndexSpace::QuantumNumbers physical_particle_attributes(
      IndexSpace::QuantumNumbers qn) const;

  /// @return the bits of \p qn not used to specify states of physical
  /// particles.
  IndexSpace::QuantumNumbers other_attributes(
      IndexSpace::QuantumNumbers qn) const;

  /// @brief returns the list of _basis_ IndexSpace::Type objects

  /// A base IndexSpace::Type object has 1 bit in its bitstring.
  /// @sa IndexBasisRegistry::is_base
  /// @return (memoized) set of base IndexSpace::Type objects, sorted in
  /// increasing order
  const std::vector<IndexSpace::Type>& base_space_types() const {
    std::scoped_lock guard{mtx_memoized_};
    if (!base_space_types_) {
      const SpacesView space_entries = spaces();
      auto types =
          space_entries |
          ranges::views::transform([](const auto& s) { return s.type(); }) |
          ranges::views::filter([](const auto& t) { return is_base(t); }) |
          ranges::views::unique | ranges::to_vector;
      ranges::sort(types, [](auto t1, auto t2) { return t1 < t2; });
      base_space_types_ =
          std::make_shared<std::vector<IndexSpace::Type>>(std::move(types));
    }
    return *base_space_types_;
  }

  /// @brief returns the list of _basis_ IndexSpace objects

  /// A base IndexSpace object has 1 bit in its type() bitstring.
  /// @sa IndexBasisRegistry::is_base
  /// @return (memoized) set of base IndexSpace objects, sorted in the order of
  /// increasing type()
  const std::vector<IndexSpace>& base_spaces() const {
    std::scoped_lock guard{mtx_memoized_};
    if (!base_spaces_) {
      const SpacesView space_entries = spaces();
      auto spaces = space_entries |
                    ranges::views::filter(
                        [this](const auto& s) { return this->is_base(s); }) |
                    ranges::views::unique | ranges::to_vector;
      ranges::sort(spaces,
                   [](auto s1, auto s2) { return s1.type() < s2.type(); });
      base_spaces_ =
          std::make_shared<std::vector<IndexSpace>>(std::move(spaces));
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
    return is_base(IS.type()) && !other_attributes(IS.qns());
  }

  /// @brief equivalent to `is_base(retrieve(space_key))`
  /// @param space_key space key
  /// @return true if the space registered with \p space_key is in the basis
  /// @sa base_spaces
  template <basic_string_convertible S>
  bool is_base(S&& space_key) const {
    return this->is_base(retrieve_or_throw(
        space_key, "IndexBasisRegistry::is_base(s): space with key s="));
  }

  /// @brief checks if an IndexSpace::Type is in the basis
  /// @param t IndexSpace::Type
  /// @return true if @p t is in the basis
  /// @sa space_type_basis
  static bool is_base(const IndexSpace::Type& t) {
    return std::has_single_bit(static_cast<std::uint32_t>(t.to_int32()));
  }

  /// @}

  /// @brief an @c IndexSpace is occupied with respect to the fermi vacuum or a
  /// subset of that space
  /// @note only makes sense to ask this if in a SingleProduct vacuum context.
  bool is_pure_occupied(const IndexSpace& IS) const {
    if (!IS) {
      return false;
    }
    // this introduces icky dependence on bit footprint of vacuum_occupied_space
    // Reduce qns to physical-particle (e.g. spin) attributes before the
    // occupancy lookup: non-particle traits (LCAO/factorization bits) do not
    // participate in occupancy. The non-throwing lookup returns a null
    // vacuum-occupied type for a non-physical auxiliary space (e.g. a
    // density-fitting or batching index that carries no particle/spin
    // character), so such a space is reported not pure-occupied rather than
    // triggering a throw.
    auto const vacocc_type =
        vacuum_occupied_type_or_null(physical_particle_attributes(IS.qns()));
    if (IS.type().to_int32() <= vacocc_type.to_int32()) {
      return true;
    } else {
      return false;
    }
  }

  /// @brief equivalent to `is_pure_occupied(retrieve(space_key))`
  /// @param space_key space key
  /// @return `is_pure_occupied(retrieve(space_key))`
  /// @sa base_spaces
  template <basic_string_convertible S>
  bool is_pure_occupied(S&& space_key) const {
    return this->is_pure_occupied(retrieve_or_throw(
        space_key,
        "IndexBasisRegistry::is_pure_occupied(s): space with key s="));
  }

  /// @brief all states are unoccupied in the fermi vacuum
  /// @note again, this only makes sense to ask if in a SingleProduct vacuum
  /// context.
  bool is_pure_unoccupied(const IndexSpace& IS) const {
    if (!IS) {
      return false;
    } else {
      // Q: would be better to express as IS is a subspace of
      // vacuum_unoccupied_space? Then subspaces that are part of neither (this
      // is not supposed to happen) will not cause an issue here
      // non-throwing lookup: a non-physical auxiliary space yields null
      // vacuum-occupied type, so this reports it as (trivially) unoccupied
      return !IS.type().intersection(
          vacuum_occupied_type_or_null(physical_particle_attributes(IS.qns())));
    }
  }

  /// @brief equivalent to `is_pure_unoccupied(retrieve(space_key))`
  /// @param space_key space key
  /// @return `is_pure_unoccupied(retrieve(space_key))`
  /// @sa base_spaces
  template <basic_string_convertible S>
  bool is_pure_unoccupied(S&& space_key) const {
    return this->is_pure_unoccupied(retrieve_or_throw(
        space_key,
        "IndexBasisRegistry::is_pure_unoccupied(s): space with key s="));
  }

  /// @brief some states are fermi vacuum occupied
  bool contains_occupied(const IndexSpace& IS) const {
    // non-throwing lookup on physical-particle attributes: a non-physical
    // auxiliary space (no particle/spin character) contains no occupied states
    return IS.type().intersection(vacuum_occupied_type_or_null(
               physical_particle_attributes(IS.qns()))) !=
           IndexSpace::Type::null;
  }

  /// @brief equivalent to `contains_occupied(retrieve(space_key))`
  /// @param space_key space key
  /// @return `contains_occupied(retrieve(space_key))`
  /// @sa base_spaces
  template <basic_string_convertible S>
  bool contains_occupied(S&& space_key) const {
    return this->contains_occupied(retrieve_or_throw(
        space_key,
        "IndexBasisRegistry::contains_occupied(s): space with key s="));
  }

  /// @brief some states are fermi vacuum unoccupied
  bool contains_unoccupied(const IndexSpace& IS) const {
    // non-throwing lookup on physical-particle attributes: a non-physical
    // auxiliary space (no particle/spin character) contains no unoccupied
    // states
    return IS.type().intersection(vacuum_unoccupied_type_or_null(
               physical_particle_attributes(IS.qns()))) !=
           IndexSpace::Type::null;
  }

  /// @brief equivalent to `contains_occupied(retrieve(space_key))`
  /// @param space_key space key
  /// @return `contains_occupied(retrieve(space_key))`
  /// @sa base_spaces
  template <basic_string_convertible S>
  bool contains_unoccupied(S&& space_key) const {
    return this->contains_unoccupied(retrieve_or_throw(
        space_key,
        "IndexBasisRegistry::contains_unoccupied(s): space with key s="));
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
  IndexBasisRegistry& vacuum_occupied_space(IndexSpace::Type t) {
    throw_if_missing(t, "vacuum_occupied_space");
    std::get<0>(vacocc_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  ///                nonzero occupancy in the reference wave function
  /// @return reference to `this`
  IndexBasisRegistry& vacuum_occupied_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "vacuum_occupied_space");
    std::get<1>(vacocc_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `vacuum_occupied_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexBasisRegistry& vacuum_occupied_space(const IndexSpace& s) {
    return vacuum_occupied_space(s.type());
  }

  /// equivalent to `vacuum_occupied_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <basic_string_convertible S>
  IndexBasisRegistry& vacuum_occupied_space(S&& l) {
    return vacuum_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in vacuum state for any set of quantum numbers
  /// @throw Exception if @p nulltype_ok is false and
  /// vacuum_occupied_space had not been specified
  const IndexSpace::Type& vacuum_occupied_space(
      bool nulltype_ok = false) const {
    if (!std::get<0>(vacocc_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw Exception(
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
  const IndexSpace& vacuum_occupied_space(IndexSpace::QuantumNumbers qn) const {
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
  IndexBasisRegistry& reference_occupied_space(IndexSpace::Type t) {
    throw_if_missing(t, "reference_occupied_space");
    std::get<0>(refocc_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  /// nonzero occupancy in
  ///          the reference wave function
  /// @return reference to `this`
  IndexBasisRegistry& reference_occupied_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "reference_occupied_space");
    std::get<1>(refocc_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `reference_occupied_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexBasisRegistry& reference_occupied_space(const IndexSpace& s) {
    return reference_occupied_space(s.type());
  }

  /// equivalent to `reference_occupied_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <basic_string_convertible S>
  IndexBasisRegistry& reference_occupied_space(S&& l) {
    return reference_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in reference state for any set of quantum
  /// numbers
  /// @throw Exception if @p nulltype_ok is false and
  /// reference_occupied_space had not been specified
  const IndexSpace::Type& reference_occupied_space(
      bool nulltype_ok = false) const {
    if (!std::get<0>(refocc_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw Exception(
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
      IndexSpace::QuantumNumbers qn) const {
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
  IndexBasisRegistry& complete_space(IndexSpace::Type s) {
    throw_if_missing(s, "complete_space");
    std::get<0>(complete_) = s;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  /// nonzero occupancy in
  ///          the reference wave function
  IndexBasisRegistry& complete_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "complete_space");
    std::get<1>(complete_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `complete_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexBasisRegistry& complete_space(const IndexSpace& s) {
    return complete_space(s.type());
  }

  /// equivalent to `complete_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <basic_string_convertible S>
  IndexBasisRegistry& complete_space(S&& l) {
    return complete_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the complete Hilbert space for any set of quantum numbers
  /// @throw Exception if @p nulltype_ok is false and complete_space
  /// had not been specified
  const IndexSpace::Type& complete_space(bool nulltype_ok = false) const {
    if (!std::get<0>(complete_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw Exception(
          "complete space has not been specified, call "
          "complete_space(IndexSpace::Type)");
    } else
      return std::get<0>(complete_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the complete Hilbert space for the given set of quantum numbers
  const IndexSpace& complete_space(IndexSpace::QuantumNumbers qn) const {
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
      IndexSpace::QuantumNumbers qn) const {
    auto complete_type = this->complete_space(qn).type();
    auto vacocc_type = this->vacuum_occupied_space(qn).type();
    auto vacuocc_type =
        complete_type.xOr(vacocc_type).intersection(complete_type);
    return this->retrieve(vacuocc_type, qn);
  }

  /// @return the space that is reference-occupied but vacuum-unoccupied: the
  /// partially occupied space of a Vacuum::MultiProduct reference (null if
  /// there is none)
  const IndexSpace& active_space(IndexSpace::QuantumNumbers qn) const {
    return this->intersection(this->reference_occupied_space(qn),
                              this->vacuum_unoccupied_space(qn));
  }

  /// @name specifies in which space holes can be created successfully from the
  /// reference wave function
  /// @note convenience for making operators
  /// @{

  /// @param t an IndexSpace::Type specifying where holes can be created;
  ///          to specify hole space per specific QN set use the other
  ///          overload
  IndexBasisRegistry& hole_space(IndexSpace::Type t) {
    throw_if_missing(t, "hole_space");
    std::get<0>(hole_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which holes
  /// can be created
  IndexBasisRegistry& hole_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "hole_space");
    std::get<1>(hole_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `hole_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexBasisRegistry& hole_space(const IndexSpace& s) {
    return hole_space(s.type());
  }

  /// equivalent to `hole_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <basic_string_convertible S>
  IndexBasisRegistry& hole_space(S&& l) {
    return hole_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which holes can be created
  /// @throw Exception if @p nulltype_ok is false and
  /// hole_space had not been specified
  const IndexSpace::Type& hole_space(bool nulltype_ok = false) const {
    if (!std::get<0>(hole_space_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw Exception(
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
  const IndexSpace& hole_space(IndexSpace::QuantumNumbers qn) const {
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
  IndexBasisRegistry& particle_space(IndexSpace::Type t) {
    throw_if_missing(t, "particle_space");
    std::get<0>(particle_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which
  /// particles can be created
  IndexBasisRegistry& particle_space(
      container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "particle_space");
    std::get<1>(particle_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `particle_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexBasisRegistry& particle_space(const IndexSpace& s) {
    return particle_space(s.type());
  }

  /// equivalent to `particle_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <basic_string_convertible S>
  IndexBasisRegistry& particle_space(S&& l) {
    return particle_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which particles can be created
  /// @throw Exception if @p nulltype_ok is false and
  /// particle_space had not been specified
  const IndexSpace::Type& particle_space(bool nulltype_ok = false) const {
    if (!std::get<0>(particle_space_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw Exception(
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
  const IndexSpace& particle_space(IndexSpace::QuantumNumbers qn) const {
    auto it = std::get<1>(particle_space_).find(qn);
    if (it != std::get<1>(particle_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->particle_space(), qn);
    }
  }

  /// @}

 private:
  table_type bases_;
  std::size_t named_count_ = 0;  // the number of named basis instances

  /// @throw Exception, naming @p caller, unless @p label is a valid label of
  /// a space or of a named basis instance, one that indices can be parsed
  /// with (see io::serialization::v1::is_base_key())
  static void validate_label(std::wstring_view label, std::string_view caller);

  static bool is_space(const table_type::value_type& e) {
    return !e.second.has_basis_instance();
  }

  /// the IndexSpace of a non-const table element, for metadata writes (the
  /// element is not const, so the cast is well-defined; metadata is excluded
  /// from IndexBasis equality and from the table key)
  static IndexSpace& mutable_space(IndexBasis& b) {
    return const_cast<IndexSpace&>(b.space());
  }

  /// @return the space of the space entry @p entry and the space copies of
  /// every named basis instance of that space, for metadata writes
  container::svector<std::reference_wrapper<IndexSpace>> space_copies_of(
      IndexBasis& entry) {
    SEQUANT_ASSERT(!entry.has_basis_instance());
    container::svector<std::reference_wrapper<IndexSpace>> result;
    result.emplace_back(mutable_space(entry));
    for (auto& [label, basis] : bases_)
      if (basis.has_basis_instance() && basis.space() == entry.space())
        result.emplace_back(mutable_space(basis));
    return result;
  }

  /// @return the entry registered under @p label
  /// @throw IndexSpace::bad_key if no entry is registered under @p label
  template <basic_string_convertible S>
  IndexBasis& entry_or_throw(S&& label) {
    auto it = bases_.find(IndexSpace::reduce_key(to_basic_string_view(label)));
    if (it == bases_.end()) throw IndexSpace::bad_key(label);
    return it->second;
  }

  /// @return the dimension (the integral argument) and the Field among
  /// @p args, std::nullopt for each that is not given
  template <typename... Args>
  static std::pair<std::optional<unsigned long>, std::optional<Field>>
  parse_size_and_field(const Args&... args) {
    std::pair<std::optional<unsigned long>, std::optional<Field>> result;
    auto h_args = boost::hana::make_tuple(args...);

    auto h_ints = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::traits::is_integral(boost::hana::decltype_(arg));
    });
    constexpr auto nints = boost::hana::size(h_ints);
    static_assert(
        nints == boost::hana::size_c<0> || nints == boost::hana::size_c<1>,
        "IndexBasisRegistry::add: only one integral argument is allowed");
    if constexpr (nints == boost::hana::size_c<1>) {
      result.first = boost::hana::at_c<0>(h_ints);
    }

    auto h_field = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::type_c<decltype(arg)> == boost::hana::type_c<Field>;
    });
    constexpr auto nfields = boost::hana::size(h_field);
    static_assert(
        nfields == boost::hana::size_c<0> || nfields == boost::hana::size_c<1>,
        "IndexBasisRegistry::add: only one Field argument is allowed");
    if constexpr (nfields == boost::hana::size_c<1>) {
      result.second = boost::hana::at_c<0>(h_field);
    }

    return result;
  }

  bitset_t physical_particle_attribute_mask_ = bitset::null;

  // memoized data
  mutable std::shared_ptr<std::vector<IndexSpace::Type>> base_space_types_;
  mutable std::shared_ptr<std::vector<IndexSpace>> base_spaces_;
  mutable std::recursive_mutex
      mtx_memoized_;  // guards every access to the memoized data
  IndexBasisRegistry& clear_memoized_data_and_return_this() {
    std::scoped_lock guard{mtx_memoized_};
    base_space_types_.reset();
    base_spaces_.reset();
    return *this;
  }

  /// @brief non-throwing counterpart of `vacuum_occupied_space(qn).type()`
  /// @param qn quantum numbers (typically already reduced to
  ///        physical-particle attributes)
  /// @return the vacuum-occupied type registered for @p qn (per-qns override if
  ///         present, else the default), or IndexSpace::Type::null when no
  ///         occupied space is registered at @p qn -- e.g. a non-physical
  ///         auxiliary space (density-fitting, batching) that carries no
  ///         particle/spin character, hence no occupancy.
  IndexSpace::Type vacuum_occupied_type_or_null(
      IndexSpace::QuantumNumbers qn) const {
    auto const& qn2type = std::get<1>(vacocc_);
    auto const it = qn2type.find(qn);
    IndexSpace::Type const t =
        (it != qn2type.end()) ? it->second : std::get<0>(vacocc_);
    if (!t) return IndexSpace::Type::null;
    return retrieve_ptr(t, qn) ? t : IndexSpace::Type::null;
  }

  /// @brief non-throwing counterpart of `complete_space(qn).type()`
  /// @param qn quantum numbers
  /// @return the complete-space type registered for @p qn, or
  ///         IndexSpace::Type::null when none is registered at @p qn
  IndexSpace::Type complete_type_or_null(IndexSpace::QuantumNumbers qn) const {
    auto const& qn2type = std::get<1>(complete_);
    auto const it = qn2type.find(qn);
    IndexSpace::Type const t =
        (it != qn2type.end()) ? it->second : std::get<0>(complete_);
    if (!t) return IndexSpace::Type::null;
    return retrieve_ptr(t, qn) ? t : IndexSpace::Type::null;
  }

  /// @brief non-throwing counterpart of `vacuum_unoccupied_space(qn).type()`
  /// @param qn quantum numbers (typically already reduced to
  ///        physical-particle attributes)
  /// @return the vacuum-unoccupied type (complete minus vacuum-occupied) for
  ///         @p qn, or IndexSpace::Type::null when no complete space is
  ///         registered at @p qn (non-physical auxiliary space)
  IndexSpace::Type vacuum_unoccupied_type_or_null(
      IndexSpace::QuantumNumbers qn) const {
    auto const complete_t = complete_type_or_null(qn);
    if (!complete_t) return IndexSpace::Type::null;
    auto const vacocc_t = vacuum_occupied_type_or_null(qn);
    return complete_t.xOr(vacocc_t).intersection(complete_t);
  }

  /// @return the IndexSpace registered under @p space_key
  /// @throw Exception if @p space_key is not registered; its message is
  /// @p message_prefix followed by the key
  template <basic_string_convertible S>
  const IndexSpace& retrieve_or_throw(const S& space_key,
                                      std::string_view message_prefix) const {
    const auto* ptr = this->retrieve_ptr(space_key);
    if (!ptr)
      throw Exception(std::string(message_prefix) + toUtf8(space_key) +
                      " must be added to the registry first");
    return *ptr;
  }

  /// @brief find an IndexSpace from its attr. return nullspace if not present.
  /// @param attr the attribute of the IndexSpace
  const IndexSpace& find_by_attr(const IndexSpace::Attr& attr) const {
    const auto* ptr = retrieve_ptr(attr);
    return ptr ? *ptr : IndexSpace::null;
  }

  void throw_if_missing(const IndexSpace::Type& t,
                        const IndexSpace::QuantumNumbers& qn,
                        std::string call_context = "") {
    if (retrieve_ptr(t, qn)) return;
    throw Exception(
        call_context + ": missing { IndexSpace::Type=" + to_string(t) +
        " , IndexSpace::QuantumNumbers=" + to_string(qn) + " } combination");
  }

  // same as above, but ignoring qn
  void throw_if_missing(const IndexSpace::Type& t,
                        std::string call_context = "") {
    for (auto&& space : spaces()) {
      if (space.type() == t) {
        return;
      }
    }
    throw Exception(call_context + ": missing { IndexSpace::Type=" +
                    to_string(t) + " , any IndexSpace::QuantumNumbers } space");
  }

  void throw_if_missing_any(const container::map<IndexSpace::QuantumNumbers,
                                                 IndexSpace::Type>& qn2type,
                            std::string call_context = "") {
    container::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type_found;
    for (auto&& space : spaces()) {
      for (auto&& [qn, t] : qn2type) {
        if (space.type() == t && space.qns() == qn) {
          [[maybe_unused]] auto [it, inserted] =
              qn2type_found.try_emplace(qn, t);
          SEQUANT_ASSERT(inserted);
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
    throw Exception(errmsg);
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
                      "IndexBasisRegistry::add{,_union,_intersect}: unknown "
                      "attribute tag");
      }
    });
  }

  /// @brief computes the dimension of the space

  /// for a base space return its extent, for a composite space compute as a sum
  /// of extents of base subspaces
  /// @param space_attr the IndexSpace attribute
  /// @return the dimension of the space
  unsigned long compute_dimension(const IndexSpace::Attr& space_attr) const {
    if (is_base(space_attr.type())) {
      return this->retrieve(space_attr).dimension();
    } else {
      // compute_dimension is used when populating the registry
      // so don't use base_spaces() here
      const SpacesView space_entries = spaces();
      unsigned long size = ranges::accumulate(
          space_entries | ranges::views::filter([this, &space_attr](auto& s) {
            return s.qns() == space_attr.qns() && this->is_base(s.type()) &&
                   space_attr.type().intersection(s.type());
          }),
          0ul, [](unsigned long size, const IndexSpace& s) {
            return size + s.dimension();
          });
      return size;
    }
  }

  Field compute_field(const IndexSpace::Attr& space_attr) const {
    if (is_base(space_attr.type())) {
      return this->retrieve(space_attr).field();
    }

    // compute_field is used when populating the registry
    // so don't use base_spaces() here

    const SpacesView space_entries = spaces();
    bool contains_complex = std::ranges::any_of(
        space_entries |
            std::ranges::views::filter([this, &space_attr](auto& s) {
              return s.qns() == space_attr.qns() && this->is_base(s.type()) &&
                     space_attr.type().intersection(s.type());
            }),
        [](const IndexSpace& s) { return s.field() == Field::Complex; });

    return contains_complex ? Field::Complex : Field::Real;
  }

  /// registries are equal if they have equal entries (spaces and named basis
  /// instances, under equal labels), of equal dimension, extent, metric and
  /// field, and specify the same physical-particle
  /// attributes and vacuum-occupied, reference-occupied, complete, hole and
  /// particle spaces
  friend bool operator==(const IndexBasisRegistry& isr1,
                         const IndexBasisRegistry& isr2) {
    // IndexBasis equality ignores the metadata
    return std::ranges::equal(
               isr1.bases_, isr2.bases_,
               [](const auto& e1, const auto& e2) {
                 return e1.first == e2.first && e1.second == e2.second &&
                        e1.second.space().dimension() ==
                            e2.second.space().dimension() &&
                        e1.second.extent() == e2.second.extent() &&
                        e1.second.metric() == e2.second.metric() &&
                        e1.second.field() == e2.second.field();
               }) &&
           isr1.physical_particle_attribute_mask_ ==
               isr2.physical_particle_attribute_mask_ &&
           isr1.vacocc_ == isr2.vacocc_ && isr1.refocc_ == isr2.refocc_ &&
           isr1.complete_ == isr2.complete_ &&
           isr1.hole_space_ == isr2.hole_space_ &&
           isr1.particle_space_ == isr2.particle_space_;
  }
};  // class IndexBasisRegistry

using IndexSpaceRegistry [[deprecated("use IndexBasisRegistry")]] =
    IndexBasisRegistry;

}  // namespace sequant
#endif  // SEQUANT_INDEX_BASIS_REGISTRY_HPP
