//
// Created by Eduard Valeyev on 2025-24-07.
//

#ifndef SEQUANT_TENSOR_NETWORK_V3_H
#define SEQUANT_TENSOR_NETWORK_V3_H

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network/canonicals.hpp>
#include <SeQuant/core/tensor_network/slot.hpp>
#include <SeQuant/core/tensor_network/vertex.hpp>

#include <range/v3/range/traits.hpp>

#include <cassert>
#include <cstdlib>
#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// forward declarations
namespace bliss {
class Graph;
}

namespace sequant {

/// @brief A (non-directed) graph view of a sequence of AbstractTensor objects

/// @note The main role of this is to canonize itself. Since Tensor objects can
/// be connected by multiple Index'es (thus edges are colored), what is
/// canonized is actually the graph of indices (roughly the dual of the tensor
/// graph), with Tensor objects represented by one or more vertices.
class TensorNetworkV3 {
 public:
  /// @return the implementation version of TN
  constexpr static int version() { return 3; }

  // for unit testing only
  friend class TensorNetworkV3Accessor;

  using Origin = SlotType;

  class Vertex {
   public:
    Vertex(Origin origin, std::size_t terminal_idx, std::size_t index_slot,
           Symmetry terminal_symm);

    Origin getOrigin() const;
    std::size_t getTerminalIndex() const;
    std::size_t getIndexSlot() const;
    Symmetry getTerminalSymmetry() const;

    bool operator<(const Vertex &rhs) const;
    bool operator==(const Vertex &rhs) const;

   private:
    Origin origin;
    std::size_t terminal_idx;
    std::size_t index_slot;
    Symmetry terminal_symm;
  };

  // clang-format off
  /// @brief a (hyper)edge in a tensor network

  /// Edge in a TensorNetworkV3 = the Index annotating it +
  /// a list of vertices corresponding to the Tensor index slots it connects
  /// @note this is move-only since using pointers to refer to Index objects
  // clang-format on
  class Edge {
   public:
    Edge() = default;
    Edge(const Edge &) = delete;
    Edge(Edge &&) = default;
    Edge &operator=(const Edge &) = delete;
    Edge &operator=(Edge &&) = default;
    explicit Edge(const Vertex &vertex) : vertices{vertex} {}
    explicit Edge(std::initializer_list<Vertex> vertices) {
      ranges::for_each(vertices,
                       [this](const Vertex &v) { this->connect_to(v); });
    }
    Edge(const Vertex &vertex, const Index *index)
        : vertices{vertex}, index(index) {}
    Edge(std::initializer_list<Vertex> vertices, const Index *index)
        : Edge(vertices) {
      this->index = index;
    }

    Edge &connect_to(const Vertex &vertex) {
      // free Edge
      if (vertices.empty()) {
        vertices.emplace(vertex);
      } else {
        // - cannot connect braket slots to aux slots
        auto &first = *(vertices.begin());
        if ((first.getOrigin() == Origin::Aux &&
             vertex.getOrigin() != Origin::Aux) ||
            (first.getOrigin() != Origin::Aux &&
             vertex.getOrigin() == Origin::Aux)) {
          throw std::invalid_argument(
              "TensorNetworkV3::Edge::connect_to: aux slot cannot be connected "
              "to a non-aux slot");
        }
        // - can connect bra slot to ket slot, and vice versa, unless there is
        // no distinction between primal and dual spaces
        if (get_default_context().braket_symmetry() != BraKetSymmetry::symm) {
          if (first.getOrigin() == Origin::Bra &&
              vertex.getOrigin() != Origin::Ket) {
            throw std::invalid_argument(
                "TensorNetworkV3::Edge::connect_to: bra slot can only be "
                "connected "
                "to a ket slot if default context's braket_symmetry() != "
                "BraKetSymmetry::symm");
          }
          if (first.getOrigin() == Origin::Ket &&
              vertex.getOrigin() != Origin::Bra) {
            throw std::invalid_argument(
                "TensorNetworkV3::Edge::connect_to: ket slot can only be "
                "connected "
                "to a bra slot if default context's braket_symmetry() != "
                "BraKetSymmetry::symm");
          }
        }
        add_vertex(vertex);
      }

      return *this;
    }

    bool operator<(const Edge &other) const {
      if (vertex_count() != other.vertex_count()) {
        // Ensure external indices (edges that are only attached to a tensor on
        // one side) always come before internal ones
        return vertex_count() < other.vertex_count();
      }

      if (!(vertices == other.vertices)) {
        return vertices < other.vertices;
      }

      assert(index && other.index);
      return index->space() < other.index->space();
    }

    bool operator==(const Edge &other) const {
      return vertices == other.vertices;
    }

    /// @param i vertex ordinal
    /// @return const reference to the `i`th Vertex object
    /// @pre `this->size() > i`
    const Vertex &vertex(std::size_t i) const {
      assert(vertices.size() > i);
      return *(vertices.begin() + i);
    }

    /// @return the number of attached terminals (0 or more)
    std::size_t vertex_count() const { return vertices.size(); }

    const Index &idx() const {
      assert(index);
      return *index;
    }

   private:
    container::set<Vertex> vertices;
    const Index *index = nullptr;

    /// @param vertex a vertex to be added
    /// @throw std::invalid_argument if @p vertex is already connected by this
    /// Edge
    void add_vertex(const Vertex &vertex) {
      auto [it, inserted] = this->vertices.emplace(vertex);
      if (!inserted)
        throw std::invalid_argument(
            "TensorNetworkV3::Edge::add_vertex(v): v is already connected by "
            "this Edge");
    }
  };

  struct Graph {
    /// The type used to encode the color of a vertex. The restriction of this
    /// being as 32-bit integer comes from how BLISS is trying to convert these
    /// into RGB values.
    using VertexColor = std::uint32_t;

    std::unique_ptr<bliss::Graph> bliss_graph;
    std::vector<std::wstring> vertex_labels;
    std::vector<std::optional<std::wstring>> vertex_texlabels;
    std::vector<VertexColor> vertex_colors;
    std::vector<VertexType> vertex_types;
    container::map<Index, std::size_t> idx_to_vertex;

    Graph() = default;

    std::size_t vertex_to_index_idx(std::size_t vertex) const;
    /// maps vertex ordinal to tensor cluster ordinal
    /// @note usable as bliss::Graph::DotOptions::vertex_to_subgraph
    std::optional<std::size_t> vertex_to_tensor_idx(std::size_t vertex) const;
  };

  TensorNetworkV3(const Expr &expr) {
    if (expr.size() > 0) {
      for (const ExprPtr &subexpr : expr) {
        add_expr(*subexpr);
      }
    } else {
      add_expr(expr);
    }

    init_edges();
  }

  TensorNetworkV3(const ExprPtr &expr) : TensorNetworkV3(*expr) {}

  template <
      typename ExprPtrRange,
      typename = std::enable_if_t<!std::is_base_of_v<ExprPtr, ExprPtrRange> &&
                                  !std::is_base_of_v<Expr, ExprPtrRange>>>
  TensorNetworkV3(const ExprPtrRange &exprptr_range) {
    static_assert(
        std::is_base_of_v<ExprPtr, ranges::range_value_t<ExprPtrRange>>);
    for (const ExprPtr &current : exprptr_range) {
      add_expr(*current);
    }

    init_edges();
  }

  TensorNetworkV3(TensorNetworkV3 &&) noexcept;
  TensorNetworkV3 &operator=(TensorNetworkV3 &&) noexcept;

  /// copy constructor
  /// @warning does not copy edges
  TensorNetworkV3(const TensorNetworkV3 &other);

  /// copy assignment
  /// @warning does not copy edges
  TensorNetworkV3 &operator=(const TensorNetworkV3 &other) noexcept;

  /// @return const reference to the sequence of tensors
  /// @note after invoking TensorNetwork::canonicalize() the order of
  /// tensors may be different from that provided as input; use
  /// tensor_input_ordinals() to obtain the input ordinals of
  /// the tensors in the result
  const auto &tensors() const { return tensors_; }

  const auto &tensor_input_ordinals() const { return tensor_input_ordinals_; }

  using NamedIndexSet = container::set<Index, Index::FullLabelCompare>;

  /// @param cardinal_tensor_labels move all tensors with these labels to the
  /// front before canonicalizing indices
  /// @param fast if true (default), does fast canonicalization that is only
  /// optimal if all tensors are distinct; set to false to perform complete
  /// canonicalization
  /// @param named_indices specifies the indices that cannot be renamed, i.e.
  /// their labels are meaningful; default is nullptr, which results in external
  /// indices treated as named indices
  /// @return byproduct of canonicalization (e.g. phase); if none, returns
  /// nullptr
  ExprPtr canonicalize(
      const container::vector<std::wstring> &cardinal_tensor_labels = {},
      bool fast = true, const NamedIndexSet *named_indices = nullptr);

  /// metadata produced by canonicalize_slots()
  struct SlotCanonicalizationMetadata {
    /// list of named indices
    NamedIndexSet named_indices;

    /// type of less-than comparison function for named indices, receives {Index
    /// ptr, its slot type}
    using named_index_compare_t =
        std::function<bool(const std::pair<const Index *, IndexSlotType> &,
                           const std::pair<const Index *, IndexSlotType> &)>;

    /// less-than comparison function for named indices, used for
    /// coarse-grained sorting of named indices,
    /// before sorting to canonical order
    named_index_compare_t named_index_compare;

    /// list of named indices, sorted first by named_index_compare,
    /// then by canonical order; iterators point to named_indices
    container::svector<NamedIndexSet::const_iterator> named_indices_canonical;

    /// canonicalized colored graph, use graph->cmp to compare against another
    /// to detect equivalence
    std::shared_ptr<bliss::Graph> graph;

    [[nodiscard]] size_t hash_value() const;

    [[nodiscard]] inline auto get_index_view() const {
      return named_indices_canonical  //
             | ranges::views::indirect;
    }

    template <typename Cont>
    auto get_indices() const {
      return get_index_view() | ranges::to<Cont>;
    }

    /// if tensor network contains tensors with antisymmetric bra/ket this
    /// reports the phase change due to permutation of slots relative to their
    /// input order
    std::int8_t phase = +1;  // +1 or -1
  };

  /// Like canonicalize(), but only use graph-based canonicalization to
  /// produce canonical list of slots occupied by named indices.
  /// This is sufficient to be able to match 2 tensor networks that
  /// differ in anonymous and named indices.
  /// @param cardinal_tensor_labels move all tensors with these labels to the
  /// front before canonicalizing indices
  /// @param named_indices specifies the indices that cannot be renamed, i.e.
  /// their labels are meaningful; default is nullptr, which results in external
  /// indices treated as named indices
  /// @param named_index_compare less-than comparison function for
  /// named indices, used for coarse-grained sorting of named indices,
  /// before sorting to canonical order; the default is to sort
  /// by Index::space()
  /// @return the computed canonicalization metadata
  SlotCanonicalizationMetadata canonicalize_slots(
      const container::vector<std::wstring> &cardinal_tensor_labels = {},
      const NamedIndexSet *named_indices = nullptr,
      SlotCanonicalizationMetadata::named_index_compare_t named_index_compare =
          default_idxptr_slottype_lesscompare{});

  /// Factorizes tensor network
  /// @return sequence of binary products; each element encodes the tensors to
  /// be
  ///         multiplied (values >0 refer to the tensors in tensors(),
  ///         values <0 refer to the elements of this sequence. E.g. sequences
  ///         @c {{0,1},{-1,2},{-2,3}} , @c {{0,2},{1,3},{-1,-2}} , @c
  ///         {{3,1},{2,-1},{0,-2}} encode the following respective
  ///         factorizations @c (((T0*T1)*T2)*T3) , @c ((T0*T2)*(T1*T3)) , and
  ///         @c (((T3*T1)*T2)*T0) .
  container::svector<std::pair<long, long>> factorize();

  /// accessor for the Edge object sequence
  /// @return const reference to the sequence container of Edge objects, sorted
  /// by their Index's full label
  /// @sa Edge
  const auto &edges() const {
    assert(have_edges_);
    return edges_;
  }

  /// @brief Returns a range of external indices, i.e. those indices that do not
  /// connect tensors
  const auto &ext_indices() const {
    assert(have_edges_);
    return ext_indices_;
  }

  /// options for generating Graph from an object of this type
  struct CreateGraphOptions {
    /// pointer to the set of named indices (ordinarily,
    /// this includes all external indices);
    /// default is nullptr, which means use all external indices for
    /// named indices
    const NamedIndexSet *named_indices = nullptr;

    /// if false, will use same color for all
    /// named indices that have same Index::color(), else will use distinct
    /// color for each
    bool distinct_named_indices = false;

    /// if false, will not generate the labels
    bool make_labels = true;

    /// if false, will not generate the TeX labels
    bool make_texlabels = true;

    /// if false, will not generate the Index->vertex map
    bool make_idx_to_vertex = false;
  };

  // clang-format off
  /// @brief converts the network into a Bliss graph whose vertices are indices
  /// and tensor vertex representations
  /// @param[in] options the options for generating the graph
  /// @return The created Graph object

  /// @note Rules for constructing the graph:
  ///   - symmetries are encoded by topology and color
  ///   - vertex is introduced for an index slot or a bundle thereof, an index, and tensors core
  ///     - bundles of slots include:
  ///        - bra (bundle of 1 or more bra index slots)
  ///        - ket
  ///        - particle bundle (column of slots in covariant notation, e.g. a bundle of a bra index slot and ket index slot)
  ///        - protoindex bundle (bundle of index slots attached to an Index)
  ///     - N.B. lack of symmetry between slots and indices (can represent bundles of slots, but not bundles of indices) is due to the fact that index is already a plain index or an index with protoindices (i.e., a collection of plain indices)
  ///     - to create a bundle of objects create the new "bundle" vertex and connect it to each object's vertex
  ///       - if set of n objects has symmetry with respect to permutation described by a particular irrep of S_n the colors of the objects' vertices must be the same (else they are distinguishable)
  ///       - for 2 objects it is not necessary to introduce a bundle vertex, but bundling is often done even for n=1 for the sake of consistency
  ///     - vertices that can be swapped should have the same color
  ///       - indices are colored by their space (not by label/ordinal) + the colors of their protoindices
  ///       - slots and their bundles use custom colors (see below)
  ///
  ///   Consider Tensor as an example:
  ///     - each index slot of a tensor is a vertex (b+k+a such vertices in
  ///     an order-{b,k,a} tensor)
  ///     - bra/ket slot vertices of an antisymmetric/symmetric tensor are
  ///     bundled into bra/ket vertex (2 such vertices); they are subsequently bundleed into a braket vertex (1 such vertex)
  ///     - if Tensor has asymmetric bras/kets each matching
  ///     (corresponding to same particle ordinal) bra/ket slot vertex pair is
  ///     bundleed into a braket vertex
  ///     (max(b,k) such indices)
  ///     - braket vertices bundleed into tensor core vertex
  ///     - bra/ket slot vertices have same color for antisymmetric/symmetric
  ///     tensors
  ///     - for asymmetric bra/ket bra+ket+braket bundles must have same color
  ///     if tensor is particle-symmetric, else
  ///     - for tensors with bra<->ket symmetry matching bra and ket slot
  ///     vertices have identical colors.
  // clang-format on
  Graph create_graph(const CreateGraphOptions &options = {
                         .named_indices = nullptr,
                         .distinct_named_indices = false,
                         .make_labels = true,
                         .make_texlabels = true,
                         .make_idx_to_vertex = false}) const;

 private:
  /// list of tensors
  /// - before canonicalize(): input
  /// - after canonicalize(): canonical
  container::svector<AbstractTensorPtr> tensors_;
  /// input ordinals of tensors in tensors_
  container::svector<std::size_t> tensor_input_ordinals_;

  container::vector<Edge> edges_;
  bool have_edges_ = false;
  /// ext indices do not connect tensors
  /// sorted by full label of the corresponding value (Index)
  /// N.B. this may contain some indices in pure_proto_indices_ if there are
  /// external indices that depend on them
  NamedIndexSet ext_indices_;
  /// some proto indices may not be in edges_ if they appear exclusively among
  /// proto indices
  /// @note these will need to be processed separately from the rest
  /// to appear as vertices on the graph
  NamedIndexSet pure_proto_indices_;

  /// initializes edges_, ext_indices_, and pure_proto_indices_
  void init_edges();

  /// Canonicalizes the network graph representation using colored graph
  /// canonicalization
  /// @return The byproduct of canonicalization
  /// @note this produces canonical representation that is invariant with
  /// respect to the renaming of named indices
  [[nodiscard]] ExprPtr canonicalize_graph(const NamedIndexSet &named_indices);

  /// Canonicalizes every individual tensor for itself, taking into account only
  /// tensor blocks
  /// @returns The byproduct of the canonicalizations
  ExprPtr canonicalize_individual_tensor_blocks(
      const NamedIndexSet &named_indices);

  /// Canonicalizes every individual tensor for itself
  /// @returns The byproduct of the canonicalizations
  ExprPtr canonicalize_individual_tensors(const NamedIndexSet &named_indices);

  ExprPtr do_individual_canonicalization(
      const TensorCanonicalizer &canonicalizer);

  void add_expr(const Expr &expr) {
    ExprPtr clone = expr.clone();

    auto tensor_ptr = std::dynamic_pointer_cast<AbstractTensor>(clone);
    if (!tensor_ptr) {
      throw std::invalid_argument(
          "TensorNetworkV3::TensorNetworkV3: tried to add non-tensor to "
          "network");
    }

    tensors_.emplace_back(std::move(tensor_ptr));
    tensor_input_ordinals_.push_back(tensor_input_ordinals_.size());
  }
};

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(
    std::basic_ostream<CharT, Traits> &stream, SlotType origin) {
  switch (origin) {
    case SlotType::Bra:
      stream << "Bra";
      break;
    case SlotType::Ket:
      stream << "Ket";
      break;
    case SlotType::Aux:
      stream << "Aux";
      break;
  }
  return stream;
}

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_H
