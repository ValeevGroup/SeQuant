//
// Created by Eduard Valeyev on 2019-02-02.
//

#ifndef SEQUANT_TENSOR_NETWORK_V2_H
#define SEQUANT_TENSOR_NETWORK_V2_H

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/vertex_type.hpp>

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
class TensorNetworkV2 {
 public:
  friend class TensorNetworkV2Accessor;

  constexpr static size_t max_rank = 256;

  enum class Origin {
    Bra = 1,
    Ket,
    Aux,
  };

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
  /// @brief Edge in a TensorNetworkV2 = the Index annotating it + a pair of indices to identify which Tensor terminals it's connected to

  /// @note tensor terminals in a sequence of tensors are indexed as follows:
  /// - >0 for bra terminals (i.e. "+7" indicated connection to a bra terminal
  /// of 7th tensor object in the sequence)
  /// - <0 for ket terminals
  /// - 0 if free (not attached to any tensor objects)
  /// - position records the terminal's location in the sequence of bra/ket
  /// terminals (always 0 for symmetric/antisymmetric tensors) Terminal indices
  /// are sorted by the tensor index (i.e. by the absolute value of the terminal
  /// index), followed by position
  // clang-format on
  class Edge {
   public:
    Edge() = default;
    explicit Edge(Vertex vertex) : first(std::move(vertex)), second() {}
    Edge(Vertex vertex, Index index)
        : first(std::move(vertex)), second(), index(std::move(index)) {}

    Edge &connect_to(Vertex vertex) {
      assert(!second.has_value());

      if (!first.has_value()) {
        // unconnected Edge
        first = std::move(vertex);
      } else {
        second = std::move(vertex);
        if (second < first) {
          // Ensure first <= second
          std::swap(first, second);
        }
      }
      return *this;
    }

    bool operator<(const Edge &other) const {
      if (vertex_count() != other.vertex_count()) {
        // Ensure external indices (edges that are only attached to a tensor on
        // one side) always come before internal ones
        return vertex_count() < other.vertex_count();
      }

      if (!(first == other.first)) {
        return first < other.first;
      }

      if (second < other.second) {
        return second < other.second;
      }

      return index.space() < other.index.space();
    }

    bool operator==(const Edge &other) const {
      return first == other.first && second == other.second;
    }

    const Vertex &first_vertex() const {
      assert(first.has_value());
      return first.value();
    }
    const Vertex &second_vertex() const {
      assert(second.has_value());
      return second.value();
    }

    /// @return the number of attached terminals (0, 1, or 2)
    std::size_t vertex_count() const {
      return second.has_value() ? 2 : (first.has_value() ? 1 : 0);
    }

    const Index &idx() const { return index; }

   private:
    std::optional<Vertex> first;
    std::optional<Vertex> second;
    Index index;
  };

  struct Graph {
    std::unique_ptr<bliss::Graph> bliss_graph;
    std::vector<std::wstring> vertex_labels;
    std::vector<std::size_t> vertex_colors;
    std::vector<VertexType> vertex_types;

    Graph() = default;

    std::size_t vertex_to_index_idx(std::size_t vertex) const;
    std::size_t vertex_to_tensor_idx(std::size_t vertex) const;
  };

  TensorNetworkV2(const Expr &expr) {
    if (expr.size() > 0) {
      for (const ExprPtr &subexpr : expr) {
        add_expr(*subexpr);
      }
    } else {
      add_expr(expr);
    }

    init_edges();
  }

  TensorNetworkV2(const ExprPtr &expr) : TensorNetworkV2(*expr) {}

  template <
      typename ExprPtrRange,
      typename = std::enable_if_t<!std::is_base_of_v<ExprPtr, ExprPtrRange> &&
                                  !std::is_base_of_v<Expr, ExprPtrRange>>>
  TensorNetworkV2(const ExprPtrRange &exprptr_range) {
    static_assert(
        std::is_base_of_v<ExprPtr, typename ExprPtrRange::value_type>);
    for (const ExprPtr &current : exprptr_range) {
      add_expr(*current);
    }

    init_edges();
  }

  /// @return const reference to the sequence of tensors
  /// @note the order of tensors may be different from that provided as input
  const auto &tensors() const { return tensors_; }

  using named_indices_t = container::set<Index, Index::LabelCompare>;

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
      bool fast = true, const named_indices_t *named_indices = nullptr);

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

  /// @brief converts the network into a Bliss graph whose vertices are indices
  /// and tensor vertex representations
  /// @param[in] named_indices pointer to the set of named indices (ordinarily,
  /// this includes all external indices);
  ///            default is nullptr, which means use all external indices for
  ///            named indices
  /// @return The created Graph object

  /// @note Rules for constructing the graph:
  ///   - Indices with protoindices are connected to their protoindices,
  ///   either directly or (if protoindices are symmetric) via a protoindex
  ///   vertex.
  ///   - Indices are colored by their space, which in general encodes also
  ///   the space of the protoindices.
  ///   - An anti/symmetric n-body tensor has 2 terminals, each connected to
  ///   each other + to n index vertices.
  ///   - A nonsymmetric n-body tensor has n terminals, each connected to 2
  ///   indices and 1 tensor vertex which is connected to all n terminal
  ///   indices.
  ///   - tensor vertices are colored by the label+rank+symmetry of the
  ///   tensor; terminal vertices are colored by the color of its tensor,
  ///     with the color of symm/antisymm terminals augmented by the
  ///     terminal's type (bra/ket).
  Graph create_graph(const named_indices_t *named_indices = nullptr) const;

 private:
  // source tensors and indices
  container::svector<AbstractTensorPtr> tensors_;

  container::vector<Edge> edges_;
  bool have_edges_ = false;
  // ext indices do not connect tensors
  // sorted by *label* (not full label) of the corresponding value (Index)
  // this ensures that proto indices are not considered and all internal indices
  // have unique labels (not full labels)
  named_indices_t ext_indices_;

  /// initializes edges_ and ext_indices_
  void init_edges();

  /// Canonicalizes the network graph representation
  /// Note: The explicit order of tensors and labelling of indices
  /// remains undefined.
  void canonicalize_graph(const named_indices_t &named_indices);

  /// Canonicalizes every individual tensor for itself, taking into account only
  /// tensor blocks
  /// @returns The byproduct of the canonicalizations
  ExprPtr canonicalize_individual_tensor_blocks(
      const named_indices_t &named_indices);

  /// Canonicalizes every individual tensor for itself
  /// @returns The byproduct of the canonicalizations
  ExprPtr canonicalize_individual_tensors(const named_indices_t &named_indices);

  ExprPtr do_individual_canonicalization(
      const TensorCanonicalizer &canonicalizer);

  void add_expr(const Expr &expr) {
    ExprPtr clone = expr.clone();

    auto tensor_ptr = std::dynamic_pointer_cast<AbstractTensor>(clone);
    if (!tensor_ptr) {
      throw std::invalid_argument(
          "TensorNetworkV2::TensorNetworkV2: tried to add non-tensor to "
          "network");
    }

    tensors_.push_back(std::move(tensor_ptr));
  }
};

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &operator<<(
    std::basic_ostream<CharT, Traits> &stream, TensorNetworkV2::Origin origin) {
  switch (origin) {
    case TensorNetworkV2::Origin::Bra:
      stream << "Bra";
      break;
    case TensorNetworkV2::Origin::Ket:
      stream << "Ket";
      break;
    case TensorNetworkV2::Origin::Aux:
      stream << "Aux";
      break;
  }
  return stream;
}

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_H