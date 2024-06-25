#ifndef SEQUANT_WICK_GRAPH_H
#define SEQUANT_WICK_GRAPH_H

#include "abstract_tensor.hpp"
#include "container.hpp"
#include "tensor.hpp"
#include "vertex_type.hpp"

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
class WickGraph {
 public:
  constexpr static size_t max_rank = 256;

  // clang-format off
  /// @brief Edge in a TensorNetwork = the Index annotating it + a pair of indices to identify which Tensor terminals it's connected to

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
    explicit Edge(int terminal_idx, int position = 0)
        : first_(0), second_(terminal_idx), second_position_(position) {}
    Edge(int terminal_idx, const Index *idxptr, int position = 0)
        : first_(0),
          second_(terminal_idx),
          idxptr_(idxptr),
          second_position_(position) {}
    //    Edge(const Edge&) = default;
    //    Edge(Edge&&) = default;
    //    Edge& operator=(const Edge&) = default;
    //    Edge& operator=(Edge&&) = default;

    Edge &connect_to(int terminal_idx, int position = 0) {
      assert(first_ == 0 || second_ == 0);  // not connected yet
      assert(terminal_idx != 0);            // valid idx
      if (second_ == 0) {                   // unconnected Edge
        second_ = terminal_idx;
        second_position_ = position;
      } else if (std::abs(second_) <
                 std::abs(terminal_idx)) {  // connected to 2 Edges? ensure
                                            // first_ < second_
        assert(first_ == 0);                // there are slots left
        first_ = second_;
        first_position_ = second_position_;
        second_ = terminal_idx;
        second_position_ = position;
      } else {  // put into first slot
        first_ = terminal_idx;
        first_position_ = position;
      }
      return *this;
    }

    bool operator<(const Edge &other) const {
      if (std::abs(first_) == std::abs(other.first_)) {
        if (first_position_ == other.first_position_) {
          if (std::abs(second_) == std::abs(other.second_)) {
            return second_position_ < other.second_position_;
          } else {
            return std::abs(second_) < std::abs(other.second_);
          }
        } else {
          return first_position_ < other.first_position_;
        }
      } else {
        return std::abs(first_) < std::abs(other.first_);
      }
    }

    bool operator==(const Edge &other) const {
      return std::abs(first_) == std::abs(other.first_) &&
             std::abs(second_) == std::abs(other.second_) &&
             first_position_ == other.first_position_ &&
             second_position_ == other.second_position_;
    }

    auto first() const { return first_; }
    auto second() const { return second_; }
    auto first_position() const { return first_position_; }
    auto second_position() const { return second_position_; }

    /// @return the number of attached terminals (0, 1, or 2)
    auto size() const { return (first_ != 0) ? 2 : ((second_ != 0) ? 1 : 0); }

    const Index &idx() const {
      assert(idxptr_ != nullptr);
      return *idxptr_;
    }

   private:
    // if only connected to 1 terminal, this is always 0
    // otherwise first_ <= second_
    int first_ = 0;
    int second_ = 0;
    const Index *idxptr_ = nullptr;
    int first_position_ = 0;
    int second_position_ = 0;
  };

 public:
  /// @throw std::logic_error if exprptr_range contains a non-tensor
  /// @note uses RTTI
  template <typename ExprPtrRange>
  WickGraph(ExprPtrRange &exprptr_range) {
    for (auto &&ex : exprptr_range) {
      auto t = std::dynamic_pointer_cast<AbstractTensor>(ex);
      if (t) {
        tensors_.emplace_back(t);
      } else {
        throw std::logic_error(
            "TensorNetwork::TensorNetwork: non-tensors in the given expression "
            "range");
      }
    }
  }

  /// @return const reference to the sequence of tensors
  /// @note the order of tensors may be different from that provided as input
  const auto &tensors() const { return tensors_; }

  using named_indices_t = container::set<Index, Index::LabelCompare>;

 private:
  // source tensors and indices
  container::svector<AbstractTensorPtr> tensors_;

  struct FullLabelCompare {
    using is_transparent = void;
    bool operator()(const Edge &first, const Edge &second) const {
      return first.idx().full_label() < second.idx().full_label();
    }
    bool operator()(const Edge &first, const std::wstring_view &second) const {
      return first.idx().full_label() < second;
    }
    bool operator()(const std::wstring_view &first, const Edge &second) const {
      return first < second.idx().full_label();
    }
  };
  // Index -> Edge, sorted by full label
  mutable container::set<Edge, FullLabelCompare> edges_;
  // set to true by init_edges();
  mutable bool have_edges_ = false;
  // ext indices do not connect tensors
  // sorted by *label* (not full label) of the corresponding value (Index)
  // this ensures that proto indices are not considered and all internal indices
  // have unique labels (not full labels)
  mutable named_indices_t ext_indices_;

  // replacements of anonymous indices produced by the last call to
  // canonicalize()
  container::map<Index, Index> idxrepl_;

  /// initializes edges_ and ext_indices_
  void init_edges() const;

 public:
  /// accessor for the Edge object sequence
  /// @return const reference to the sequence container of Edge objects, sorted
  /// by their Index's full label
  /// @sa Edge
  const auto &edges() const {
    init_edges();
    return edges_;
  }

  /// @brief Returns a range of external indices, i.e. those indices that do not
  /// connect tensors

  /// @note The external indices are sorted by *label* (not full label) of the
  /// corresponding value (Index)
  const auto &ext_indices() const {
    if (edges_.empty()) init_edges();
    return ext_indices_;
  }

 public:
  /// @brief converts the network into a Bliss graph whose vertices are indices
  /// and tensor vertex representations
  /// @param[in] named_indices pointer to the set of named indices (ordinarily,
  /// this includes all external indices);
  ///            default is nullptr, which means use all external indices for
  ///            named indices
  /// @return {shared_ptr to Graph, vector of vertex labels, vector of vertex
  /// colors, vector of vertex types}

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
  std::tuple<std::shared_ptr<bliss::Graph>, std::vector<std::wstring>,
             std::vector<std::size_t>, std::vector<VertexType>>
  make_bliss_graph(const named_indices_t *named_indices = nullptr) const;
};

}  // namespace sequant

#endif  // SEQUANT_WICK_GRAPH_H
