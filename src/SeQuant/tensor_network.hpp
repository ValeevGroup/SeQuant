//
// Created by Eduard Valeyev on 2019-02-02.
//

#ifndef SEQUANT_TENSOR_NETWORK_H
#define SEQUANT_TENSOR_NETWORK_H

#include "../SeQuant/container.hpp"
#include "../SeQuant/tensor.hpp"

// forward declarations
namespace bliss {
class Graph;
}

namespace sequant {

/// @brief A (non-directed) graph view of a set of Tensor (or Tensor like) objects

/// @tparam Tensor_ Tensor or another type for which @c bra(t) and @c ket(t)
///         are valid expressions and evaluate to a range of Index objects.
/// @note The main role of this is to canonize itself. Since Tensor_ objects can be connected
/// by multiple Index'es (thus edges are colored), what is canonized is actually
/// the graph of indices (i.e. the dual of the tensor graph), with Tensor_ objects
/// represented by one or more vertices.
template <typename Tensor_ = Tensor>
class TensorNetwork {
 public:

  /// @brief Edge in a TensorNetwork = Index annotating it + pair of indices to identify which Tensor_ terminals it's connected to

  /// @note tensor terminals in a sequence of tensors are indexed as follows:
  /// - >0 for bra terminals (i.e. "+7" indicated connection to a bra terminal
  /// of 7th tensor object in the sequence)
  /// - <0 for ket terminals
  /// - 0 if free (not attached to any tensor objects)
  /// - position records the terminals location in the sequence of bra/ket
  /// terminals (always 0 for symmetric/antisymmetric tensors) Terminal indices
  /// are sorted by the tensor index (i.e. by the absolute value of the terminal
  /// index), followed by position
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
      assert(terminal_idx != 0);  // valid idx
      if (second_ == 0) {
        second_ = terminal_idx;
        second_position_ = position;
      } else if (std::abs(second_) < std::abs(terminal_idx)) {
        assert(first_ == 0);  // there are slots left
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
    /// this assumes antisymmetric tensors ... for nonsymmetric tensors need to
    /// know not only tensor index but also terminal index within bra/ket
    int first_ = 0;
    int second_ = 0;
    const Index *idxptr_ = nullptr;
    int first_position_ = 0;
    int second_position_ = 0;
  };

  enum class VertexType {
    Index,
    SPBundle,
    TensorBra,
    TensorKet,
    TensorBraKet,
    TensorCore
  };

 public:
  template <typename ExprPtrRange>
  TensorNetwork(ExprPtrRange &exprptr_range) {
    const bool contains_a_nontensor = ranges::any_of(
        exprptr_range,
        [](const ExprPtr &exprptr) { return !exprptr->is<Tensor_>(); });
    if (contains_a_nontensor)
      throw std::logic_error(
          "TensorNetwork(exprptr_range): exprptr_range contains a non-Tensor");

    auto tsrptr_range =
        exprptr_range | ranges::view::transform([](const ExprPtr &ex) {
          return std::static_pointer_cast<Tensor_>(ex);
        });
    assert(ranges::size(tsrptr_range) == ranges::size(exprptr_range));

    for (auto&&t : tsrptr_range) {
      tensors_.push_back(t);
    }
  }

  /// @return sequence container of tensors
  const auto &tensors() const { return tensors_; }

  /// @param cardinal_tensor_labels move all tensors with these labels to the
  /// front before canonicalizing indices
  /// @param fast if true (default), does fast canonicalization that is only
  /// optimal if all tensors are distinct; set to false to perform complete
  /// canonicalization
  ExprPtr canonicalize(
      const container::vector<std::wstring> &cardinal_tensor_labels = {},
      bool fast = true);

 private:
  using Tensor_Ptr = std::shared_ptr<Tensor_>;

  // source tensors and indices
  container::svector<Tensor_Ptr> tensors_;

  struct FullLabelCompare {
    using is_transparent = void;
    bool operator()(const Edge &first,
                    const Edge &second) const {
      return first.idx().full_label() < second.idx().full_label();
    }
    bool operator()(const Edge &first,
                    const std::wstring_view &second) const {
      return first.idx().full_label() < second;
    }
    bool operator()(const std::wstring_view &first,
                    const Edge &second) const {
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
  mutable container::set<Index, Index::LabelCompare> ext_indices_;

  /// initializes edges_ and ext_indices_
  void init_edges() const;

 public:
  /// accessor for Edge objects
  /// @return a set of Edge objects, sorted by their Index's full label
  const auto& edges() const {
    init_edges();
    return edges_;
  }

 private:
  /// @brief converts the network into a graph whose vertices are indices and
  /// tensor vertex representations
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
  make_bliss_graph() const;
};

}  // namespace sequant

#include "tensor_network.impl.hpp"

#endif  // SEQUANT_TENSOR_NETWORK_H
