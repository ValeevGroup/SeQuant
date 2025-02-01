//
// Created by Eduard Valeyev on 2019-02-02.
//

#ifndef SEQUANT_TENSOR_NETWORK_H
#define SEQUANT_TENSOR_NETWORK_H

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/canonicals.hpp>
#include <SeQuant/core/tensor_network/slot.hpp>
#include <SeQuant/core/tensor_network/vertex.hpp>

#include <cassert>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
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
/// @warning the terminology is a mix at best, e.g. terminal vs. slot, etc.
class TensorNetwork {
 public:
  constexpr static size_t max_rank = 256;

  // clang-format off
  /// @brief Edge in a TensorNetwork = the Index annotating it + a pair of indices to identify which Tensor terminals it's connected to

  /// @note tensor terminals in a sequence of tensors are indexed as follows:
  /// - >0 for bra terminals (i.e. "+7" indicated connection to a bra terminal
  /// of 7th tensor object in the sequence)
  /// - <0 for ket terminals
  /// - 0 if free (not attached to any tensor objects)
  /// - position records the terminal's location in the sequence of bra/ket/aux
  /// terminals (always 0 for symmetric/antisymmetric tensors) Terminal indices
  /// are sorted by the tensor index (i.e. by the absolute value of the terminal
  /// index), followed by position
  // clang-format on
  class Edge {
   public:
    struct Terminal {
      int tensor_ord = -1;
      TensorIndexSlotType slot_type = TensorIndexSlotType::Invalid;
      // index slots are grouped according to degrees of freedom and/or
      // symmetry. E.g. bra/ket slots for same particle of a nonsymmetric tensor
      // are grouped together. Bra and ket slots of a symmetric/antisymmetric
      // tensor are also grouped into their own slot groups. Each aux slot for
      // now is its own slot group. The slot groups are indexed 0, 1, ...
      int slot_group_ord = -1;

      Terminal() noexcept {};
      Terminal(int tensor_ord, TensorIndexSlotType slot_type,
               int slot_group_ord) noexcept
          : tensor_ord(tensor_ord),
            slot_type(slot_type),
            slot_group_ord(slot_group_ord) {
        assert(tensor_ord >= 0 && slot_type != TensorIndexSlotType::Invalid &&
               slot_group_ord >= 0);
      }

      friend bool operator==(const Terminal &a, const Terminal &b) {
        return a.tensor_ord == b.tensor_ord && a.slot_type == b.slot_type &&
               a.slot_group_ord == b.slot_group_ord;
      }
      friend bool operator<(const Terminal &a, const Terminal &b) {
        return std::tie(a.tensor_ord, a.slot_type, a.slot_group_ord) <
               std::tie(b.tensor_ord, b.slot_type, b.slot_group_ord);
      }

      explicit operator bool() const { return tensor_ord >= 0; }
      bool null() const { return tensor_ord < 0; }
      bool nonnull() const { return tensor_ord >= 0; }
    };

    Edge() = default;
    explicit Edge(const Terminal &t) : second_(t) {}
    Edge(const Terminal &t, const Index *idxptr)
        : second_(t), idxptr_(idxptr) {}
    //    Edge(const Edge&) = default;
    //    Edge(Edge&&) = default;
    //    Edge& operator=(const Edge&) = default;
    //    Edge& operator=(Edge&&) = default;

    Edge &connect_to(const Terminal &t) {
      assert(first_.null() || second_.null());  // not fully connected yet
      if (second_.null()) {
        assert(first_.null());  // unconnected Edge
        second_ = t;
      } else {
        // - cannot connect braket slot to aux slot
        switch (t.slot_type) {
          case TensorIndexSlotType::Aux:
            if (second_.slot_type != TensorIndexSlotType::Aux) {
              throw std::logic_error(
                  "TensorNetwork::Edge::connect_to: aux slot cannot be "
                  "connected to a non-aux slot");
            }
            break;
            // - can connect bra slot to ket slot, and vice versa
          case TensorIndexSlotType::Bra:
            if (second_.slot_type != TensorIndexSlotType::Ket) {
              throw std::logic_error(
                  "TensorNetwork::Edge::connect_to: bra slot can only be "
                  "connected to a ket slot");
            }
            break;
          case TensorIndexSlotType::Ket:
            if (second_.slot_type != TensorIndexSlotType::Bra) {
              throw std::logic_error(
                  "TensorNetwork::Edge::connect_to: ket slot can only be "
                  "connected to a bra slot");
            }
            break;
          default:
            throw std::logic_error(
                "TensorNetwork::Edge::connect_to: invalid slot");
        }

        first_ = t;
        // ensure first_ < second_
        if (second_ < first_) {
          std::swap(first_, second_);
        }
      }
      return *this;
    }

    bool operator<(const Edge &other) const {
      return std::tie(first_, second_) < std::tie(other.first_, other.second_);
    }

    friend bool operator==(const Edge &a, const Edge &b) {
      return a.first_ == b.first_ && a.second_ == b.second_;
    }

    const auto &first() const {
      assert(first_.nonnull());
      return first_;
    }
    const auto &second() const {
      assert(second_.nonnull());
      return second_;
    }
    /// access terminals by index, nonnull terminals first
    /// @param i the ordinal index, `i` must be 0 or 1
    /// @return if `i==0` return first(), if nonnull, else second(), if nonnull,
    /// else null;  if `i==1` return second(), if nonnull, else else null;
    const auto &operator[](std::size_t i) const {
      assert(i == 0 || i == 1);
      if (i == 0) {
        if (first_.nonnull())
          return first_;
        else if (second_.nonnull())
          return second_;
        else
          return null_terminal_;
      } else {  // i == 1
        if (second_.nonnull())
          return second_;
        else
          return null_terminal_;
      }
    }

    /// @return the number of attached terminals (0, 1, or 2)
    auto size() const {
      return first_.nonnull() ? 2 : (second_.nonnull() ? 1 : 0);
    }

    const Index &idx() const {
      assert(idxptr_ != nullptr);
      return *idxptr_;
    }

   private:
    // if only connected to 1 terminal, this is always null
    Terminal first_;
    // invariant: first_.tensor_order <= second_.tensor_order
    Terminal second_;
    const Index *idxptr_ = nullptr;

    static inline Terminal null_terminal_ = {};
  };

  static inline auto edge2index_ = [](const Edge &e) -> const Index & {
    return e.idx();
  };

 public:
  /// @throw std::logic_error if exprptr_range contains a non-tensor
  /// @note uses RTTI
  template <typename ExprPtrRange>
  TensorNetwork(ExprPtrRange &exprptr_range) {
    if (exprptr_range.size() > 0) {
      for (auto &&ex : exprptr_range) {
        auto t = std::dynamic_pointer_cast<AbstractTensor>(ex);
        std::size_t count = 0;
        if (t) {
          tensors_.emplace_back(t);
          tensor_input_ordinals_.emplace_back(count++);
        } else {
          throw std::logic_error(
              "TensorNetwork::TensorNetwork: non-tensors in the given "
              "expression "
              "range");
        }
      }
      return;
    } else {
      if constexpr (Expr::is_shared_ptr_of_expr<ExprPtrRange>::value) {
        if (auto tensor =
                std::dynamic_pointer_cast<AbstractTensor>(exprptr_range)) {
          tensors_.emplace_back(tensor);
          tensor_input_ordinals_.emplace_back(0);
          return;
        }
      }
    }
    throw std::logic_error(
        "TensorNetwork::TensorNetwork: non-tensors in the given expression "
        "range");
  }

  /// @return const reference to the sequence of tensors
  /// @note after invoking TensorNetwork::canonicalize() the order of
  /// tensors may be different from that provided as input; use
  /// tensor_input_ordinals() to obtain the input ordinals of
  /// the tensors in the result
  const auto &tensors() const { return tensors_; }

  const auto &tensor_input_ordinals() const { return tensor_input_ordinals_; }

  using named_indices_t = container::set<Index, Index::FullLabelCompare>;

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

  /// metadata produced by canonicalize_slots()
  struct SlotCanonicalizationMetadata {
    /// list of named indices
    named_indices_t named_indices;

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
    container::svector<named_indices_t::const_iterator> named_indices_canonical;

    /// canonicalized colored graph, use graph->cmp to compare against another
    /// to detect equivalence
    std::shared_ptr<bliss::Graph> graph;
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
      const named_indices_t *named_indices = nullptr,
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

 private:
  /// list of tensors
  /// - before canonicalize(): input
  /// - after canonicalize(): canonical
  container::svector<AbstractTensorPtr> tensors_;
  /// input ordinals of tensors in tensors_
  container::svector<std::size_t> tensor_input_ordinals_;

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
  using edges_t = container::set<Edge, FullLabelCompare>;
  mutable edges_t edges_;
  // set to true by init_edges();
  mutable bool have_edges_ = false;
  // ext indices do not connect tensors
  // sorted by *label* (not full label) of the corresponding value (Index)
  // this ensures that proto indices are not considered and all internal indices
  // have unique labels (not full labels)
  mutable named_indices_t ext_indices_;
  /// some proto indices may not be in edges_ if they appear exclusively among
  /// proto indices
  /// @note these will need to be processed separately from the rest
  /// to appear as vertices on the graph
  mutable named_indices_t pure_proto_indices_;
  /// grand list of all indices is view of concatenated ranges of indices in
  /// edges_ and pure_proto_indices_
  mutable ranges::concat_view<
      ranges::transform_view<ranges::ref_view<edges_t>, decltype(edge2index_)>,
      ranges::ref_view<named_indices_t>>
      grand_index_list_;

  // replacements of anonymous indices produced by the last call to
  // canonicalize()
  container::map<Index, Index> idxrepl_;

  /// initializes edges_, ext_indices_, and pure_proto_indices_
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

  /// accessor for the list of anonymous index replacements performed by the
  /// last call to canonicalize()
  /// @return replacements of anonymous indices performed by the last call to
  /// canonicalize()
  const auto &idxrepl() const { return idxrepl_; };

 public:
  /// @brief converts the network into a Bliss graph whose vertices are indices
  /// and tensor vertex representations
  /// @param[in] named_indices pointer to the set of named indices (ordinarily,
  /// this includes all external indices);
  ///            default is nullptr, which means use all external indices for
  ///            named indices
  /// @param[in] distinct_named_indices if false, will use same color for all
  /// named indices that have same Index::color(), else will use distinct color
  /// for each
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
  make_bliss_graph(const named_indices_t *named_indices = nullptr,
                   bool distinct_named_indices = true) const;
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_H
