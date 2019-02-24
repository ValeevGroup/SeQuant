//
// Created by Eduard Valeyev on 2019-02-02.
//

#ifndef SEQUANT2_TENSOR_NETWORK_H
#define SEQUANT2_TENSOR_NETWORK_H

#include "../../external/bliss/graph.hh"
#include "../SeQuant2/container.hpp"
#include "../SeQuant2/tensor.hpp"

namespace sequant2 {

/// @brief A (non-directed) graph view of a set of Tensor objects

/// The main role of this is to canonize itself. Since Tensors can be connected
/// by multiple Index'es (thus edges are colored), what is canonized is actually
/// the graph of indices (i.e. the dual of the tensor graph), with Tensors
/// represented by one or more vertices.
class TensorNetwork {

  /// @brief Pair of indices to Tensor terminals

  /// @note Tensor terminals in a sequence of Tensors are indexed as follows:
  /// - >0 for bra terminals (i.e. "+7" indicated connection to a bra terminal of 7th Tensor object in the sequence)
  /// - <0 for ket terminals
  /// - 0 if free (not attached to any Tensor objects)
  /// Terminal indices are sorted by the tensor index (i.e. by the absolute value of the terminal index)
  class TensorTerminalPair {
   public:
    TensorTerminalPair() = default;
    explicit TensorTerminalPair(int terminal_idx) : first_(0), second_(terminal_idx) {}
    TensorTerminalPair(int terminal_idx, const Index *idxptr) : first_(0), second_(terminal_idx), idxptr_(idxptr) {}
//    TensorTerminalPair(const TensorTerminalPair&) = default;
//    TensorTerminalPair(TensorTerminalPair&&) = default;
//    TensorTerminalPair& operator=(const TensorTerminalPair&) = default;
//    TensorTerminalPair& operator=(TensorTerminalPair&&) = default;

    TensorTerminalPair &add(int terminal_idx) {
      assert(terminal_idx != 0);  // valid idx
      if (second_ == 0) {
        second_ = terminal_idx;
      } else if (std::abs(second_) < std::abs(terminal_idx)) {
        assert(first_ == 0);  // there are slots left
        first_ = second_;
        second_ = terminal_idx;
      } else {  // put into first slot
        first_ = terminal_idx;
      }
      return *this;
    }

    bool operator<(const TensorTerminalPair &other) const {
      if (std::abs(first_) == std::abs(other.first_)) {
        return std::abs(second_) < std::abs(other.second_);
      } else {
        return std::abs(first_) < std::abs(other.first_);
      }
    }

    bool operator==(const TensorTerminalPair &other) const {
      return std::abs(first_) == std::abs(other.first_) && std::abs(second_) == std::abs(other.second_);
    }

    auto first() const { return first_; }
    auto second() const { return second_; }

    /// @return the number of attached terminals (0, 1, or 2)
    auto size() const {
      return (first_ != 0) ? 2 : ((second_ != 0) ? 1 : 0);
    }

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
  };

 public:
  template<typename ExprPtrRange>
  TensorNetwork(ExprPtrRange &exprptr_range) {
    const bool contains_a_nontensor = ranges::any_of(
        exprptr_range,
        [](const ExprPtr &exprptr) {
          return !exprptr->is<Tensor>();
        });
    if (contains_a_nontensor)
      throw std::logic_error("TensorNetwork(exprptr_range): exprptr_range contains a non-Tensor");

    auto tsrptr_range = exprptr_range | ranges::view::transform([](const ExprPtr &ex) {
      return std::static_pointer_cast<Tensor>(ex);
    });
    assert(ranges::size(tsrptr_range) == ranges::size(exprptr_range));

    for (const auto &t: tsrptr_range) {
      tensors_.push_back(t);
    }
  }

  /// @return sequence container of tensors
  const auto &tensors() const { return tensors_; }

  /// @param cardinal_tensor_labels move all tensors with these labels to the
  /// front before canonicalizing indices
  /// @param fast if true (default), does fast canonicalization that is only
  /// optimal if all tensors are distinct;
  ///        set to false to perform complete canonicalization
  ExprPtr canonicalize(
      const container::vector<std::wstring> &cardinal_tensor_labels = {},
      bool fast = true) {
    ExprPtr canon_biproduct = ex<Constant>(1);
    container::svector<TensorTerminalPair>
        idx_terminals_sorted;  // to avoid memory allocs

    if (debug_canonicalize) {
      std::wcout << "TensorNetwork::canonicalize: input tensors\n";
      size_t cnt = 0;
      ranges::for_each(tensors_, [&](const TensorPtr &t) {
        std::wcout << "tensor " << cnt++ << ": " << t->to_latex() << std::endl;
      });
    }

    if (!fast) {
      // TODO implement rigorous approach:
      // - canonize indices
      // - canonize tensors using canonical list of indices
      // Algorithm sketch:
      // - to canonize indices make a graph whose vertices are the indices as
      // well as the tensors and their terminals.
      //   - Indices with protoindices are connected to their protoindices,
      //   either directly or (if protoindices are symmetric) via a protoindex
      //   vertex.
      //   - Indices are colored by their space, which in general encodes also
      //   the space of the protoindices.
      //   - An anti/symmetric n-body tensor has 2 terminals, each connected to
      //   each other + to n index vertices.
      //   - A nonsymmetric n-body tensor has n terminals, each connected to 2
      //   indices and 1 tensor vertex which is connected to all n terminal
      //   indices.
      //   - Tensor vertices are colored by the label+rank+symmetry of the
      //   tensor; terminal vertices are colored by the color of its tensor,
      //     with the color of symm/antisymm terminals augmented by the
      //     terminal's type (bra/ket).
      // - canonize the graph using nauty

      // make the graph
      // bliss::Graph g = make_bliss_graph();
      abort();

      // canonize the graph
      bliss::Stats stats;
      //      const unsigned int* cl = g.canonical_form(stats, nullptr,
      //      nullptr);

      // reindex internal indices and reorder tensors using the computed
      // canonical order
      abort();
    } else {
      // simpler approach that will work perfectly as long as tensors are
      // distinguishable

      // - resort tensors (already done in Product::canonicalize but to make
      // this standalone do this again)
      using std::begin;
      using std::end;
      std::stable_sort(
          begin(tensors_), end(tensors_),
          [&cardinal_tensor_labels](const TensorPtr &first,
                                    const TensorPtr &second) {
            const auto cardinal_tensor_labels_end = end(cardinal_tensor_labels);
            const auto first_cardinal_it =
                std::find(begin(cardinal_tensor_labels),
                          end(cardinal_tensor_labels), first->label());
            const auto second_cardinal_it =
                std::find(begin(cardinal_tensor_labels),
                          end(cardinal_tensor_labels), second->label());
            const auto first_is_cardinal =
                first_cardinal_it != cardinal_tensor_labels_end;
            const auto second_is_cardinal =
                second_cardinal_it != cardinal_tensor_labels_end;
            if (first_is_cardinal && second_is_cardinal) {
              if (first_cardinal_it == second_cardinal_it)
                return *first < *second;
              else
                return first_cardinal_it < second_cardinal_it;
            } else if (first_is_cardinal)
              return true;
            else if (second_is_cardinal)
              return false;
            else  // neither is cardinal
              return *first < *second;
          });

      // - reindex internal indices using ordering of TensorTerminalPair as the
      // canonical definition of the internal index list
      init_indices();
      {
        auto int_idx_validator = [this](const Index &idx) {
          return this->ext_indices_.find(idx) == this->ext_indices_.end();
        };
        IndexFactory idxfac(int_idx_validator, 1);  // start reindexing from 1
        container::map<Index, Index> idxrepl;
        // resort indices_ by TensorTerminalPair ... this automatically puts
        // external indices first
        idx_terminals_sorted.resize(indices_.size());
        std::partial_sort_copy(begin(indices_), end(indices_),
                               begin(idx_terminals_sorted),
                               end(idx_terminals_sorted));

        // make index replacement list for internal indices only
        const auto num_ext_indices = ext_indices_.size();
        std::for_each(
            begin(idx_terminals_sorted) + num_ext_indices,
            end(idx_terminals_sorted),
            [&idxrepl, &idxfac](const auto &terminals) {
              const auto &idx = terminals.idx();
              if (terminals.size() == 2) {  // internal index?
                idxrepl.emplace(std::make_pair(idx, idxfac.make(idx)));
              }
            });
        if (debug_canonicalize) {
          for (const auto &idxpair : idxrepl) {
            std::wcout << "TensorNetwork::canonicalize: replacing "
                       << to_latex(idxpair.first) << " with "
                       << to_latex(idxpair.second) << std::endl;
          }
        }

        // transform indices
        const bool tag_transformed_indices =
            true;  // to replace indices when maps image and domain overlap, tag
                   // transformed indices
#ifndef NDEBUG
        // assert that tensors_ indices are not tagged if going to tag indices
        if (tag_transformed_indices) {
          for (auto &tensor : tensors_) {
            assert(ranges::none_of(
                tensor->const_braket(),
                [](const Index &idx) { return idx.tag().has_value(); }));
          }
        }
#endif
        bool pass_mutated = false;
        bool mutated = false;
        do {
          pass_mutated = false;
          for (auto &tensor : tensors_) {
            pass_mutated |=
                tensor->transform_indices(idxrepl, tag_transformed_indices);
          }
          mutated |= pass_mutated;
        } while (pass_mutated);  // transform till stops changing

        // if any replacements made, untag transformed indices as needed
        if (tag_transformed_indices) {
          for (auto &tensor : tensors_) {
            tensor->reset_tags();
          }
        }
      }
    }

    // - re-canonize tensors
    {
      DefaultTensorCanonicalizer tensor_canonizer(ext_indices_);
      for (auto &tensor : tensors_) {
        auto bp = tensor_canonizer.apply(*tensor);
        if (bp) *canon_biproduct *= *bp;
      }
    }
    indices_.clear();
    ext_indices_.clear();

    assert(canon_biproduct->is<Constant>());
    return (canon_biproduct->as<Constant>().value() == 1.) ? nullptr
                                                           : canon_biproduct;
  }

 private:
  // source tensors and indices
  container::svector<TensorPtr> tensors_;

  struct FullLabelCompare {
    using is_transparent = void;
    bool operator()(const TensorTerminalPair &first,
                    const TensorTerminalPair &second) const {
      return first.idx().full_label() < second.idx().full_label();
    }
    bool operator()(const TensorTerminalPair &first,
                    const std::wstring_view &second) const {
      return first.idx().full_label() < second;
    }
    bool operator()(const std::wstring_view &first,
                    const TensorTerminalPair &second) const {
      return first < second.idx().full_label();
    }
  };
  // Index -> TensorTerminalPair, sorted by labels
  container::set<TensorTerminalPair, FullLabelCompare> indices_;
  // ext indices do not connect tensors
  // sorted by *label* (not full label) of the corresponding value (Index)
  // this ensures that proto indices are not considered and all internal indices
  // have unique labels (not full labels)
  container::set<Index, Index::LabelCompare> ext_indices_;

  void init_indices() {
    auto idx_insert = [this](const Index &idx, int tensor_idx) {
      decltype(indices_) &indices = this->indices_;
      auto it = indices.find(idx.full_label());
      if (it == indices.end()) {
        indices.emplace(TensorTerminalPair(tensor_idx, &idx));
      } else {
        const_cast<TensorTerminalPair &>(*it).add(tensor_idx);
      }
    };

    /// this assumes antisymmetric tensors
    int t_idx = 1;
    for (const auto &t: tensors_) {
      for (const auto &idx: t->bra()) {
        idx_insert(idx, t_idx);
      }
      for (const auto &idx: t->ket()) {
        idx_insert(idx, -t_idx);
      }
      ++t_idx;
    }

    // extract external indices
    for (const auto &terminals : indices_) {
      assert(terminals.size() != 0);
      if (terminals.size() == 1) {  // external?
        auto insertion_result = ext_indices_.emplace(terminals.idx());
        assert(insertion_result.second);
      }
    }
  }
};

}  // namespace sequant2

#endif  // SEQUANT2_TENSOR_NETWORK_H
