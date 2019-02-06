//
// Created by Eduard Valeyev on 2019-02-02.
//

#ifndef SEQUANT2_TENSOR_NETWORK_H
#define SEQUANT2_TENSOR_NETWORK_H

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
      } else if (std::abs(second_) > std::abs(terminal_idx)) {
        assert(first_ != 0);  // there are slots left
        first_ = second_;
        second_ = terminal_idx;
      } else { // put into first slot
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

  ExprPtr canonicalize() {
    ExprPtr canon_biproduct = make<Constant>(1);
    container::svector<TensorTerminalPair> idx_terminals_sorted;  // to avoid memory allocs

    // TODO implement rigorous approach:
    // - canonize indices
    // - canonize tensors using canonical list of indices

    // simpler approach that will work perfectly as long as tensors are distinguishable
    // - resort tensors (already done in Product::canonicalize but to make this standalone do this again)
    using std::begin;
    using std::end;
    std::stable_sort(begin(tensors_), end(tensors_), [](const TensorPtr &first, const TensorPtr &second) {
      return *first < *second;
    });
    // - reindex internal indices using ordering of TensorTerminalPair as the canonical definition of the internal index list
    init_indices();
    {
      IndexFactory idxfac;
      std::map<Index, Index> idxrepl;
      // resort indices_ by TensorTerminalPair ... this automatically puts external indices first
      idx_terminals_sorted.resize(indices_.size());
      std::partial_sort_copy(begin(indices_), end(indices_), begin(idx_terminals_sorted), end(idx_terminals_sorted));

      // make index replacement list for internal indices only
      const auto num_ext_indices = ext_indices_.size();
      std::for_each(begin(idx_terminals_sorted) + num_ext_indices,
                    end(idx_terminals_sorted),
                    [&idxrepl, &idxfac](const auto &terminals) {
                      const auto &idx = terminals.idx();
                      if (terminals.size() == 2) {  // internal index?
                        idxrepl.emplace(std::make_pair(idx, idxfac.make(idx.space())));
                      }
                    });
      if (debug_canonicalize) {
        for (const auto &idxpair: idxrepl) {
          std::wcout << "TensorNetwork::canonicalize: replacing " << idxpair.first.label() << " with "
                     << idxpair.second.label() << std::endl;
        }
      }
      for (auto &tensor: tensors_) {
        tensor->transform_indices(idxrepl);
      }
    }
    // - re-canonize tensors
    {
      DefaultTensorCanonicalizer tensor_canonizer(ext_indices_);
      for (auto &tensor: tensors_) {
        auto bp = tensor_canonizer.apply(*tensor);
        if (bp) *canon_biproduct *= *bp;
      }
    }

    indices_.clear();
    ext_indices_.clear();

    assert(canon_biproduct->is<Constant>());
    return (canon_biproduct->as<Constant>().value() == 1.) ? nullptr : canon_biproduct;
  }

 private:
  // source tensors and indices
  container::svector<TensorPtr> tensors_;

  struct LabelComparer {
    using is_transparent = void;
    bool operator()(const TensorTerminalPair &first, const TensorTerminalPair &second) const {
      return first.idx().label() < second.idx().label();
    }
    bool operator()(const TensorTerminalPair &first, const std::wstring_view &second) const {
      return first.idx().label() < second;
    }
    bool operator()(const std::wstring_view &first, const TensorTerminalPair &second) const {
      return first < second.idx().label();
    }
  };
  // Index -> TensorTerminalPair, sorted by labels
  std::set<TensorTerminalPair, LabelComparer> indices_;
  // ext indices do no connect tensors
  std::set<Index> ext_indices_;

  void init_indices() {
    auto idx_insert = [this](const Index &idx, int tensor_idx) {
      decltype(indices_) &indices = this->indices_;
      auto it = indices.find(idx.label());
      if (it == indices.end()) {
        indices.emplace(TensorTerminalPair(tensor_idx, &idx));
      } else {
        const_cast<TensorTerminalPair &>(*it).add(tensor_idx);
      }
    };

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
    for (const auto &terminals: indices_) {
      assert(terminals.size() != 0);
      if (terminals.size() == 1) { // external?
        ext_indices_.emplace(terminals.idx());
      }
    }

  }
};

}  // namespace sequant2

#endif  // SEQUANT2_TENSOR_NETWORK_H
