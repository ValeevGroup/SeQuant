//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT_WICK_HPP
#define SEQUANT_WICK_HPP

#include <bitset>
#include <mutex>
#include <utility>

#include "op.hpp"
#include "ranges.hpp"
#include "runtime.hpp"
#include "tensor.hpp"

namespace sequant {

/// Applies Wick's theorem to a sequence of normal-ordered operators.
///
/// @tparam S particle statistics
template <Statistics S>
class WickTheorem {
 public:
  // see
  // https://softwareengineering.stackexchange.com/questions/257705/unit-test-private-method-in-c-using-a-friend-class?answertab=votes#tab-top
  template <class T>
  struct access_by;
  template <class T>
  friend struct access_by;

  static constexpr const Statistics statistics = S;

  WickTheorem(WickTheorem &&) = default;
  WickTheorem &operator=(WickTheorem &&) = default;

 private:
  // make copying private to be able to adjust state
  WickTheorem(const WickTheorem &) = default;
  WickTheorem &operator=(const WickTheorem &) = default;

 public:
  explicit WickTheorem(const NormalOperatorSequence<S> &input) : input_(input) {
    assert(input.size() <= max_input_size);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
    if constexpr (statistics == Statistics::BoseEinstein) {
      assert(input.empty() || input.vacuum() == Vacuum::Physical);
    }
  }

  explicit WickTheorem(ExprPtr expr_input) : expr_input_(expr_input) {}

  /// constructs WickTheorem from @c other with expression input set to @c
  /// expr_input
  WickTheorem(ExprPtr expr_input, const WickTheorem &other)
      : WickTheorem(other) {
    // copy ctor does not do anything useful, so this is OK
    expr_input_ = expr_input;
    reset_stats();
  }

  /// Controls whether next call to compute() will full contractions only or all
  /// (including partial) contractions. By default compute() generates full
  /// contractions only.
  /// @param sf if false, will evaluate all contractions.
  /// @return reference to @c *this , for daisy-chaining
  WickTheorem &full_contractions(bool fc) {
    full_contractions_ = fc;
    return *this;
  }
  /// Controls whether next call to compute() will assume spin-free or
  /// spin-orbital normal-ordered operators By default compute() assumes
  /// spin-orbital operators.
  /// @param sf if true, will complete full contractions only.
  WickTheorem &spinfree(bool sf) {
    spinfree_ = sf;
    return *this;
  }
  /// Controls whether:
  /// - Op's of the same type within each NormalOperator are
  /// assumed topologically equivalent, and
  /// - (<b>if an Expr is given as input</b>) NormalOperator objects attached to
  /// Tensor objects are
  ///   considered equivalent. (If a NormalOperatorSequence is given as input,
  ///   such use of topological equivalence if enabled by invoking
  ///   set_op_partitions() ).
  ///
  /// This is useful to to eliminate the topologically-equivalent contractions
  /// when fully-contracted result (i.e. the vacuum average) is sought.
  /// By default the use of topology is not enabled.
  /// @param sf if true, will utilize the topology to minimize work.
  /// @warning currently is only supported if full contractions are requested
  /// @sa set_op_partitions()
  WickTheorem &use_topology(bool ut) {
    assert(full_contractions_);
    use_topology_ = ut;
    return *this;
  }

  /// Specifies the external indices; by default assume all indices are summed
  /// over
  /// @param ext_inds external (nonsummed) indices
  template <typename IndexContainer>
  WickTheorem &set_external_indices(IndexContainer &&external_indices) {
    if constexpr (std::is_convertible_v<IndexContainer,
                                        decltype(external_indices_)>)
      external_indices_ = std::forward<IndexContainer>(external_indices);
    else {
      ranges::for_each(std::forward<IndexContainer>(external_indices),
                       [this](auto &&v) {
                         auto result = this->external_indices_.emplace(v);
                         assert(result.second);
                       });
    }
    return *this;
  }

  /// @name operator connectivity specifiers
  ///
  /// Ensures that the given pairs of normal operators are connected; by default
  /// will not constrain connectivity
  /// @param op_index_pairs the list of pairs of op indices to be connected in
  /// the result
  ///
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  ///@{

  /// @tparam IndexPairContainer a sequence of std::pair<Integer,Integer>
  template <typename IndexPairContainer>
  WickTheorem &set_nop_connections(IndexPairContainer &&op_index_pairs) {
    if (expr_input_ == nullptr || !nop_connections_input_.empty()) {
      for (const auto &opidx_pair : op_index_pairs) {
        if (opidx_pair.first < 0 || opidx_pair.first >= input_.size()) {
          throw std::invalid_argument(
              "WickTheorem::set_nop_connections: nop index out of range");
        }
        if (opidx_pair.second < 0 || opidx_pair.second >= input_.size()) {
          throw std::invalid_argument(
              "WickTheorem::set_nop_connections: nop index out of range");
        }
      }
      if (op_index_pairs.size() != 0) {
        nop_connections_.resize(input_.size());
        for (auto &v : nop_connections_) {
          v.set();
        }
        for (const auto &opidx_pair : op_index_pairs) {
          nop_connections_[opidx_pair.first].reset(opidx_pair.second);
          nop_connections_[opidx_pair.second].reset(opidx_pair.first);
        }
      }
      nop_connections_input_.clear();
    } else {
      ranges::for_each(op_index_pairs, [this](const auto &idxpair) {
        nop_connections_input_.push_back(idxpair);
      });
    }

    return *this;
  }

  /// @tparam Integer an integral type
  template <typename Integer = long>
  WickTheorem &set_nop_connections(
      std::initializer_list<std::pair<Integer, Integer>> op_index_pairs) {
    return this->set_nop_connections<const decltype(op_index_pairs) &>(
        op_index_pairs);
  }
  ///@}

  /// @name topological partition specifiers
  ///
  /// Specifies topological partition of normal operators; free (non-connected)
  /// operators in the same partition are considered topologically equivalent,
  /// hence if only full contractions are needed only contractions to the first
  /// available operator in a partition is needed (multiplied by the degeneracy)
  /// @param op_partitions list of operator partitions
  /// @note if this partitions are not given, every operator is assumed to be in
  /// its own partition
  /// @internal this performs only the first phase of initialization of
  /// nop_topological_partition_
  ///           since the number of operators is not guaranteed to be known
  ///           until compute() is called (e.g. if op product is given as an
  ///           Expr ). The second initialization phase (i.e. resizing to make
  ///           sure every operator is assigned to a partition, including those
  ///           not mentioned in @c op_partitions) will be completed in
  ///           compute_nopseq().
  ///
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  ///@{

  /// @tparam IndexListContainer a sequence of sequences of Integer types
  template <typename IndexListContainer>
  WickTheorem &set_op_partitions(IndexListContainer &&op_partitions) {
    using std::size;
    size_t partition_cnt = 0;
    auto current_nops = size(nop_topological_partition_);
    for (auto &&partition : op_partitions) {
      for (auto &&op_idx : partition) {
        assert(op_idx >= 0);
        if (op_idx >= current_nops) {
          current_nops = upsize_nop_topological_partition(op_idx + 1);
        }
        nop_topological_partition_[op_idx] = partition_cnt + 1;
      }
      ++partition_cnt;
    }
    return *this;
  }

  /// @tparam Integer an integral type
  template <typename Integer = long>
  WickTheorem &set_op_partitions(
      std::initializer_list<std::initializer_list<Integer>> op_partitions) {
    return this->set_op_partitions<const decltype(op_partitions) &>(
        op_partitions);
  }
  ///@}

  /// Computes and returns the result
  /// @param count_only if true, will return the total number of terms, as a
  /// Constant.
  /// @return the result of applying Wick's theorem; either a Constant, a
  /// Product, or a Sum
  /// @warning this is not reentrant, but is optionally threaded internally
  ExprPtr compute(const bool count_only = false);

  /// Collects compute statistics
  class Stats {
   public:
    Stats() : num_attempted_contractions(0), num_useful_contractions(0) {}
    Stats(const Stats &other) noexcept {
      num_attempted_contractions.store(other.num_attempted_contractions.load());
      num_useful_contractions.store(other.num_useful_contractions.load());
    }
    Stats &operator=(const Stats &other) noexcept {
      num_attempted_contractions.store(other.num_attempted_contractions.load());
      num_useful_contractions.store(other.num_useful_contractions.load());
      return *this;
    }

    void reset() {
      num_attempted_contractions = 0;
      num_useful_contractions = 0;
    }

    Stats &operator+=(const Stats &other) {
      num_attempted_contractions += other.num_attempted_contractions;
      num_useful_contractions += other.num_useful_contractions;
      return *this;
    }

    std::atomic<size_t> num_attempted_contractions;
    std::atomic<size_t> num_useful_contractions;
  };

  /// Statistics accessor
  /// @return statistics of compute calls since creation, or since the last call
  /// to reset_stats()
  const Stats &stats() const { return stats_; }

  /// Statistics accessor
  /// @return statistics of compute calls since creation, or since the last call
  /// to reset_stats()
  Stats &stats() { return stats_; }

  /// Statistics reset
  void reset_stats() { stats_.reset(); }

 private:
  static constexpr size_t max_input_size =
      32;  // max # of operators in the input sequence

  // if nonnull, apply wick to the whole expression recursively, else input_ is
  // set this is mutated by compute
  mutable ExprPtr expr_input_;

  mutable NormalOperatorSequence<S> input_;
  bool full_contractions_ = true;
  bool spinfree_ = false;
  bool use_topology_ = false;
  mutable Stats stats_;

  container::set<Index> external_indices_;
  // for each operator specifies the reverse bitmask of target connections
  // (0 = must connect)
  container::svector<std::bitset<max_input_size>> nop_connections_;
  container::vector<std::pair<size_t, size_t>>
      nop_connections_input_;  // only used to cache input to
                               // set_nop_connections_

  // for each operator specifies its topological partition (0 = topologically
  // unique)
  mutable container::svector<size_t> nop_topological_partition_;

  /// upsizes nop_topological_partition_, filling new entries with zeroes
  /// noop if current size > new_size
  /// @return the (updated) size of nop_topological_partition_
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  size_t upsize_nop_topological_partition(size_t new_size) const {
    using std::size;
    const auto current_size = size(nop_topological_partition_);
    if (new_size > current_size) {
      nop_topological_partition_.resize(new_size);
      for (size_t i = current_size; i != new_size; ++i)
        nop_topological_partition_[i] = 0;
      return new_size;
    } else
      return current_size;
  }

  /// Evaluates wick_ theorem for a single NormalOperatorSequence
  /// @return the result of applying Wick's theorem
  ExprPtr compute_nopseq(const bool count_only) const {
    if (spinfree_)
      throw std::logic_error(
          "WickTheorem::compute: spinfree=true not yet supported");
    // process cached nop_connections_input_, if needed
    if (!nop_connections_input_.empty())
      const_cast<WickTheorem<S> &>(*this).set_nop_connections(
          nop_connections_input_);
    // size nop_topological_partition_ to match input_, if needed
    upsize_nop_topological_partition(input_.size());
    // now compute
    auto result = compute_nontensor_wick(count_only);
    return std::move(result);
  }

  /// carries state down the stack of recursive calls
  struct NontensorWickState {
    NontensorWickState(
        const NormalOperatorSequence<S> &opseq,
        /// TODO rename op -> nop to distinguish Op and NormalOperator
        const container::svector<size_t> &op_toppart)
        : opseq(opseq),
          opseq_size(opseq.opsize()),
          level(0),
          left_op_offset(0),
          count_only(false),
          count(0),
          nop_connections(opseq.size()),
          nop_adjacency_matrix(opseq.size() * (opseq.size() - 1) / 2, 0),
          nop_nconnections(opseq.size(), 0),
          nop_topological_partition(op_toppart) {
      init_topological_partitions();
      init_input_index_columns();
    }

    NontensorWickState(const NontensorWickState &) = delete;
    NontensorWickState(NontensorWickState &&) = delete;
    NontensorWickState &operator=(const NontensorWickState &) = delete;
    NontensorWickState &operator=(NontensorWickState &&) = delete;

    NormalOperatorSequence<S> opseq;  //!< current state of operator sequence
    std::size_t opseq_size;           //!< current size of opseq
    Product sp;                       //!< current prefactor
    int level;                        //!< level in recursive wick call stack
    size_t left_op_offset;  //!< where to start looking for contractions
    bool count_only;  //!< if true, only track the total number of summands in
                      //!< the result (i.e. 1 (the normal product) + the number
                      //!< of contractions (if normal wick result is wanted) or
                      //!< the number of complete constractions (if want
                      //!< complete contractions only)
    std::atomic<size_t> count;  //!< if count_only is true, will count the
                                //!< total number of terms
    container::svector<std::bitset<max_input_size>>
        nop_connections;  //!< bitmask of connections for each nop (1 =
                          //!< connected)
    container::svector<size_t>
        nop_adjacency_matrix;  //!< number of connections between each nop, only
                               //!< lower triangle is kept

    container::svector<std::pair<Index, Index>>
        input_partner_indices;  //!< list of {cre,ann} pairs of Index objects in
                                //!< the input whose corresponding Op<S> objects
                                //!< act on the same particle

    /// "merges" partner index pair from input_index_columns with contracted
    /// Index pairs in this->sp
    auto make_target_partner_indices() const {
      // copy all pairs in the input product
      container::svector<std::pair<Index, Index>> result(input_partner_indices);
      // for every contraction so far encountered ...
      for (auto &contr : sp) {
        // N.B. sp is composed of 1-particle contractions
        assert(contr->template is<Tensor>() &&
               contr->template as<Tensor>().rank() == 1 &&
               contr->template as<Tensor>().label() ==
                   sequant::overlap_label());
        const auto &contr_t = contr->template as<Tensor>();
        const auto &bra_idx = contr_t.bra().at(0);
        const auto &ket_idx = contr_t.ket().at(0);
        // ... if both bra and ket indices were in the input list, "merge" their
        // pairs
        auto bra_it = ranges::find_if(result, [&bra_idx](const auto &p) {
          assert(p.first != bra_idx);
          return p.second == bra_idx;
        });
        if (bra_it != result.end()) {
          auto ket_it = ranges::find_if(result, [&ket_idx](const auto &p) {
            assert(p.second != ket_idx);
            return p.first == ket_idx;
          });
          if (ket_it != result.end()) {
            assert(bra_it != ket_it);
            if (ket_it > bra_it) {
              bra_it->second = std::move(ket_it->second);
              result.erase(ket_it);
            } else {
              ket_it->first = std::move(bra_it->first);
              result.erase(bra_it);
            }
          }
        }
      }
      return result;
    }

    /// for each operator specifies how many connections it currently has
    container::svector<size_t> nop_nconnections;
    /// maps op to its topological partition index (1-based, 0 = no partition)
    container::svector<size_t> nop_topological_partition;
    /// current state of partitions (will only match nop_topological_partition
    /// before any contractions have occurred)
    /// - when an operator is connected it's removed from the partition
    /// - when it is disconnected fully it's re-added to the partition
    container::vector<container::set<size_t>> topological_partitions;

    // populates partitions using the data from nop_topological_partition
    void init_topological_partitions() {
      // partition indices in nop_topological_partition are 1-based
      const auto npartitions = *ranges::max_element(nop_topological_partition);
      topological_partitions.resize(npartitions);
      size_t op_cnt = 0;
      ranges::for_each(
          nop_topological_partition, [this, &op_cnt](size_t toppart_idx) {
            if (toppart_idx > 0) {  // in a partition
              topological_partitions.at(toppart_idx - 1).insert(op_cnt);
            }
            ++op_cnt;
          });
      // assert that we don't have empty partitions due to invalid contents of
      // nop_topological_partition
      assert(ranges::any_of(topological_partitions, [](auto &&partition) {
               return partition.empty();
             }) == false);
    }

    // populates target_particle_ops
    void init_input_index_columns() {
      // for each NormalOperator
      for (auto &nop : opseq) {
        using ranges::views::reverse;
        using ranges::views::zip;

        // zip in reverse order to handle non-number-conserving ops (i.e.
        // nop.creators().size() != nop.annihilators().size()) reverse after
        // insertion to restore canonical partner index order (in the order of
        // particle indices in the normal operators) 7/18/2022 N.B. reverse(zip)
        // for some reason is broken, hence the ugliness
        std::size_t ninserted = 0;
        for (auto &&[cre, ann] :
             zip(reverse(nop.creators()), reverse(nop.annihilators()))) {
          input_partner_indices.emplace_back(cre.index(), ann.index());
          ++ninserted;
        }
        std::reverse(input_partner_indices.rbegin(),
                     input_partner_indices.rbegin() + ninserted);
      }
    }

    template <typename T>
    static auto lowtri_idx(T i, T j) {
      assert(i != j);
      auto ii = std::max(i, j);
      auto jj = std::min(i, j);
      return ii * (ii - 1) / 2 + jj;
    }

    /// @brief Updates connectivity if contraction satisfies target connectivity

    /// If the target connectivity will be violated by this contraction, keep
    /// the state unchanged and return false
    template <typename Cursor>
    inline bool connect(const container::svector<std::bitset<max_input_size>>
                            &target_nop_connections,
                        const Cursor &op1_cursor, const Cursor &op2_cursor) {
      auto update_topology = [this](size_t op_idx) {
        const auto nconnections = nop_nconnections[op_idx];
        // if using topological partitions for normal ops, and this operator is
        // in one of them, remove it on first connection
        if (!topological_partitions.empty()) {
          auto partition_idx = nop_topological_partition[op_idx];
          if (nconnections == 0 && partition_idx > 0) {
            --partition_idx;  // to 0-based
            assert(topological_partitions.at(partition_idx).size() > 0);
            auto removed = topological_partitions[partition_idx].erase(op_idx);
            assert(removed);
          }
        }
        ++nop_nconnections[op_idx];
      };

      // local vars
      const auto op1_idx = op1_cursor.range_ordinal();
      const auto op2_idx = op2_cursor.range_ordinal();
      if (target_nop_connections
              .empty()) {  // if no constraints, all is fair game
        update_topology(op1_idx);
        update_topology(op2_idx);
        return true;
      }
      const auto op1_op2_connected = nop_connections[op1_idx].test(op2_idx);

      // update the connectivity
      if (!op1_op2_connected) {
        nop_connections[op1_idx].set(op2_idx);
        nop_connections[op2_idx].set(op1_idx);
      }

      // test if op1 has enough remaining indices to satisfy
      const auto nidx_op1_remain =
          op1_cursor.range_iter()->size() -
          1;  // how many indices op1 has minus this index
      const auto nidx_op1_needs =
          (nop_connections[op1_idx] | target_nop_connections[op1_idx])
              .flip()
              .count();
      if (nidx_op1_needs > nidx_op1_remain) {
        if (!op1_op2_connected) {
          nop_connections[op1_idx].reset(op2_idx);
          nop_connections[op2_idx].reset(op1_idx);
        }
        return false;
      }

      // test if op2 has enough remaining indices to satisfy
      const auto nidx_op2_remain =
          op2_cursor.range_iter()->size() -
          1;  // how many indices op2 has minus this index
      const auto nidx_op2_needs =
          (nop_connections[op2_idx] | target_nop_connections[op2_idx])
              .flip()
              .count();
      if (nidx_op2_needs > nidx_op2_remain) {
        if (!op1_op2_connected) {
          nop_connections[op1_idx].reset(op2_idx);
          nop_connections[op2_idx].reset(op1_idx);
        }
        return false;
      }

      nop_adjacency_matrix[lowtri_idx(op1_idx, op2_idx)] += 1;
      update_topology(op1_idx);
      update_topology(op2_idx);

      return true;
    }
    /// @brief Updates connectivity when contraction is reversed
    template <typename Cursor>
    inline void disconnect(const container::svector<std::bitset<max_input_size>>
                               &target_nop_connections,
                           const Cursor &op1_cursor, const Cursor &op2_cursor) {
      auto update_topology = [this](size_t op_idx) {
        assert(nop_nconnections.at(op_idx) > 0);
        const auto nconnections = --nop_nconnections[op_idx];
        if (!topological_partitions.empty()) {
          auto partition_idx = nop_topological_partition[op_idx];
          if (nconnections == 0 && partition_idx > 0) {
            --partition_idx;  // to 0-based
            auto inserted =
                topological_partitions.at(partition_idx).insert(op_idx);
            assert(inserted.second);
          }
        }
      };

      // local vars
      const auto op1_idx = op1_cursor.range_ordinal();
      const auto op2_idx = op2_cursor.range_ordinal();
      update_topology(op1_idx);
      update_topology(op2_idx);
      if (target_nop_connections.empty())  // if no constraints, we don't keep
                                           // track of individual connections
        return;
      assert(nop_connections[op1_idx].test(op2_idx));

      auto &adjval = nop_adjacency_matrix[lowtri_idx(op1_idx, op2_idx)];
      assert(adjval > 0);
      adjval -= 1;
      if (adjval == 0) {
        nop_connections[op1_idx].reset(op2_idx);
        nop_connections[op2_idx].reset(op1_idx);
      }
    }
  };  // NontensorWickState

  /// Applies most naive version of Wick's theorem, where the sign rule involves
  /// counting Ops
  /// @return the result
  ExprPtr compute_nontensor_wick(const bool count_only) const {
    std::vector<std::pair<Product, std::shared_ptr<NormalOperator<S>>>>
        result;      //!< current value of the result
    std::mutex mtx;  // used in critical sections updating the result
    auto result_plus_mutex = std::make_pair(&result, &mtx);
    NontensorWickState state(input_, nop_topological_partition_);
    state.count_only = count_only;
    // TODO extract index->particle maps

    if (Logger::get_instance().wick_contract) {
      std::wcout << "nop topological partitions: {\n";
      for (auto &&toppart : state.topological_partitions) {
        std::wcout << "{" << std::endl;
        for (auto &&op_idx : toppart) {
          std::wcout << op_idx << std::endl;
        }
        std::wcout << "}" << std::endl;
      }
      std::wcout << "}" << std::endl;
    }

    recursive_nontensor_wick(result_plus_mutex, state);

    // if computing everything, include the contraction-free term
    if (!full_contractions_) {
      if (count_only) {
        ++state.count;
      } else {
        auto [phase, normop] = normalize(input_, state.input_partner_indices);
        result_plus_mutex.first->push_back(
            std::make_pair(Product(phase, {}), std::move(normop)));
      }
    }

    // convert result to an Expr
    // if result.size() == 0, return null ptr
    ExprPtr result_expr;
    if (count_only) {  // count only? return the total number as a Constant
      assert(result.empty());
      result_expr = ex<Constant>(state.count.load());
    } else if (result.size() == 1) {  // if result.size() == 1, return Product
      auto product = std::make_shared<Product>(std::move(result.at(0).first));
      if (full_contractions_)
        assert(result.at(0).second == nullptr);
      else {
        if (result.at(0).second)
          product->append(1, std::move(result.at(0).second));
      }
      result_expr = product;
    } else if (result.size() > 1) {
      auto sum = std::make_shared<Sum>();
      for (auto &&term : result) {
        if (full_contractions_) {
          assert(term.second == nullptr);
          sum->append(ex<Product>(std::move(term.first)));
        } else {
          auto term_product = std::make_shared<Product>(std::move(term.first));
          if (term.second) {
            term_product->append(1, term.second);
          }
          sum->append(term_product);
        }
      }
      result_expr = sum;
    } else if (result_expr == nullptr)
      result_expr = ex<Constant>(0);
    return result_expr;
  }

 public:
  virtual ~WickTheorem();

 private:
  void recursive_nontensor_wick(
      std::pair<
          std::vector<std::pair<Product, std::shared_ptr<NormalOperator<S>>>> *,
          std::mutex *> &result,
      NontensorWickState &state) const {
    using opseq_view_type = flattened_rangenest<NormalOperatorSequence<S>>;
    auto opseq_view = opseq_view_type(&state.opseq);
    using std::begin;
    using std::end;

    // if full contractions needed, make contractions involving first index with
    // another index, else contract any index i with index j (i<j)
    auto left_op_offset = state.left_op_offset;
    auto op_left_iter =
        ranges::next(begin(opseq_view), left_op_offset);  // left op to contract

    // optimization: can't contract fully if first op is not a qp annihilator
    if (full_contractions_ && !is_qpannihilator(*op_left_iter, input_.vacuum()))
      return;

    const auto op_left_iter_fence =
        full_contractions_ ? ranges::next(op_left_iter) : end(opseq_view);
    for (; op_left_iter != op_left_iter_fence;
         ++op_left_iter, ++left_op_offset) {
      auto op_right_iter = ranges::next(op_left_iter);
      for (; op_right_iter != end(opseq_view);) {
        if (op_right_iter != op_left_iter &&
            ranges::get_cursor(op_right_iter).range_iter() !=
                ranges::get_cursor(op_left_iter)
                    .range_iter()  // can't contract within same normop
        ) {
          // computes topological degeneracy:
          // 0 = nonunique index
          // n>0 = unique index in a group of n indices
          auto topological_degeneracy = [&]() {
            size_t result = 1;
            if (use_topology_) {
              auto &opseq_right =
                  *(ranges::get_cursor(op_right_iter).range_iter());
              auto &op_right_it = ranges::get_cursor(op_right_iter).elem_iter();
              auto op_right_idx_in_opseq =
                  op_right_it - ranges::begin(opseq_right);
              auto &hug_right = opseq_right.hug();
              auto &group_right = hug_right->group(op_right_idx_in_opseq);
              if (group_right.second.find(op_right_idx_in_opseq) ==
                  group_right.second.begin())
                result = hug_right->group_size(op_right_idx_in_opseq);
              else
                result = 0;
            }
            // account for topologically-equivalent normal operators
            if (result > 0 && !state.topological_partitions.empty()) {
              auto opseq_right_idx =
                  ranges::get_cursor(op_right_iter)
                      .range_ordinal();  // the index of normal operator
              auto opseq_right_toppart_idx = state.nop_topological_partition.at(
                  opseq_right_idx);  // the partition to which this normal
                                     // operator belongs to (0 = none)
              if (opseq_right_toppart_idx > 0) {  // if part of a partition ...
                --opseq_right_toppart_idx;        // to 0-based
                const auto &opseq_right_toppart =
                    state.topological_partitions.at(opseq_right_toppart_idx);
                if (!opseq_right_toppart
                         .empty()) {  // ... and the partition is not empty ...
                  const auto it = opseq_right_toppart.find(opseq_right_idx);
                  // .. and not missing from the partition (because then it's
                  // topologically unique) ...
                  if (it != opseq_right_toppart.end()) {
                    // ... and first in the partition
                    if (it == opseq_right_toppart.begin()) {
                      // account for the entire partition by scaling the
                      // contribution from the first contraction from this
                      // normal operator
                      result *= opseq_right_toppart.size();
                    } else
                      result = 0;
                  }
                }
              }
            }
            return result;
          };

          // check if can contract these indices and
          // check connectivity constraints (if needed)
          int64_t top_degen;
          if (can_contract(*op_left_iter, *op_right_iter, input_.vacuum())) {
            if ((top_degen = topological_degeneracy()) > 0) {
              if (state.connect(nop_connections_,
                                ranges::get_cursor(op_right_iter),
                                ranges::get_cursor(op_left_iter))) {
                if (Logger::get_instance().wick_contract) {
                  std::wcout << "level " << state.level << ":contracting "
                             << to_latex(*op_left_iter) << " with "
                             << to_latex(*op_right_iter)
                             << " (top_degen=" << top_degen << ")" << std::endl;
                  std::wcout << " current opseq = " << to_latex(state.opseq)
                             << std::endl;
                }

                // update the phase, if needed
                int64_t phase = 1;
                if (statistics == Statistics::FermiDirac) {
                  const auto distance =
                      ranges::get_cursor(op_right_iter).ordinal() -
                      ranges::get_cursor(op_left_iter).ordinal() - 1;
                  if (distance % 2) {
                    phase *= -1;
                  }
                }

                // update the prefactor and opseq
                Product sp_copy = state.sp;
                state.sp.append(
                    top_degen * phase,
                    contract(*op_left_iter, *op_right_iter, input_.vacuum()));

                // update the stats
                ++stats_.num_attempted_contractions;

                // remove from back to front
                Op<S> right = *op_right_iter;
                ranges::get_cursor(op_right_iter).erase();
                --state.opseq_size;
                Op<S> left = *op_left_iter;
                ranges::get_cursor(op_left_iter).erase();
                --state.opseq_size;

                //            std::wcout << "  opseq after contraction = " <<
                //            to_latex(state.opseq) << std::endl;

                // if have a nonzero result ...
                if (!state.sp.empty()) {
                  // update the result if nothing left to contract and
                  if (!full_contractions_ ||
                      (full_contractions_ && state.opseq_size == 0)) {
                    if (!state.count_only) {
                      if (full_contractions_) {
                        result.second->lock();
                        //              std::wcout << "got " <<
                        //              to_latex(state.sp)
                        //              << std::endl;
                        result.first->push_back(std::make_pair(
                            std::move(state.sp.deep_copy()),
                            std::shared_ptr<NormalOperator<S>>{}));
                        //              std::wcout << "now up to " <<
                        //              result.first->size()
                        //              << " terms" << std::endl;
                        result.second->unlock();
                      } else {
                        auto [phase, op] = normalize(
                            state.opseq, state.make_target_partner_indices());
                        result.second->lock();
                        result.first->push_back(std::make_pair(
                            std::move(state.sp.deep_copy().scale(phase)),
                            op->empty() ? nullptr : std::move(op)));
                        result.second->unlock();
                      }
                    } else
                      ++state.count;

                    // update the stats: count this contraction as useful
                    ++stats_.num_useful_contractions;
                  }
                }

                if (state.opseq_size != 0) {
                  const auto current_num_useful_contractions =
                      stats_.num_useful_contractions.load();
                  ++state.level;
                  state.left_op_offset = left_op_offset;
                  recursive_nontensor_wick(result, state);
                  --state.level;
                  // this contraction is useful if it leads to useful
                  // contractions as a result
                  if (current_num_useful_contractions !=
                      stats_.num_useful_contractions.load())
                    ++stats_.num_useful_contractions;
                }

                // restore the prefactor and opseq
                state.sp = std::move(sp_copy);
                // restore from front to back
                ranges::get_cursor(op_left_iter).insert(std::move(left));
                ++state.opseq_size;
                ranges::get_cursor(op_right_iter).insert(std::move(right));
                ++state.opseq_size;
                state.disconnect(nop_connections_,
                                 ranges::get_cursor(op_right_iter),
                                 ranges::get_cursor(op_left_iter));
                //            std::wcout << "  restored opseq = " <<
                //            to_latex(state.opseq) << std::endl;
              }  // connect succeeded
            }    // topologically-unique contraction
          }      // can_contract
          ++op_right_iter;
        } else {
          ++op_right_iter;
        }
      }  // right op iter
    }    // left op iter
  }

 public:
  static bool can_contract(const Op<S> &left, const Op<S> &right,
                           Vacuum vacuum = get_default_context().vacuum()) {
    // can only do Wick's theorem for physical vacuum (or similar)
    if constexpr (statistics == Statistics::BoseEinstein)
      assert(vacuum == Vacuum::Physical);

    if (is_qpannihilator<S>(left, vacuum) && is_qpcreator<S>(right, vacuum)) {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      if (qpspace_common != IndexSpace::null_instance()) return true;
    }
    return false;
  }

  static ExprPtr contract(const Op<S> &left, const Op<S> &right,
                          Vacuum vacuum = get_default_context().vacuum()) {
    assert(can_contract(left, right, vacuum));
    //    assert(
    //        !left.index().has_proto_indices() &&
    //            !right.index().has_proto_indices());  // I don't think the
    //            logic is
    // correct for dependent indices
    if (is_pure_qpannihilator<S>(left, vacuum) &&
        is_pure_qpcreator<S>(right, vacuum))
      return make_overlap(left.index(), right.index());
    else {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      const auto index_common = Index::make_tmp_index(qpspace_common);

      // preserve bra/ket positions of left & right
      const auto left_is_ann = left.action() == Action::annihilate;
      assert(left_is_ann || right.action() == Action::annihilate);

      if (qpspace_common != left.index().space() &&
          qpspace_common !=
              right.index().space()) {  // may need 2 overlaps if neither space
        // is pure qp creator/annihilator
        auto result = std::make_shared<Product>();
        result->append(1, left_is_ann
                              ? make_overlap(left.index(), index_common)
                              : make_overlap(index_common, left.index()));
        result->append(1, left_is_ann
                              ? make_overlap(index_common, right.index())
                              : make_overlap(right.index(), index_common));
        return result;
      } else {
        return left_is_ann ? make_overlap(left.index(), right.index())
                           : make_overlap(right.index(), left.index());
      }
    }
  }

  /// @param[in,out] on input, Wick's theorem result, on output the result of
  /// reducing the overlaps
  void reduce(ExprPtr &expr) const;
};

using BWickTheorem = WickTheorem<Statistics::BoseEinstein>;
using FWickTheorem = WickTheorem<Statistics::FermiDirac>;

}  // namespace sequant

#include "wick.impl.hpp"

#endif  // SEQUANT_WICK_HPP
