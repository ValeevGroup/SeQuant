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
  static_assert(S == Statistics::FermiDirac,
                "WickTheorem not yet implemented for Bose-Einstein");

  WickTheorem(const WickTheorem &) = default;
  WickTheorem(WickTheorem &&) = default;
  WickTheorem &operator=(const WickTheorem &) = default;
  WickTheorem &operator=(WickTheorem &&) = default;

  explicit WickTheorem(const NormalOperatorSequence<S> &input) : input_(input) {
    assert(input.size() <= max_input_size);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
  }

  explicit WickTheorem(ExprPtr expr_input) : expr_input_(expr_input) {}

  /// constructs WickTheorem from @c other with expression input set to @c
  /// expr_input
  WickTheorem(ExprPtr expr_input, const WickTheorem &other)
      : WickTheorem(other) {
    // copy ctor does not do anything useful, so this is OK
    expr_input_ = expr_input;
  }

  /// Controls whether next call to compute() will full contractions only or all
  /// (including partial) contractions. By default compute() generates all
  /// contractions.
  /// @param sf if true, will complete full contractions only.
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
  WickTheorem &set_op_connections(IndexPairContainer &&op_index_pairs) {
    if (expr_input_ == nullptr || !op_connections_input_.empty()) {
      for (const auto &opidx_pair : op_index_pairs) {
        if (opidx_pair.first < 0 || opidx_pair.first >= input_.size()) {
          throw std::invalid_argument(
              "WickTheorem::set_op_connections: op index out of range");
        }
        if (opidx_pair.second < 0 || opidx_pair.second >= input_.size()) {
          throw std::invalid_argument(
              "WickTheorem::set_op_connections: op index out of range");
        }
      }
      if (op_index_pairs.size() != 0) {
        op_connections_.resize(input_.size());
        for (auto &v : op_connections_) {
          v.set();
        }
        for (const auto &opidx_pair : op_index_pairs) {
          op_connections_[opidx_pair.first].reset(opidx_pair.second);
          op_connections_[opidx_pair.second].reset(opidx_pair.first);
        }
      }
      op_connections_input_.clear();
    } else {
      ranges::for_each(op_index_pairs, [this](const auto &idxpair) {
        op_connections_input_.push_back(idxpair);
      });
    }

    return *this;
  }

  /// @tparam Integer an integral type
  template <typename Integer = long>
  WickTheorem& set_op_connections(std::initializer_list<std::pair<Integer,Integer>> op_index_pairs) {
    return this->set_op_connections<const decltype(op_index_pairs)&>(op_index_pairs);
  }
  ///@}

  /// @name topological partition specifiers
  ///
  /// Specifies topological partition of normal operators; free (non-connected) operators
  /// in the same partition are considered topologically equivalent, hence if
  /// only full contractions are needed only contractions to the first available
  /// operator in a partition is needed (multiplied by the degeneracy)
  /// @param op_partitions list of operator partitions
  /// @note if this partitions are not given, every operator is assumed to be in
  /// its own partition
  /// @internal this performs only the first phase of initialization of
  /// op_topological_partition_
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
    auto current_nops = size(op_topological_partition_);
    for (auto &&partition : op_partitions) {
      for (auto &&op_idx : partition) {
        assert(op_idx >= 0);
        if (op_idx >= current_nops) {
          current_nops = upsize_op_topological_partition(op_idx + 1);
        }
        op_topological_partition_[op_idx] = partition_cnt + 1;
      }
      ++partition_cnt;
    }
    return *this;
  }

  /// @tparam Integer an integral type
  template <typename Integer = long>
  WickTheorem& set_op_partitions(std::initializer_list<std::initializer_list<Integer>> op_partitions) {
    return this->set_op_partitions<const decltype(op_partitions)&>(op_partitions);
  }
  ///@}

  /// Computes and returns the result
  /// @param count_only if true, will return a vector of default-initialized
  /// values, useful if only interested in the total count
  /// @return the result of applying Wick's theorem
  ExprPtr compute(const bool count_only = false);

 private:
  static constexpr size_t max_input_size =
      32;  // max # of operators in the input sequence

  // if nonnull, apply wick to the whole expression recursively, else input_ is
  // set this is mutated by compute
  mutable ExprPtr expr_input_;

  mutable NormalOperatorSequence<S> input_;
  bool full_contractions_ = false;
  bool spinfree_ = false;
  bool use_topology_ = false;
  container::set<Index> external_indices_;
  // for each operator specifies the reverse bitmask of target connections
  // (0 = must connect)
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  container::svector<std::bitset<max_input_size>> op_connections_;
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  container::vector<std::pair<size_t, size_t>>
      op_connections_input_;  // only used to cache input to set_op_connections_

  // for each operator specifies its topological partition (0 = topologically
  // unique)
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  mutable container::svector<size_t> op_topological_partition_;

  /// upsizes op_topological_partition_, filling new entries with zeroes
  /// noop if current size > new_size
  /// @return the (updated) size of op_topological_partition_
  /// TODO rename op -> nop to distinguish Op and NormalOperator
  size_t upsize_op_topological_partition(size_t new_size) const {
    using std::size;
    const auto current_size = size(op_topological_partition_);
    if (new_size > current_size) {
      op_topological_partition_.resize(new_size);
      for (size_t i = current_size; i != new_size; ++i)
        op_topological_partition_[i] = 0;
      return new_size;
    } else
      return current_size;
  }

  /// Evaluates wick_ theorem for a single NormalOperatorSequence
  /// @return the result of applying Wick's theorem
  ExprPtr compute_nopseq(const bool count_only) const {
    if (!full_contractions_)
      throw std::logic_error(
          "WickTheorem::compute: full_contractions=false not yet supported");
    if (spinfree_)
      throw std::logic_error(
          "WickTheorem::compute: spinfree=true not yet supported");
    // process cached op_connections_input_, if needed
    if (!op_connections_input_.empty())
      const_cast<WickTheorem<S> &>(*this).set_op_connections(
          op_connections_input_);
    // size op_topological_partition_ to match input_, if needed
    upsize_op_topological_partition(input_.size());
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
          count_only(false),
          op_connections(opseq.size()),
          adjacency_matrix(opseq.size() * (opseq.size() - 1) / 2, 0),
          op_nconnections(opseq.size(), 0),
          op_topological_partition(op_toppart) {
      init_topological_partitions();
    }

    NontensorWickState(const NontensorWickState&) = delete;
    NontensorWickState(NontensorWickState&&) = delete;
    NontensorWickState& operator=(const NontensorWickState&) = delete;
    NontensorWickState& operator=(NontensorWickState&&) = delete;

    NormalOperatorSequence<S> opseq;  //!< current state of operator sequence
    std::size_t opseq_size;           //!< current size of opseq
    Product sp;                       //!< current prefactor
    int level;                        //!< level in recursive wick call stack
    bool count_only;                  //!< if true, only update result size
    /// TODO rename op -> nop to distinguish Op and NormalOperator
    container::svector<std::bitset<max_input_size>>
        op_connections;  //!< bitmask of connections for each op (1 = connected)
    container::svector<size_t>
        adjacency_matrix;  //!< number of connections between each normop, only
                           //!< lower triangle is kept

    /// for each operator specifies how many connections it currently has
    /// TODO rename op -> nop to distinguish Op and NormalOperator
    container::svector<size_t> op_nconnections;
    /// maps op to its topological partition index (1-based, 0 = no partition)
    /// TODO rename op -> nop to distinguish Op and NormalOperator
    container::svector<size_t> op_topological_partition;
    /// current state of partitions (will only match op_topological_partition before any contractions have occurred)
    /// - when an operator is connected it's removed from the partition
    /// - when it is disconnected fully it's re-added to the partition
    container::svector<container::set<size_t>> topological_partitions;

    // populates partitions using the data from op_topological_partition
    void init_topological_partitions() {
      // partition indices in op_topological_partition are 1-based
      const auto npartitions = *ranges::max_element(op_topological_partition);
      topological_partitions.resize(npartitions);
      size_t op_cnt = 0;
      ranges::for_each(op_topological_partition, [this,&op_cnt](size_t toppart_idx) {
        if (toppart_idx > 0) {  // in a partition
          topological_partitions.at(toppart_idx - 1).insert(op_cnt);
        }
        ++op_cnt;
      });
      // assert that we don't have empty partitions due to invalid contents of op_topological_partition
      assert(ranges::any_of(topological_partitions, [](auto&& partition){
        return partition.empty();
      }) == false);
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
                            &target_op_connections,
                        const Cursor &op1_cursor, const Cursor &op2_cursor) {
      auto result = true;

      auto update_topology = [this](size_t op_idx) {
        const auto nconnections = op_nconnections[op_idx];
        // if using topological partitions for normal ops, and this operator is in one of them, remove it on first connection
        if (!topological_partitions.empty()) {
          auto partition_idx = op_topological_partition[op_idx];
          if (nconnections == 0 && partition_idx > 0) {
            --partition_idx;  // to 0-based
            assert(topological_partitions.at(partition_idx).size() > 0);
            auto removed = topological_partitions[partition_idx].erase(op_idx);
            assert(removed);
          }
        }
        ++op_nconnections[op_idx];
      };

      // local vars
      const auto op1_idx = op1_cursor.range_ordinal();
      const auto op2_idx = op2_cursor.range_ordinal();
      if (target_op_connections
              .empty()) {  // if no constraints, all is fair game
        update_topology(op1_idx);
        update_topology(op2_idx);
        return true;
      }
      const auto op1_op2_connected = op_connections[op1_idx].test(op2_idx);

      // update the connectivity
      if (!op1_op2_connected) {
        op_connections[op1_idx].set(op2_idx);
        op_connections[op2_idx].set(op1_idx);
      }

      // test if op1 has enough remaining indices to satisfy
      const auto nidx_op1_remain =
          op1_cursor.range_iter()->size() -
          1;  // how many indices op1 has minus this index
      const auto nidx_op1_needs =
          (op_connections[op1_idx] | target_op_connections[op1_idx])
              .flip()
              .count();
      if (nidx_op1_needs > nidx_op1_remain) {
        if (!op1_op2_connected) {
          op_connections[op1_idx].reset(op2_idx);
          op_connections[op2_idx].reset(op1_idx);
        }
        return false;
      }

      // test if op2 has enough remaining indices to satisfy
      const auto nidx_op2_remain =
          op2_cursor.range_iter()->size() -
          1;  // how many indices op2 has minus this index
      const auto nidx_op2_needs =
          (op_connections[op2_idx] | target_op_connections[op2_idx])
              .flip()
              .count();
      if (nidx_op2_needs > nidx_op2_remain) {
        if (!op1_op2_connected) {
          op_connections[op1_idx].reset(op2_idx);
          op_connections[op2_idx].reset(op1_idx);
        }
        return false;
      }

      adjacency_matrix[lowtri_idx(op1_idx, op2_idx)] += 1;
      update_topology(op1_idx);
      update_topology(op2_idx);

      return true;
    }
    /// @brief Updates connectivity when contraction is reversed
    template <typename Cursor>
    inline void disconnect(const container::svector<std::bitset<max_input_size>>
                               &target_op_connections,
                           const Cursor &op1_cursor, const Cursor &op2_cursor) {
      auto update_topology = [this](size_t op_idx) {
        assert(op_nconnections.at(op_idx) > 0);
        const auto nconnections = --op_nconnections[op_idx];
        if (!topological_partitions.empty()) {
          auto partition_idx = op_topological_partition[op_idx];
          if (nconnections == 0 && partition_idx > 0) {
            --partition_idx;  // to 0-based
            auto inserted = topological_partitions.at(partition_idx).insert(op_idx);
            assert(inserted.second);
          }
        }
      };

      // local vars
      const auto op1_idx = op1_cursor.range_ordinal();
      const auto op2_idx = op2_cursor.range_ordinal();
      update_topology(op1_idx);
      update_topology(op2_idx);
      if (target_op_connections.empty())  // if no constraints, we don't keep
                                          // track of individual connections
        return;
      assert(op_connections[op1_idx].test(op2_idx));

      auto &adjval = adjacency_matrix[lowtri_idx(op1_idx, op2_idx)];
      assert(adjval > 0);
      adjval -= 1;
      if (adjval == 0) {
        op_connections[op1_idx].reset(op2_idx);
        op_connections[op2_idx].reset(op1_idx);
      }
    }
  };  // NontensorWickState

  /// Applies most naive version of Wick's theorem, where sign rule involves
  /// counting Ops
  /// @return the result, or nullptr if the result is zero
  ExprPtr compute_nontensor_wick(const bool count_only) const {
    std::vector<std::pair<Product, NormalOperator<S>>>
        result;      //!< current value of the result
    std::mutex mtx;  // used in critical sections updating the result
    auto result_plus_mutex = std::make_pair(&result, &mtx);
    NontensorWickState state(input_, op_topological_partition_);
    state.count_only = count_only;

    recursive_nontensor_wick(result_plus_mutex, state);

    // convert result to an Expr
    // if result.size() == 0, return null ptr
    ExprPtr result_expr;
    if (result.size() == 1) {  // if result.size() == 1, return Product
      result_expr = ex<Product>(std::move(result[0].first));
    } else if (result.size() > 1) {
      auto sum = std::make_shared<Sum>();
      for (auto &term : result) {
        sum->append(ex<Product>(std::move(term.first)));
      }
      result_expr = sum;
    }
    return result_expr;
  }

  void recursive_nontensor_wick(
      std::pair<std::vector<std::pair<Product, NormalOperator<S>>> *,
                std::mutex *> &result,
      NontensorWickState &state) const {
    // if full contractions needed, make contractions involving first index with
    // another index, else contract any index i with index j (i<j)
    if (full_contractions_) {
      using opseq_view_type = flattened_rangenest<NormalOperatorSequence<S>>;
      auto opseq_view = opseq_view_type(&state.opseq);
      using std::begin;
      using std::end;
      auto opseq_view_begin = begin(opseq_view);

      // optimization: can't contract fully if first op is not a qp annihilator
      if (!is_qpannihilator(*opseq_view_begin, input_.vacuum())) return;

      auto op_iter = opseq_view_begin;
      ++op_iter;
      for (; op_iter != end(opseq_view);) {
        if (op_iter != opseq_view_begin &&
            ranges::get_cursor(op_iter).range_iter() !=
                ranges::get_cursor(opseq_view_begin)
                    .range_iter()  // can't contract within same normop
        ) {
          // computes topological degeneracy:
          // 0 = nonunique index
          // n>0 = unique index in a group of n indices
          auto topological_degeneracy = [&]() {
            size_t result = 1;
            if (use_topology_) {
              auto &opseq_right = *(ranges::get_cursor(op_iter).range_iter());
              auto &op_it = ranges::get_cursor(op_iter).elem_iter();
              auto op_idx_in_opseq = op_it - ranges::begin(opseq_right);
              auto &hug = opseq_right.hug();
              auto &group = hug->group(op_idx_in_opseq);
              if (group.second.find(op_idx_in_opseq) == group.second.begin())
                result = hug->group_size(op_idx_in_opseq);
              else
                result = 0;
            }
            // account for topologically-equivalent normal operators
            if (result > 0 && !state.topological_partitions.empty()) {
              auto opseq_right_idx =
                  ranges::get_cursor(op_iter)
                      .range_ordinal();  // the index of normal operator
              auto opseq_right_toppart_idx = state.op_topological_partition.at(
                  opseq_right_idx);  // the partition to which this normal
                                     // operator belongs to (0 = none)
              if (opseq_right_toppart_idx > 0) {  // if part of a partition ...
                --opseq_right_toppart_idx;        // to 0-based
                const auto &opseq_right_toppart =
                    state.topological_partitions.at(opseq_right_toppart_idx);
                if (!opseq_right_toppart
                         .empty()) {  // ... and the partition is not empty ...
                  const auto it = opseq_right_toppart.find(opseq_right_idx);
                  // .. and not missing from the partition (because then it's topologically unique) ...
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
          size_t top_degen;
          if (can_contract(*opseq_view_begin, *op_iter, input_.vacuum()) &&
              (top_degen = topological_degeneracy()) > 0 &&
              state.connect(op_connections_, ranges::get_cursor(op_iter),
                            ranges::get_cursor(opseq_view_begin))) {
            if (Logger::get_instance().wick_contract) {
              std::wcout << "level " << state.level << ":contracting "
                         << to_latex(*opseq_view_begin) << " with "
                         << to_latex(*op_iter) << " (top_degen=" << top_degen
                         << ")" << std::endl;
              std::wcout << " current opseq = " << to_latex(state.opseq)
                         << std::endl;
            }

            // update the phase, if needed
            double phase = 1;
            if (statistics == Statistics::FermiDirac) {
              const auto distance =
                  ranges::get_cursor(op_iter).ordinal() -
                  ranges::get_cursor(opseq_view_begin).ordinal() - 1;
              if (distance % 2) {
                phase *= -1;
              }
            }

            // update the prefactor and opseq
            Product sp_copy = state.sp;
            state.sp.append(
                top_degen * phase,
                contract(*opseq_view_begin, *op_iter, input_.vacuum()));
            // remove from back to front
            Op<S> right = *op_iter;
            ranges::get_cursor(op_iter).erase();
            --state.opseq_size;
            Op<S> left = *opseq_view_begin;
            ranges::get_cursor(opseq_view_begin).erase();
            --state.opseq_size;

            //            std::wcout << "  opseq after contraction = " <<
            //            to_latex(state.opseq) << std::endl;

            // update the result if nothing left to contract and have a
            // nonzero result
            if (state.opseq_size == 0 && !state.sp.empty()) {
              result.second->lock();
              //              std::wcout << "got " << to_latex(state.sp) <<
              //              std::endl;
              if (!state.count_only)
                result.first->push_back(std::make_pair(
                    std::move(state.sp.deep_copy()), NormalOperator<S>{}));
              else
                result.first->resize(result.first->size() + 1);
              //              std::wcout << "now up to " <<
              //              result.first->size()
              //              << " terms" << std::endl;
              result.second->unlock();
            }

            if (state.opseq_size != 0) {
              ++state.level;
              recursive_nontensor_wick(result, state);
              --state.level;
            }

            // restore the prefactor and opseq
            state.sp = std::move(sp_copy);
            // restore from front to back
            ranges::get_cursor(opseq_view_begin).insert(std::move(left));
            ++state.opseq_size;
            ranges::get_cursor(op_iter).insert(std::move(right));
            ++state.opseq_size;
            state.disconnect(op_connections_, ranges::get_cursor(op_iter),
                             ranges::get_cursor(opseq_view_begin));
            //            std::wcout << "  restored opseq = " <<
            //            to_latex(state.opseq) << std::endl;
          }
          ++op_iter;
        } else
          ++op_iter;
      }
    } else
      assert(false);  // full_contraction_=false not implemented yet, should
    // result in error earlier
  }

 public:
  static bool can_contract(const Op<S> &left, const Op<S> &right,
                           Vacuum vacuum = get_default_context().vacuum()) {
    if (is_qpannihilator<S>(left, vacuum) && is_qpcreator<S>(right, vacuum)) {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      if (qpspace_common != IndexSpace::null_instance()) return true;
    }
    return false;
  }

  static std::shared_ptr<Expr> contract(
      const Op<S> &left, const Op<S> &right,
      Vacuum vacuum = get_default_context().vacuum()) {
    assert(can_contract(left, right, vacuum));
    //    assert(
    //        !left.index().has_proto_indices() &&
    //            !right.index().has_proto_indices());  // I don't think the
    //            logic is
    // correct for dependent indices
    if (is_pure_qpannihilator<S>(left, vacuum) &&
        is_pure_qpcreator<S>(right, vacuum))
      return overlap(left.index(), right.index());
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
        result->append(1, left_is_ann ? overlap(left.index(), index_common)
                                      : overlap(index_common, left.index()));
        result->append(1, left_is_ann ? overlap(index_common, right.index())
                                      : overlap(right.index(), index_common));
        return result;
      } else {
        return left_is_ann ? overlap(left.index(), right.index())
                           : overlap(right.index(), left.index());
      }
    }
  }

  /// @param[in,out] on input, Wick theorem result, on output the result of
  /// reducing the overlaps
  void reduce(ExprPtr &expr) const;
};

using BWickTheorem = WickTheorem<Statistics::BoseEinstein>;
using FWickTheorem = WickTheorem<Statistics::FermiDirac>;

}  // namespace sequant

#include "wick.impl.hpp"

#endif  // SEQUANT_WICK_HPP
