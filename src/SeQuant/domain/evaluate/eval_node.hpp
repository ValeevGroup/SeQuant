#ifndef SEQUANT_EVALUATE_EVAL_NODE_HPP
#define SEQUANT_EVALUATE_EVAL_NODE_HPP

#include "SeQuant/core/container.hpp"
#include "eval_fwd.hpp"

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <numeric>

namespace sequant::evaluate {

///
/// @brief Representation of binary evaluation of sequant::Expr objects.
///
/// The atomic SeQuant Expr(s) can be thought of either as a representation of a
/// single data-tensor (eg. btas::Tensor<double>, TA::TArrayD), or binary
/// evaluations of such data-tensors. An evaluation could be summing two
/// data-tensors or taking a product of them. A result of a binary evaluation
/// can be an input for another evaluation, i.e. evaluations can be nested.
///
/// The classes herein represent such evaluations.
///
/// @author Bimal Gaudel
/// @date Apr 12, 2020
///

// forward declare
template <typename DataTensorType>
class EvalNodeInternal;

// forward declare
template <typename DataTensorType>
class EvalNodeLeaf;

///
/// Base class for evaluation node.
/// \tparam DataTensorType Type of the data tensor. eg. TA::TArrayD from
///                        TiledArray
///
template <typename DataTensorType>
class EvalNode {
 protected:
  /// The index labels of the tensor's bra and ket in that order.
  IndexContainer indices_;

  /// A unique identifier of this evaluation.
  HashType hash_value_{0};

  /// The scalar to multiply this evaluation with.
  ScalarType scalar_{1};

  /// Hash the current node and store the value to the hash_value_.
  virtual void hash_node() = 0;

 public:
  /// Construct evaluation tree out of sequant expression.
  EvalNode(const ExprPtr& expr) {
    if (expr->is<Tensor>()) {
      EvalNodeLeaf<DataTensorType>{expr};
    } else if (expr->is<Product>() || expr->is<Sum>()) {
      EvalNodeInternal<DataTensorType>{expr};
    } else {
      throw std::domain_error("Only sum, product or tensor is allowed!");
    }
  }

  /// Make a eval tree pointer.
  EvalTreePtr<DataTensorType> make_tree_ptr(const ExprPtr& expr) {
    return std::make_shared<EvalNode<DataTensorType>>(
        EvalNode<DataTensorType>(expr));
  }

  /// Getter method of the indices_ field.
  const IndexContainer& indices() const { return indices_; }

  /// Getter method of the hash_value_ field.
  HashType hash_value() const { return hash_value_; }

  /// Getter method of the scalar_ field.
  ScalarType scalar() const { return scalar_; }

  /// Set the scalar of the node.
  void scale(ScalarType fac) { scalar_ = fac; }

  /// Check if this is one of the end tensors.
  virtual bool is_leaf() const {
    throw std::logic_error("is_leaf() not implemented in this derived class.");
  }

  /// Get operations count for this evaluation.
  /// \param ispace_size_map A map from IndexSpace type to the size of the
  ///                       space.
  ///                       e.g. IndexSpace::active_occupied -> 5
  ///                       e.g. IndexSpace::active_unoccupied -> 20
  virtual OpsCount ops_count(const container::map<IndexSpace::TypeAttr, size_t>&
                                 ispace_size_map) const {
    throw std::logic_error(
        "ops_count() not implemented in this derived class.");
  }

  /// Visit the tree by pre-order traversal.
  virtual void visit(const std::function<void(const EvalNode<DataTensorType>&)>&
                         visitor) const {
    throw std::logic_error("visit() not implemented in this derived class.");
  }

  /// Evaluate the EvalNode.
  /// \param context A map from the hash values in the evaluation tree to the
  /// data tensors stored externally outside this class.
  virtual DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const {
    throw std::logic_error("evaluate() not implemented in this derived class.");
  }

  /// Get a directed graph definitions and paths in dot language.
  /// virtual std::wstring to_digraph() const = 0;
};

/// Evaluation node internal.
template <typename DataTensorType>
class EvalNodeInternal : public EvalNode<DataTensorType> {
 private:
  /// Evaluation node to the left.
  EvalTreePtr<DataTensorType> left_{nullptr};

  /// Evaluation node to the right.
  EvalTreePtr<DataTensorType> right_{nullptr};

  /// Type of operation for this internal node.
  Operation operation_;

  /// Hash the internal node and set hash_value_.
  void hash_node() override {
    auto imed_operation = this->operation();

    boost::hash<size_t> number_hasher;

    HashType imed_hash_value;
    // hashing the type of binary operation by casting operation into size_t
    imed_hash_value = number_hasher(static_cast<size_t>(imed_operation));

    // operation is of (anti)symmetrization type then its hash value is solely
    // determined by which indices are being permuted, this information is
    // carried by the left tensor -- the tensor with label 'A'(or 'P')
    if ((imed_operation == Operation::ANTISYMMETRIZE) ||
        (imed_operation == Operation::SYMMETRIZE)) {
      boost::hash_combine(imed_hash_value, left()->hash_value());
      this->hash_value_ = imed_hash_value;
      return;
    }

    // hash the hash values of the left and right tensors
    // but do so in an order-agnostic manner
    HashType lhash_value = this->left()->hash_value();
    HashType rhash_value = this->right()->hash_value();

    if (lhash_value < rhash_value) {
      boost::hash_combine(imed_hash_value, lhash_value);
      boost::hash_combine(imed_hash_value, rhash_value);
    } else {
      boost::hash_combine(imed_hash_value, rhash_value);
      boost::hash_combine(imed_hash_value, lhash_value);
    }

    // when Operation is product ie. a tensor contraction, not only the left and
    // the right tensor's hash values are important but also the information
    // about which of their indices got contracted is needed to uniquely
    // identify such an evaluation
    if (imed_operation == Operation::PRODUCT) {
      // a lambda expression that combines the hashes of the non-contracting
      // indices of right or left tensor
      auto index_hash_combiner =
          [&](const EvalTreePtr<DataTensorType>& lr_node) {
            // iterate through the imed_ptr indices
            for (const auto& idx : this->indices()) {
              // iterate through the indices of the passed left/right tensor
              // by a counter
              for (auto i = 0; i < lr_node->indices().size(); ++i) {
                // if the index of the imed_ptr matches with the lr_tensor
                // then time to break as we found the contracted indices
                // no need to look further
                if (idx.label() == lr_node->indices()[i].label()) {
                  break;
                }
                // combine the hash of the counter of non-contracted index
                boost::hash_combine(imed_hash_value, number_hasher(i));
              }
            }
          };  // lambda function index_hash_combiner

      // when the same contraction occurs with left tensor switched to the right
      // we should obtain the same hash value. ie hashing should be agnostic to
      // the order of the two contracting tensors
      if (lhash_value < rhash_value) {
        index_hash_combiner(this->left());
        index_hash_combiner(this->right());
      } else {
        index_hash_combiner(this->right());
        index_hash_combiner(this->left());
      }
    }  // done combining hashes for product type evaluation
    this->hash_value_ = imed_hash_value;
    return;
  }

 public:
  /// Get left node.
  const EvalTreePtr<DataTensorType>& left() const { return left_; }

  /// Get right node.
  const EvalTreePtr<DataTensorType>& right() const { return right_; }

  /// Get operation type of this internal node.
  Operation operation() const { return operation_; }

  /// Internal node is not leaf.
  bool is_leaf() const override { return false; }

  /// Get operations count for internal node.
  OpsCount ops_count(const container::map<IndexSpace::TypeAttr, size_t>&
                         ispace_size_map) const override {
    // count operations on left and right node
    auto count = this->left()->ops_count(ispace_size_map) +
                 this->right()->ops_count(ispace_size_map);

    if (this->operation() == Operation::PRODUCT) {
      // the ops count for a product is the product of the sizes of the index
      // space of the unique indices in the left and right node combined
      auto& left_indices = this->left()->indices();
      auto& right_indices = this->right()->indices();

      container::set<Index> unique_indices;
      for (const auto& idx : left_indices) {
        unique_indices.insert(idx);
      }
      for (const auto& idx : right_indices) {
        unique_indices.insert(idx);
      }

      auto contractions_ops = 1;
      for (const auto& idx : unique_indices) {
        contractions_ops *= (ispace_size_map.find(idx.space().type()))->second;
      }
      count += contractions_ops;
    }  // operation type PRODUCT

    return count;
  }

  /// Visit internal node.
  void visit(const std::function<void(const EvalNode<DataTensorType>&)>&
                 visitor) const override {
    // visit node
    visitor(*this);
    // visit left
    this->left()->visit(visitor);
    // visit right
    this->right()->visit(visitor);
  }

  /// Evaluate internal node.
  DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const override {
    auto operation = this->operation();

    if (operation == Operation::ANTISYMMETRIZE) {
      auto right = this->right()->evaluate(context);
      auto bra_rank = this->indices().size() / 2;
      auto ket_rank = this->indices().size() - bra_rank;
      return antisymmetrize(right, bra_rank, ket_rank);
    }

    // generates tiledarray annotation based on a node's index labels
    // @note this wouldn't be necessary if the tensor algebra library
    // would support std::string_view as annotations
    auto TA_annotation = [this](decltype(this->indices())& indices) {
      std::string annot = "";
      for (const auto& idx : indices)
        annot += std::string(idx.label().begin(), idx.label().end()) + ", ";

      annot.erase(annot.size() - 2);  // remove trailing L", "
      return annot;
    };

    auto left_annot = TA_annotation(this->left()->indices());
    auto right_annot = TA_annotation(this->right()->indices());
    auto this_annot = TA_annotation(this->indices());

    if (operation == Operation::SUM) {
      // sum left and right evaluated tensors
      // using tiled array syntax
      DataTensorType result;
      result(this_annot) =
          left()->scalar() * left()->evaluate(context)(left_annot) +
          right()->scalar() * right()->evaluate(context)(right_annot);
      return result;
    } else if (operation == Operation::PRODUCT) {
      // contract left and right evaluated tensors
      // using tiled array syntax
      DataTensorType result;
      result(this_annot) =
          left()->scalar() * left()->evaluate(context)(left_annot) *
          right()->scalar() * right()->evaluate(context)(right_annot);
      return result;
    } else {
      throw std::domain_error("Operation: " + std::to_string(operation) +
                              " not supported!");
    }
  }

  /// Construct internal node.
  EvalNodeInternal(const ExprPtr& expr) {
    if (expr->is<Sum>()) {
      EvalNodeInternal(expr->as<Sum>());
    } else if (expr->is<Product>()) {
      EvalNodeInternal(expr->as<Product>());
    } else {
      throw std::logic_error(
          "Only Sum or Product is valid for internal node constructor.");
    }
  }

 private:
  /// Private constructor for internal node.
  /// \param sum sequant Sum.
  EvalNodeInternal(const Sum& sum) {
    auto sum_accumulator = [this](const EvalTreePtr<DataTensorType>& lexpr,
                                  const ExprPtr& summand) {
      return EvalNodeInternal<DataTensorType>(
          lexpr, this->make_tree_ptr(summand), Operation::SUM);
    };

    *this =
        std::accumulate(sum.begin() + 1, sum.end(),
                        this->make_tree_ptr(sum.summand(0)), sum_accumulator);
  }

  /// Private constructor for internal node.
  /// \param prod sequant Product.
  EvalNodeInternal(const Product& prod) {
    auto prod_accumulator = [this](const EvalTreePtr<DataTensorType>& lexpr,
                                   const ExprPtr& factor) {
      return EvalNodeInternal<DataTensorType>(
          lexpr, this->make_tree_ptr(factor), Operation::PRODUCT);
    };

    auto& fac0 = prod.factor(0);
    if (auto label = fac0->is<Tensor>() ? fac0->as<Tensor>().label() : L"";
        label == L"A" || label == L"P") {
      // (anti-)symmetrization tensor encountered
      auto right = std::make_shared<Product>(prod.begin() + 1, prod.end());
      right->scale(prod.scalar());
      //
      EvalNodeInternal(
          this->make_tree_ptr(fac0), this->make_tree_ptr(right),
          label == L"A" ? Operation::ANTISYMMETRIZE : Operation::SYMMETRIZE);
    } else {
      auto init = this->make_tree_ptr(fac0);
      init->scale(prod.scalar().real());
      *this =
          std::accumulate(prod.begin() + 1, prod.end(), init, prod_accumulator);
    }
  }

  /// Private constructor for internal node.
  EvalNodeInternal(const EvalTreePtr<DataTensorType>& left,
                   const EvalTreePtr<DataTensorType>& right,
                   Operation operation)
      : operation_{operation}, left_{left}, right_{right} {
    if (operation == Operation::SUM) {
      if (left()->indices().size() != right()->indices().size())
        throw std::domain_error(
            "While summing two tensors, their number of indices must match!");
      // indices for the intermediate
      this->indices_ =
          !left()->is_leaf()
              ? left()->indices()
              : !right()->is_leaf() ? right()->indices() : IndexContainer{};
      // couldn't decide on the indices of the result because
      // both nodes were leaves and leaf indices may not be sorted
      if (this->indices().empty()) {
        this->indices_ = left()->indices();  // or right()->indices()
        std::sort(this->indices_.begin(), this->indices_.end());
      }
    } else if (operation == Operation::PRODUCT) {
      auto lindices = left()->indices();
      auto rindices = right()->indices();
      // indices from leaves may not be sorted
      if (!left()->is_leaf()) std::sort(lindices.begin(), lindices().end());
      if (!right()->is_leaf()) std::sort(rindices.begin(), rindices().end());
      // set the non-contracting indices
      std::set_symmetric_difference(lindices.begin(), lindices.end(),
                                    rindices.begin(), rindices.end(),
                                    std::back_inserter(this->indices_));
    } else if ((operation == Operation::ANTISYMMETRIZE) ||
               (operation == Operation::SYMMETRIZE)) {
      this->indices_ = left()->indices();
      // std::sort(this->indices_.begin(), this->indices_.end());
    } else {
      throw std::logic_error("Invalid Operation while forming intermediate");
    }
    // set the hash_value_
    hash_node();
  }

  /// Antisymmetrize DataTensorType
  /// \param ta_tensor TiledArray tensor.
  /// \param bra_rank rank of the tensor bra.
  /// \param ket_rank rank of the tensor ket.
  DataTensorType antisymmetrize(const DataTensorType& ta_tensor,
                                size_t bra_rank, size_t ket_rank) {
    container::svector<size_t> bra_indices(bra_rank), ket_indices(ket_rank);
    // {0, 1, .. bra_rank - 1}
    std::iota(bra_indices.begin(), bra_indices.end(), 0);
    // {0, 1, .. ket_rank - 1}
    std::iota(ket_indices.begin(), ket_indices.end(), 0);

    // generates a string annotation
    // input: vector<size_t>{10, 14, 19}
    // output:             "10,14,19"
    auto ords_to_csv_str = [](const container::svector<size_t> ords) {
      std::string str = "";
      for (auto ii : ords) {
        std::to_string(ii) + ",";
      }
      str.pop_back();  // remove the trailing comma ","
      return str;
    };

    // combine index vectors for bra and ket so that
    // a vector<size_t> can be generated for annotating
    // the whole tensor. Of course, the resulting vector<size_t>
    // has to be converted to string using 'ords_to_csv_str'
    // input: {1, 2, 0}, {0, 1, 2}
    // output: {1, 2, 0, 3, 4, 5}
    auto combine_ords = [](const std::vector<size_t>& ords1,
                           const std::vector<size_t>& ords2) {
      std::vector<size_t> combined(ords1.size() + ords2.size());
      std::copy(ords1.begin(), ords1.end(), combined.begin());
      std::transform(ords2.begin(), ords2.end(),
                     combined.begin() + ords1.size(),
                     [&ords1](size_t x) { return x + ords1.size(); });
      return combined;
    };

    DataTensorType result(ta_tensor.world(), ta_tensor.trange());
    result.fill(0.);

    // lhs_annot is always result of
    // ords_to_csv_str( 0, 1, ..., ta_tensor.rank()-1 )
    // ie "0,1,2,...ta_tensor.rank()-1"
    auto lhs_annot = ords_to_csv_str(combine_ords(bra_indices, ket_indices));

    // iter through the permutations of bra
    for (const auto& bp : _phase_perm(bra_indices)) {
      // iter through the permutations of ket
      for (const auto& kp : _phase_perm(ket_indices)) {
        // bra + ket permutation as a whole is even or odd?
        auto phase = std::get<0>(bp) * std::get<0>(kp);
        // regular TA scaling operation
        result(lhs_annot) = phase * ta_tensor(ords_to_csv_str(combine_ords(
                                        std::get<1>(bp), std::get<1>(kp))));
      }
    }
    // done antisymmetrizing
    return result;
  }

  ///
  /// Compute permutation of ordinals with phase.
  /// Phase is 1 if the the permutation is even -1 if it's odd.
  /// \param ords The ordinals to compute permutations on.
  /// \param begin The position index on @c ords to start computation
  /// \param swaps_count Number of swaps peformed before func call.
  /// \return vector of all permutations of ords
  /// TODO use coroutines to flatten the function
  static container::svector<std::tuple<int, container::svector<size_t>>>
  _phase_perm(container::svector<size_t>& ords, size_t begin = 0,
              size_t swaps_count = 0) {
    if (begin + 1 == ords.size()) {
      // found a new permutation
      // even permutation phase: +1 odd: -1.
      int phase = (swaps_count % 2 == 0) ? 1 : -1;

      return container::svector<std::tuple<int, container::svector<size_t>>>{
          std::make_tuple(phase, ords)};
    }

    // recursively call the function to compute the permutations
    // @note swaps_count is incremented only if the swapped indices are not
    // the same ie. to swap an element with itself is not counted.
    container::svector<std::tuple<int, container::svector<size_t>>> result;
    for (auto ii = begin; ii < ords.size(); ++ii) {
      std::swap(ords[begin], ords[ii]);

      auto recursed_result = _phase_perm(
          ords, begin + 1, ii == begin ? swaps_count : swaps_count + 1);

      for (auto& p : recursed_result) result.push_back(p);

      // undo the swap so the side effect is nullified for next call
      std::swap(ords[begin], ords[ii]);
    }
    return result;
  }
};  // class EvalNodeInternal

/// Evaluation node leaf.
template <typename DataTensorType>
class EvalNodeLeaf : public EvalNode<DataTensorType> {
 private:
  /// sequant ExprPtr to sequant Tensor corresponding to this leaf.
  ExprPtr expr_{nullptr};

  /// Swapped bra-ket indices?
  bool swapped_bra_ket_{false};

  void hash_node() override {
    //
    auto& tnsr = expr_->as<Tensor>();
    boost::hash<std::wstring_view> label_hasher;
    HashType leaf_hash_value;

    leaf_hash_value = label_hasher(tnsr.label());

    // cast IndexSpace attribute to size_t
    using ispace_cast = size_t;

    boost::hash<ispace_cast> number_hasher;

    using index_container_type =
        decltype(tnsr.bra());  // const svector<Index, ..>

    // hash and combine index space type from an Index container
    auto hash_idx_space = [&leaf_hash_value,
                           &number_hasher](index_container_type& idx_vec) {
      for (const auto& idx : idx_vec)
        boost::hash_combine(
            leaf_hash_value,
            number_hasher(static_cast<ispace_cast>(idx.space().type())));
    };

    if (swapped_bra_ket_) {
      hash_idx_space(tnsr.ket());
      hash_idx_space(tnsr.bra());
    } else {
      hash_idx_space(tnsr.bra());
      hash_idx_space(tnsr.ket());
    }

    // done hashing
    this->hash_value_ = leaf_hash_value;
    return;
  }

 public:
  /// Return true for a leaf node.
  bool is_leaf() const override { return true; }

  /// Operation count for leaf node is zero.
  OpsCount ops_count(const container::map<IndexSpace::TypeAttr, size_t>&
                         ispace_size_map) const override {
    return 0;
  }

  /// Visit node.
  void visit(const std::function<void(const EvalNode<DataTensorType>&)>&
                 visitor) const override {
    visitor(*this);
  }

  /// Evaluate the node.
  DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const override {
    if (auto label = expr_->as<Tensor>().label();
        (label == L"A" || label == L"P"))
      throw std::logic_error(
          "(anti-)symmetrization tensors cannot be evaluated from here!");

    auto found_it = context.find(this->hash_value());
    if (found_it != context.end()) {
      return *(found_it->second);
    }

    std::wstring error_msg_os = L"";

    error_msg_os += L"EvalNodeLeaf::evaluate(): ";
    error_msg_os += L"did not find such tensor in context (expr=\"";
    error_msg_os += expr_->as<Tensor>().to_latex() + L"\")";

    throw std::logic_error(
        std::string(error_msg_os.begin(), error_msg_os.end()));
  }

  /// Construct leaf node.
  /// \param tensor_expr sequant Expr to sequant Tensor
  EvalNodeLeaf(const ExprPtr& tensor_expr) : expr_{tensor_expr} {
    //
    auto& tnsr = expr_->as<Tensor>();
    swapped_bra_ket_ = need_bra_ket_swap();
    // set indices_
    if (swapped_bra_ket_) {
      for (const auto& idx : tnsr.ket()) this->indices_.emplace_back(idx);
      for (const auto& idx : tnsr.bra()) this->indices_.emplace_back(idx);
    } else {
      for (const auto& idx : tnsr.bra()) this->indices_.emplace_back(idx);
      for (const auto& idx : tnsr.ket()) this->indices_.emplace_back(idx);
    }
    // set hash_value_
    hash_node();
  }

  std::shared_ptr<EvalNodeLeaf> make_eval_node_leaf(
      const ExprPtr& tensor_expr) {
    return std::make_shared<EvalNodeLeaf>(tensor_expr);
  }

 private:
  /// Should we swap the whole bra indices with the whole ket indices?
  ///
  /// When a sequant Tensor represents a real valued data tensor extra
  /// canonicalization is done.
  /// @return True if swapping is necessary.
  bool need_bra_ket_swap() {
    auto& tnsr = expr_->as<Tensor>();
    for (auto i = 0; i < (tnsr.bra_rank() < tnsr.ket_rank() ? tnsr.bra_rank()
                                                            : tnsr.ket_rank());
         ++i) {
      auto bindex_space = tnsr.bra().at(i).space();
      auto kindex_space = tnsr.ket().at(i).space();
      // same index space of both indices in corresponding positions in bra
      // and ket: continue
      if (bindex_space == kindex_space) continue;

      // if index space of ket index is smaller than that of bra index
      // swapping is necessary
      return kindex_space < bindex_space;
    }

    // survived this far implies no swapping necessary
    return false;
  }
};  // class EvalNodeLeaf

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVALUATE_EVAL_NODE_HPP
