#ifndef SEQUANT_EVALUATE_EVAL_TREE_HPP
#define SEQUANT_EVALUATE_EVAL_TREE_HPP

#include "eval_fwd.hpp"
#include "eval_tree_node.hpp"

#include <SeQuant/core/tensor.hpp>

#include <functional>
#include <memory>
#include <numeric>

namespace sequant::evaluate {
///
/// \brief Representation of binary evaluation of sequant::Expr objects.
///
/// The SeQuant Expr(s) can be thought of either as a representation of a
/// single data-tensor (eg. btas::Tensor<double>, TA::TArrayD), or binary
/// evaluations of such data-tensors. An evaluation could be summing two
/// data-tensors, taking a product of them antisymmetrizing tensors. A result of
/// a binary evaluation can be an input for another evaluation, i.e. evaluations
/// can be nested.
///
/// \author Bimal Gaudel
/// \date Apr 21, 2020
///
class EvalTree {
 private:
  /// The root node of the evaluation tree.
  EvalNodePtr root{nullptr};

 public:
  /// Return the hash value of the tree.
  HashType hash_value();

  /// Compute operation counts for the evaluation.
  /// \param ispace_size_map A map from IndexSpace type to the size of the
  ///                       space.
  ///                       e.g. IndexSpace::active_occupied -> 5
  ///                       e.g. IndexSpace::active_unoccupied -> 20
  OpsCount ops_count(
      const container::map<IndexSpace::TypeAttr, size_t>& ispace_size_map);

  /// Evaluate the tree in a given context.
  /// \tparam DataTensorType Type of the backend data tensor. eg. TA::TArrayD
  /// while using TiledArray
  /// \param  context A map that maps hash values of (at least) all the leaf
  /// nodes in the tree to the DataTensorType tensor.
  /// \return Result of evaluation that is of DataTensorType.
  template <typename DataTensorType>
  DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const {
    return _evaluate(root, context);
  }

  /// \brief Evaluate: Use lambda to generate unknown Tensor objects
  /// \param eval_tensor External lambda function that generates tensor
  /// \return Result after evaluating all tensor operations
  template <typename DataTensorType>
  DataTensorType evaluate_and_make(
      const std::function<DataTensorType(const Tensor&)>& eval_tensor)
  const {
    return _evaluate_and_make(root, eval_tensor);
  }

  template<typename DataTensorType>
  DataTensorType symmetrize(DataTensorType& tensor, int rank){
    return _symmetrize(tensor, rank);
  }

  /// Visit the tree by pre-order traversal.
  /// ie. Node is visited first followed by left node and then right node.
  void visit(const std::function<void(const EvalNodePtr&)>& visitor);

  /// Swap labels of leaf nodes made out of such sequant tensor object whose
  /// to_latex() result matches that of @c expr.
  bool swap_labels(const ExprPtr& expr);

  /// Construct eval tree from sequant expression.
  ///
  /// \param canonize_leaf_braket When true (default) swaps bra indices and with
  ///        the ket indices and returns while constructing leaf node from
  ///        sequant tensor as soon as IndexSpace attribute in ket is lower than
  ///        that in bra at corresponding positions.
  explicit EvalTree(const ExprPtr& expr, bool canonize_leaf_braket = false);

  /// \brief Antisymmetrizes a data tensor with same bra/ket ranks
  /// \details The function performs the action of Antisymmetrizer operator on a tensor generating
  /// a sum of all (n!)^2 permutations where n is the bra/ket rank.
  /// e.g: A_{ab}^{ij}t_{ij}^{ab} = t_{ij}^{ab} - t_{ji}^{ab} - t_{ij}^{ba} + t_{ji}^{ba}
  /// \param ta_tensor Data tensor to be antisymmetrized
  /// \param bra_rank
  /// \param ket_rank
  /// \result the antisymmetrized data tensor
  template <typename DataTensorType>
  static DataTensorType antisymmetrize_tensor(
      const DataTensorType& ta_tensor,
      size_t bra_rank,
      size_t ket_rank){
    return _antisymmetrize(ta_tensor, bra_rank, ket_rank);
  }

 private:
  /// Build EvalTreeNode pointer from sequant expression of Sum, Product or
  /// Tensor type.
  static EvalNodePtr build_expr(const ExprPtr& expr, bool canonize_leaf_braket);

  /// Build EvalTreeNode pointer from sequant expression of Sum type.
  static EvalNodePtr build_sum(const ExprPtr& expr, bool canonize_leaf_braket);

  /// Build EvalTreeNode pointer from sequant expression of Product type.
  static EvalNodePtr build_prod(const ExprPtr& expr, bool canonize_leaf_braket);

  /// visit each node
  static void _visit(const EvalNodePtr& node,
                     const std::function<void(const EvalNodePtr&)>& visitor);

  /// Get operation counts of a node.
  static OpsCount _ops_count(
      const EvalNodePtr& node,
      const container::map<IndexSpace::TypeAttr, size_t>& ispace_size_map);

  /// Swap bra ket labels of those leaf tensors for which the @c predicate
  /// function evaluates true.
  static bool _swap_braket_labels(
      const EvalNodePtr& node,
      const std::function<bool(const EvalTreeLeafNode&)>& predicate);

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
              size_t swaps_count = 0);

  /// \brief Antisymmetrize DataTensorType
  /// \details The function performs the action of Antisymmetrizer operator on a tensor generating
  /// a sum of all (n!)^2 permutations where n is the bra/ket rank.
  /// e.g: A_{ab}^{ij}t_{ij}^{ab} = t_{ij}^{ab} - t_{ji}^{ab} - t_{ij}^{ba} + t_{ji}^{ba}
  /// \tparam DataTensorType Backend data tensor type eg. TA::TArrayD
  /// \param ta_tensor TiledArray tensor.
  /// \param bra_rank rank of the tensor bra.
  /// \param ket_rank rank of the tensor ket.
  /// \param scal ScalarType factor to scale the result. 1 by default.
  /// \return An antisymmetrized data tensor
  template <typename DataTensorType>
  static DataTensorType _antisymmetrize(const DataTensorType& ta_tensor,
                                        size_t bra_rank, size_t ket_rank,
                                        ScalarType scal = 1);

  /// \brief Symmetrize DataTensorType
  /// \details The function performs the action of Symmetrizer on a tensor generating
  /// a sum of all n! permutations where n is the bra/ket rank.
  /// e.g: S_{ab}^{ij}t_{ij}^{ab} = t_{ij}^{ab} + t_{ji}^{ba}
  /// \tparam DataTensorType Backend data tensor type eg. TA::TArrayD
  /// \param ta_tensor TiledArray tensor.
  /// \param braket_rank rank of the tensor bra.
  /// \param scal ScalarType factor to scale the result. 1 by default.
  /// \return A symmetrized data tensor
  template <typename DataTensorType>
  static DataTensorType _symmetrize(const DataTensorType& ta_tensor,
                                    size_t braket_rank, ScalarType scal = 1);

  /// Evaluate the tree in a given context.
  /// \tparam DataTensorType Type of the backend data tensor. eg. TA::TArrayD
  /// while using TiledArray
  /// \param  context A map that maps hash values of (at least) all the leaf
  /// nodes in the tree to the DataTensorType tensor.
  /// \return Result of evaluation that is of DataTensorType.
  template <typename DataTensorType>
  static DataTensorType _evaluate(
      const EvalNodePtr& node,
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context);

  /// \brief Evaluate a node in absence of a context
  /// \details If no context is provided, use a lambda function to generate
  /// a data tensor corresponding to every sequant::Tensor in the node
  /// \param node root node of evaluation tree
  /// \param eval_tensor A lambda function that generates a data tensor
  /// \return Result after evaluating all tensor operations
  template <typename DataTensorType>
  static DataTensorType _evaluate_and_make(
      const EvalNodePtr& node,
      const std::function<DataTensorType(const Tensor&)>& eval_tensor);

};  // class EvalTree

/// \brief Antisymmetrize DataTensorType
/// \tparam DataTensorType Backend data tensor type eg. TA::TArrayD
/// \param ta_tensor TiledArray tensor.
/// \param bra_rank rank of the tensor bra.
/// \param ket_rank rank of the tensor ket.
/// \param scal ScalarType factor to scale the result. 1 by default.
/// \return An antisymmetrized data tensor
template <typename DataTensorType>
DataTensorType EvalTree::_antisymmetrize(const DataTensorType& ta_tensor,
                                         size_t bra_rank, size_t ket_rank,
                                         ScalarType scal) {
  using ordinal_indices = container::svector<size_t>;
  ordinal_indices bra_indices(bra_rank), ket_indices(ket_rank);
  // {0, 1, .. bra_rank - 1}
  std::iota(bra_indices.begin(), bra_indices.end(), 0);
  // {0, 1, .. ket_rank - 1}
  std::iota(ket_indices.begin(), ket_indices.end(), 0);

  // generates a string annotation
  // input: vector<size_t>{10, 14, 19}
  // output:             "10,14,19"
  auto ords_to_csv_str = [](const auto& ords) {
    std::string str;
    for (auto ii : ords) {
      str += std::to_string(ii) + ",";
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
  auto combine_ords = [](const ordinal_indices& ords1,
                         const ordinal_indices& ords2) {
    ordinal_indices combined(ords1.size() + ords2.size());
    std::copy(ords1.begin(), ords1.end(), combined.begin());
    std::transform(ords2.begin(), ords2.end(), combined.begin() + ords1.size(),
                   [&ords1](size_t x) { return x + ords1.size(); });
    return combined;
  };

  DataTensorType result(ta_tensor.world(), ta_tensor.trange());
  result.fill(0);

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

      auto rhs_annot =
          ords_to_csv_str(combine_ords(std::get<1>(bp), std::get<1>(kp)));
      // regular TA scaling operation
      result(lhs_annot) += phase * ta_tensor(rhs_annot);
    }
  }
  // done antisymmetrizing

  // scale if needed
  if (scal != 1) result(lhs_annot) = scal * result(lhs_annot);

  return result;

}  // function _antisymmetrize
/// \brief Symmetrize DataTensorType
/// \tparam DataTensorType Backend data tensor type eg. TA::TArrayD
/// \param ta_tensor TiledArray tensor.
/// \param braket_rank rank of the tensor bra.
/// \param scal ScalarType factor to scale the result. 1 by default.
/// \return A symmetrized data tensor
template <typename DataTensorType>
DataTensorType EvalTree::_symmetrize(const DataTensorType& ta_tensor,
                                     size_t braket_rank, ScalarType scal) {
  // generates a string annotation
  // input: vector<size_t>{10, 14, 19}
  // output:             "10,14,19"
  auto ords_to_csv_str = [](const auto& ords, auto add) {
    std::string str;
    for (auto ii : ords) {
      str += std::to_string(ii + add) + ",";
    }
    str.pop_back();  // remove the trailing comma ","
    return str;
  };

  DataTensorType result(ta_tensor.world(), ta_tensor.trange());
  result.fill(0);

  // generate a vector and fill it with int starting with 0
  auto perm_vec = container::svector<size_t>(braket_rank);
  std::iota(perm_vec.begin(), perm_vec.end(), 0);

  // fix indices of the result tensor
  auto lhs_annot = ords_to_csv_str(perm_vec, 0) + "," +
                   ords_to_csv_str(perm_vec, braket_rank);

  // keep adding the permuted tensors to the result
  do {
    auto rhs_annot = ords_to_csv_str(perm_vec, 0) + "," +
                     ords_to_csv_str(perm_vec, braket_rank);
    result(lhs_annot) += ta_tensor(rhs_annot);

  } while (std::next_permutation(perm_vec.begin(), perm_vec.end()));

  result(lhs_annot) = scal * result(lhs_annot);
  return result;
}

/// Evaluate the tree in a given context.
/// \param node Root node of evaluation tree
/// \param context A map that maps hash values of (at least) all the leaf
/// nodes in the tree to the DataTensorType tensor.
/// \return Result of evaluation that is of DataTensorType.
template <typename DataTensorType>
DataTensorType EvalTree::_evaluate(
    const EvalNodePtr& node,
    const container::map<HashType, std::shared_ptr<DataTensorType>>& context) {
  if (node->is_leaf()) {
    auto leaf_node = std::dynamic_pointer_cast<EvalTreeLeafNode>(node);
    if (auto label = leaf_node->expr()->as<Tensor>().label();
        (label == L"A" || label == L"S")) {
      throw std::logic_error(
          "(anti-)symmetrization tensors cannot be evaluated from here!");
    }

    auto found_it = context.find(node->hash_value());
    if (found_it != context.end()) {
      return *(found_it->second);
    } else {
      std::wstring error_msg_os;

      error_msg_os += L"EvalNodeLeaf::evaluate(): ";
      error_msg_os += L"did not find such tensor in context (expr=\"";
      error_msg_os += leaf_node->expr()->as<Tensor>().to_latex() + L"\")";

      throw std::logic_error(
          std::string(error_msg_os.begin(), error_msg_os.end()));
    }
  }  // done leaf evaluation

  //
  // non-leaf evaluation
  //
  auto intrnl_node = std::dynamic_pointer_cast<EvalTreeInternalNode>(node);

  auto opr = intrnl_node->operation();
  if (opr == Operation::ANTISYMMETRIZE) {
    auto bra_rank = intrnl_node->indices().size() / 2;
    auto ket_rank = intrnl_node->indices().size() - bra_rank;
    return _antisymmetrize(_evaluate(intrnl_node->right(), context), bra_rank,
                           ket_rank, intrnl_node->right()->scalar());
  }  // anitsymmetrization type evaluation done

  if (opr == Operation::SYMMETRIZE) {
    auto braket_rank = intrnl_node->indices().size();
    if (braket_rank % 2 != 0)
      throw std::logic_error("Can not symmetrize odd-ordered tensor!");

    braket_rank /= 2;
    return _symmetrize(_evaluate(intrnl_node->right(), context), braket_rank,
                       intrnl_node->right()->scalar());
  }  // symmetrization type evaluation done

  // generates tiledarray annotation based on a node's index labels
  // @note this wouldn't be necessary if the tensor algebra library
  // would support std::string_view as annotations
  auto TA_annotation =
      [&intrnl_node](decltype(intrnl_node->indices())& indices) {
        std::string annot;
        for (const auto& idx : indices)
          annot += std::string(idx.label().begin(), idx.label().end()) + ", ";

        annot.erase(annot.size() - 2);  // remove trailing ", "
        return annot;
      };

  auto left_annot = TA_annotation(intrnl_node->left()->indices());
  auto right_annot = TA_annotation(intrnl_node->right()->indices());
  auto this_annot = TA_annotation(intrnl_node->indices());

  DataTensorType result;
  if (opr == Operation::SUM) {
    // sum left and right evaluated tensors
    // using tiled array syntax
    result(this_annot) =
        intrnl_node->left()->scalar() *
            _evaluate(intrnl_node->left(), context)(left_annot) +
        intrnl_node->right()->scalar() *
            _evaluate(intrnl_node->right(), context)(right_annot);
    //
  } else if (opr == Operation::PRODUCT) {
    // contract left and right evaluated tensors
    // using tiled array syntax
    result(this_annot) = intrnl_node->left()->scalar() *
                         _evaluate(intrnl_node->left(), context)(left_annot) *
                         intrnl_node->right()->scalar() *
                         _evaluate(intrnl_node->right(), context)(right_annot);
  } else {
    throw std::domain_error("Operation: " + std::to_string((size_t)opr) +
                            " not supported!");
  }  // sum and product type evaluation
  result.world().gop.fence();
  return result;

}  // function _evaluate

/// @brief evaluates a tree, uses a tensor generator lambda if the tensor object is missing from context
/// @param node Root node of evaluation tree
/// @param context Mutable map of hash values and Tensors
/// @param eval_tensor lambda function that returns a Tensor object
/// @return Result after evaluating all tensor operations
template <typename DataTensorType>
DataTensorType EvalTree::_evaluate_and_make(
    const EvalNodePtr& node,
    const std::function<DataTensorType(const Tensor&)>& eval_tensor) {

  // If node is a leaf, return the object from context
  if (node->is_leaf()) {
    auto leaf_node = std::dynamic_pointer_cast<EvalTreeLeafNode>(node);
    if (auto label = leaf_node->expr()->as<Tensor>().label();
        (label == L"A" || label == L"S")) {
      throw std::logic_error(
          "(anti-)symmetrization tensors cannot be evaluated from here!");
    }
    return eval_tensor(leaf_node->expr()->as<Tensor>());
  }  // done leaf evaluation

  //
  // non-leaf evaluation
  //
  auto intrnl_node = std::dynamic_pointer_cast<EvalTreeInternalNode>(node);

  auto opr = intrnl_node->operation();
  if (opr == Operation::ANTISYMMETRIZE) {
    auto bra_rank = intrnl_node->indices().size() / 2;
    auto ket_rank = intrnl_node->indices().size() - bra_rank;
    return _antisymmetrize(_evaluate_and_make(intrnl_node->right(), eval_tensor), bra_rank,
                           ket_rank, intrnl_node->right()->scalar());
  }  // anitsymmetrization type evaluation done

  if (opr == Operation::SYMMETRIZE) {
    auto braket_rank = intrnl_node->indices().size();
    if (braket_rank % 2 != 0)
      throw std::logic_error("Can not symmetrize odd-ordered tensor!");

    braket_rank /= 2;
    return _symmetrize(_evaluate_and_make(intrnl_node->right(), eval_tensor), braket_rank,
                       intrnl_node->right()->scalar());
  }  // symmetrization type evaluation done

  // generates tiledarray annotation based on a node's index labels
  // @note this wouldn't be necessary if the tensor algebra library
  // would support std::string_view as annotations
  auto TA_annotation =
      [](decltype(intrnl_node->indices())& indices) {
        std::string annot;
        for (const auto& idx : indices)
          annot += std::string(idx.label().begin(), idx.label().end()) + ", ";

        annot.erase(annot.size() - 2);  // remove trailing ", "
        return annot;
      };

  auto left_annot = TA_annotation(intrnl_node->left()->indices());
  auto right_annot = TA_annotation(intrnl_node->right()->indices());
  auto this_annot = TA_annotation(intrnl_node->indices());

  DataTensorType result;
  if (opr == Operation::SUM) {
    // Evaluate terms as intermediates
    result(this_annot) =
        intrnl_node->left()->scalar() *
            _evaluate_and_make(intrnl_node->left(), eval_tensor)(left_annot);

     result(this_annot) = result(this_annot) +
            intrnl_node->right()->scalar() *
                _evaluate_and_make(intrnl_node->right(), eval_tensor)(right_annot);
  } else if (opr == Operation::PRODUCT) {
    // contract left and right evaluated tensors using TA syntax
    result(this_annot) = intrnl_node->left()->scalar() *
        _evaluate_and_make(intrnl_node->left(), eval_tensor)(left_annot) *
        intrnl_node->right()->scalar() *
        _evaluate_and_make(intrnl_node->right(), eval_tensor)(right_annot);
  } else {
    throw std::domain_error("Operation: " + std::to_string((size_t)opr) +
        " not supported!");
  }  // sum and product type evaluation
  result.world().gop.fence();
  return result;

}  // function _evaluate_and_make

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVALUATE_EVAL_TREE_HPP
