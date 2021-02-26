#include "eval_tree_node.hpp"
#include "eval_fwd.hpp"

#include <SeQuant/core/tensor.hpp>

#include <algorithm>

namespace sequant::evaluate {

const IndexContainer& sequant::evaluate::EvalTreeNode::indices() const {
  return indices_;
}

HashType EvalTreeNode::hash_value() const { return hash_value_; }

ScalarType EvalTreeNode::scalar() const { return scalar_; }

void EvalTreeNode::scale(ScalarType scal) { scalar_ = scal; }

void EvalTreeNode::update_hash() { this->hash_value_ = this->hash_node(); }

const EvalNodePtr& EvalTreeInternalNode::left() const { return left_; }

const EvalNodePtr& EvalTreeInternalNode::right() const { return right_; }

Operation EvalTreeInternalNode::operation() const { return operation_; }

bool EvalTreeInternalNode::is_leaf() const { return false; }

EvalTreeInternalNode::EvalTreeInternalNode(const EvalNodePtr& left_ptr,
                                           const EvalNodePtr& right_ptr,
                                           Operation op)
    : left_{left_ptr}, right_{right_ptr}, operation_{op} {
  auto& left = this->left();
  auto& right = this->right();

  if (op == Operation::SUM) {  // building sum type evaluation
    if (left->indices().size() != right->indices().size())
      throw std::logic_error(
          "Not matching number of indices for summation type evaluation.");
    // indices for the intermediate
    this->indices_ = !left->is_leaf() ? left->indices()
                                      : !right->is_leaf() ? right->indices()
                                                          : IndexContainer{};
    // couldn't decide on the indices of the result because
    // both nodes were leaves and leaf indices may not be sorted
    if (this->indices().empty()) {
      this->indices_ = left->indices();  // or right()->indices()
      std::sort(this->indices_.begin(), this->indices_.end());
    }
  } else if (op == Operation::PRODUCT) {  // building product type evaluation
    auto lindices = left->indices();
    auto rindices = right->indices();
    // indices from leaves may not be sorted
    if (left->is_leaf()) std::sort(lindices.begin(), lindices.end());
    if (right->is_leaf()) std::sort(rindices.begin(), rindices.end());
    // set the non-contracting indices
    std::set_symmetric_difference(lindices.begin(), lindices.end(),
                                  rindices.begin(), rindices.end(),
                                  std::back_inserter(this->indices_));
  } else if ((op == Operation::ANTISYMMETRIZE) ||
             (op == Operation::SYMMETRIZE)) {
    this->indices_ = left->indices();
    std::sort(this->indices_.begin(), this->indices_.end());
  } else {
    throw std::logic_error("Invalid Operation while forming intermediate");
  }

  // set the hash value
  this->update_hash();
}

HashType EvalTreeInternalNode::hash_node() const {
  auto imed_operation = this->operation();

  // number_hasher is used to hash index positions and operation type
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
    return imed_hash_value;
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
    auto index_hash_combiner = [&number_hasher, &imed_hash_value,
                                this](const EvalNodePtr& lr_node) {
      // iterate through this node's indices
      for (const auto& idx : this->indices()) {
        // iterate through the indices of the passed left/right node
        for (auto i = 0; i < lr_node->indices().size(); ++i) {
          if (idx.label() == lr_node->indices()[i].label()) {
            // if the indexes match then time to break as
            // we found the contracted indices
            // no need to look further for this same label
            // in the left/right node
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

  return imed_hash_value;
}

std::wstring EvalTreeInternalNode::to_latex() const {
    std::wostringstream tex;

    tex << L"I_{"; // All intermediates have the tensor label I

    for (auto ii = 0;  ii < indices().size()/2; ++ii)
        tex << indices().at(ii).to_latex();
    tex << L"}^{";

    for (auto ii = indices().size()/2; ii < indices().size(); ++ii)
        tex << indices().at(ii).to_latex();
    tex << L"}";

    return tex.str();
}

    bool EvalTreeLeafNode::is_leaf() const { return true; }

void EvalTreeLeafNode::swap_labels() {
  // based on the previous state
  // toggle the state of the object
  // to indicate whether the bra-ket labels
  // should be swapped or unswapped
  swapped_labels_ = !swapped_labels_;
  // then set labels
  set_labels();
  this->update_hash();
}

void EvalTreeLeafNode::set_labels() {
  this->indices_.clear();
  if (auto& tnsr = expr()->as<Tensor>(); swapped_labels_) {
    // if bra-ket swapped state
    // set indices from ket first followed by bra
    for (const auto& idx : tnsr.ket()) this->indices_.emplace_back(idx);
    for (const auto& idx : tnsr.bra()) this->indices_.emplace_back(idx);
  } else {
    // else set from bra first followed by ket
    for (const auto& idx : tnsr.bra()) this->indices_.emplace_back(idx);
    for (const auto& idx : tnsr.ket()) this->indices_.emplace_back(idx);
  }
}

const ExprPtr& EvalTreeLeafNode::expr() const { return expr_; }

bool EvalTreeLeafNode::swapped_labels() const { return swapped_labels_; }

bool EvalTreeLeafNode::canonize_swap(const ExprPtr& expr) {
  auto& tnsr = expr->as<Tensor>();
  for (auto i = 0; i < (tnsr.bra_rank() < tnsr.ket_rank() ? tnsr.bra_rank()
                                                          : tnsr.ket_rank());
       ++i) {  // iterate through index positions in bra and ket
    auto bindex_space = tnsr.bra().at(i).space();
    auto kindex_space = tnsr.ket().at(i).space();
    // same index space of both indices in corresponding positions
    // of bra and ket: continue
    if (bindex_space == kindex_space) continue;

    // if index space of ket index is smaller than that of bra index
    // swapping is to be done
    return kindex_space < bindex_space;
  }

  // survived this far implies no swapping necessary
  return false;
}

HashType EvalTreeLeafNode::hash_node() const {
  auto& tnsr = expr()->as<Tensor>();

  HashType leaf_hash_value;

  boost::hash<decltype(tnsr.bra().begin()->label())> label_hasher;
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

  if (swapped_labels_) {
    hash_idx_space(tnsr.ket());
    hash_idx_space(tnsr.bra());
  } else {
    hash_idx_space(tnsr.bra());
    hash_idx_space(tnsr.ket());
  }
  // done hashing
  return leaf_hash_value;
}

EvalTreeLeafNode::EvalTreeLeafNode(const ExprPtr& tnsr_expr,
                                   bool canonize_braket)
    : expr_{tnsr_expr} {
  swapped_labels_ = canonize_braket && canonize_swap(expr());

  set_labels();

  this->update_hash();
}

std::wstring EvalTreeLeafNode::to_latex() const {
        return expr()->to_latex();
}

}  // namespace sequant::evaluate
