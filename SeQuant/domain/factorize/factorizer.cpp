#include "factorizer.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_network.hpp>
// IndexSpace type based hashing of tensors for AdjacencyMatrix
#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <iomanip>
#include <ios>
#include <tuple>

#define DEBUG_PRINT 0

namespace sequant::factorize {

using pos_type = AdjacencyMatrix::pos_type;
using color_mat_type = AdjacencyMatrix::color_mat_type;

AdjacencyMatrix::AdjacencyMatrix(
    const ExprPtr& expr,
    const TensorNetwork::named_indices_t& external_indices)
    : colorMatrix_(  // allocates color matrix
          expr->size(),
          color_mat_type::value_type(
              expr->size(), color_mat_type::value_type::value_type{})) {
  // get a Product out of a TensorNetwork
  auto tn_to_prod = [](const auto& tn) {
    auto prod = std::make_shared<Product>();
    for (const auto& tnsr : tn.tensors())
      prod->append(1, std::dynamic_pointer_cast<Tensor>(tnsr));
    return prod;
  };

  // fill the color values
  for (auto ii = 0; ii < expr->size(); ++ii)
    for (auto jj = ii + 1; jj < expr->size(); ++jj) {
      // set color data
      if (are_connected(expr->at(ii), expr->at(jj))) {
        // auto prod = std::make_shared<Product>(Product{expr->at(ii),
        // expr->at(jj)});
        auto tnet =
            TensorNetwork(*(expr->at(ii)->clone() * expr->at(jj)->clone()));
        tnet.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false,
                          &external_indices);
        colorMatrix_[ii][jj] = colorMatrix_[jj][ii] =
            evaluate::EvalTree(tn_to_prod(tnet)).hash_value();
      }
    }
}

AdjacencyMatrix::AdjacencyMatrix(
    const container::svector<ExprPtr>& tensors,
    const TensorNetwork::named_indices_t& external_indices)
    : AdjacencyMatrix(
          std::make_shared<Product>(1, tensors.begin(), tensors.end()),
          external_indices) {}

bool AdjacencyMatrix::are_connected(const ExprPtr& t1, const ExprPtr& t2) {
  auto tnsr1 = t1->as<Tensor>();
  auto tnsr2 = t2->as<Tensor>();
  // iterate through the bra and the ket labels of tnsr1 and tnsr2
  // if any index label is common, they are connected.
  for (const auto& idx1 : tnsr1.const_braket())
    for (const auto& idx2 : tnsr2.const_braket()) {
      if (idx1.label() == idx2.label()) return true;
    }
  return false;
}

const color_mat_type& AdjacencyMatrix::color_mat() const {
  return colorMatrix_;
}

size_t AdjacencyMatrix::num_verts() const { return color_mat().size(); }

bool AdjacencyMatrix::are_connected(pos_type pos1, pos_type pos2) const {
  return color_mat()[pos1][pos2] != color_mat_type::value_type::value_type{};
}

color_mat_type::value_type::value_type AdjacencyMatrix::color(
    pos_type pos1, pos_type pos2) const {
  return color_mat()[pos1][pos2];
}

// expr is Product type
// tnsr is Tensor type
bool tensor_exists(const ExprPtr& expr, const ExprPtr& tnsr) {
  if (expr->is<Tensor>()) {
    return *expr == *tnsr;
  }

  for (const auto& xpr : *expr)
    if (tensor_exists(xpr, tnsr)) return true;
  return false;
};

/***********************************
 * Functions for factorization     *
 ***********************************/

// get target indices
TensorNetwork::named_indices_t target_indices(const ExprPtr& expr) {
  TensorNetwork::named_indices_t result;
  for (const auto& tnsr : *expr) {
    for (const auto& idx : tnsr->as<Tensor>().const_braket()) {
      if (result.contains(idx))
        result.erase(idx);
      else
        result.insert(idx);
    }
  }
  return result;
}

// Get positions of common type of tensors in a pair of Exprs'.
std::tuple<container::set<pos_type>, container::set<pos_type>> common_tensors(
    const ExprPtr& expr1, const ExprPtr& expr2) {
  container::set<pos_type> commonT_1, commonT_2;
  for (pos_type ii = 0; ii < expr1->size(); ++ii)
    for (pos_type jj = 0; jj < expr2->size(); ++jj) {
      //
      // NOTE: As an example: t_{i j}^{a b} has the same
      // hash value as t_{a b}^{i j}. To hash such expressions
      // differently, use EvalTree(expr, false).
      //
      if (evaluate::EvalTree(expr1->at(ii)).hash_value() ==
          evaluate::EvalTree(expr2->at(jj)).hash_value()) {
        commonT_1.insert(ii);
        commonT_2.insert(jj);
      }
    }
  return std::tuple(commonT_1, commonT_2);
}

container::map<std::tuple<pos_type, pos_type>, std::tuple<pos_type, pos_type>>
common_pairs(const AdjacencyMatrix& mat1, const AdjacencyMatrix& mat2) {
  auto result = container::map<std::tuple<pos_type, pos_type>,
                               std::tuple<pos_type, pos_type>>{};
  for (auto ii = 0; ii < mat1.num_verts(); ++ii)
    for (auto jj = ii + 1; jj < mat1.num_verts(); ++jj)
      for (auto kk = 0; kk < mat2.num_verts(); ++kk)
        for (auto ll = kk + 1; ll < mat2.num_verts(); ++ll)
          if (auto color = mat1.color(ii, jj);
              color != color_mat_type::value_type::value_type{} &&
              color == mat2.color(kk, ll)) {
            result.insert(decltype(result)::value_type{std::tuple(ii, jj),
                                                       std::tuple(kk, ll)});
            break;
          }

  return result;
}

std::tuple<container::svector<container::set<AdjacencyMatrix::pos_type>>,
           container::svector<container::set<AdjacencyMatrix::pos_type>>>
common_nets(const container::map<std::tuple<pos_type, pos_type>,
                                 std::tuple<pos_type, pos_type>>& pairs) {
  // make a copy of pairs
  auto common_p = pairs;
  // to hold processed nets
  container::svector<container::set<AdjacencyMatrix::pos_type>> common_net1,
      common_net2;
  decltype(common_p)::iterator pair_iter;
  while (!common_p.empty()) {
    pair_iter = common_p.begin();

    decltype(common_net1)::value_type net1, net2;

    net1.insert(std::get<0>(pair_iter->first));
    net1.insert(std::get<1>(pair_iter->first));

    net2.insert(std::get<0>(pair_iter->second));
    net2.insert(std::get<1>(pair_iter->second));

    pair_iter = common_p.erase(pair_iter);
    while (pair_iter != common_p.end()) {
      auto [pFirst1, pFirst2] = pair_iter->first;
      auto [pSecond1, pSecond2] = pair_iter->second;

      if (net1.contains(pFirst1) || net1.contains(pFirst2)) {
        assert(net2.contains(pSecond1) || net2.contains(pSecond2));
        pair_iter = common_p.erase(pair_iter);
      } else {
        assert(!(net2.contains(pSecond1) || net2.contains(pSecond2)));
        ++pair_iter;
      }
    }
    common_net1.push_back(net1);
    common_net2.push_back(net2);
  }
  return std::tuple(common_net1, common_net2);
}

std::tuple<ExprPtr, ExprPtr> factorize_pair(const ExprPtr& expr1,
                                            const ExprPtr& expr2) {
  // get common type of tensor's positions
  auto [commonIdx1, commonIdx2] = common_tensors(expr1, expr2);

  container::svector<ExprPtr> commonExpr1, commonExpr2;
  for (auto idx : commonIdx1) {
    commonExpr1.push_back(expr1->at(idx));
  }
  for (auto idx : commonIdx2) {
    commonExpr2.push_back(expr2->at(idx));
  }

#ifndef DEBUG_PRINT
  std::wcout << "\nprinting common tensors\n"
             << "-----------------------\n";
  for (const auto& xpr : commonExpr1) std::wcout << xpr->to_latex() << " ";
  std::wcout << "\n";
  for (const auto& xpr : commonExpr2) std::wcout << xpr->to_latex() << " ";
  std::wcout << "\n";
#endif

  // get the target indices
  auto target1 = target_indices(expr1);
  auto target2 = target_indices(expr2);

#ifndef DEBUG_PRINT
  std::wcout << "\nprinting target indices\n"
             << "-----------------------\n";
  for (const auto& idx : target1) std::wcout << idx.to_latex() << " ";
  std::wcout << "\n";
  for (const auto& idx : target2) std::wcout << idx.to_latex() << " ";
  std::wcout << "\n";
#endif

  // form adjacency matrices for common tensors
  auto adjMat1 = AdjacencyMatrix(commonExpr1, target1);
  auto adjMat2 = AdjacencyMatrix(commonExpr2, target2);

#ifndef DEBUG_PRINT
  auto print_adj_mat = [](const auto& mat, bool color = false) {
    for (auto ii = 0; ii < mat.num_verts(); ++ii) {
      for (auto jj = 0; jj < mat.num_verts(); ++jj) {
        if (color)
          std::wcout << std::setw(27) << mat.color(ii, jj);
        else
          std::wcout << mat.are_connected(ii, jj);
        std::wcout << "  ";
      }
      std::wcout << "\n";
    }
  };

  std::wcout << "\nprinting adjacency matrix1\n"
             << "--------------------------\n";
  print_adj_mat(adjMat1, true);
  std::wcout << "\nprinting adjacency matrix2\n"
             << "--------------------------\n";
  print_adj_mat(adjMat2, true);
#endif

  // find common pairs
  auto common_p = common_pairs(adjMat1, adjMat2);

#ifndef DEBUG_PRINT
  std::wcout << "\nprinting common pair indices\n"
             << "----------------------------\n";
  for (const auto& p : common_p) {
    auto [tup1, tup2] = p;
    std::wcout << "(" << std::get<0>(tup1) << ", " << std::get<1>(tup1)
               << ")  ";
    std::wcout << "(" << std::get<0>(tup2) << ", " << std::get<1>(tup2) << ")";
    std::wcout << "\n";
  }
  std::wcout << std::endl;
#endif

  // find common nets
  auto [netIdx1, netIdx2] = common_nets(common_p);

  // form selective ExprPtr out of an iterable of
  // pos_types and a given reference Expr of size at least max(iterable.values)
  auto iter_to_expr = [](const auto& iterable, const auto& expr) {
    auto result = std::make_shared<Product>();
    for (auto idx : iterable) result->append(1, expr.at(idx));
    return result;
  };

  auto subnet1 = std::make_shared<Product>();
  for (const auto& group : netIdx1)
    subnet1->append(1, iter_to_expr(group, commonExpr1));

  auto subnet2 = std::make_shared<Product>();
  for (const auto& group : netIdx2)
    subnet2->append(1, iter_to_expr(group, commonExpr2));

  // append left out tensors from the original expr to individual subnets
  auto get_left_out = [](const auto& commonidx, const ExprPtr& factorized,
                         const ExprPtr& original) {
    auto left_out = std::make_shared<Product>();

    for (pos_type ii = 0; ii < original->size(); ++ii)
      if (auto tnsr = original->at(ii); !(tensor_exists(factorized, tnsr)))
        left_out->append(1, tnsr->clone());
    return left_out;
  };

  auto left1 = get_left_out(commonIdx1, subnet1, expr1);
  auto left2 = get_left_out(commonIdx2, subnet2, expr2);

#ifndef DEBUG_PRINT
  std::wcout << "subnet1 = " << subnet1->to_latex() << "\n"
             << "subnet2 = " << subnet2->to_latex() << "\n"
             << "left1   = " << left1->to_latex() << "\n"
             << "left2   = " << left2->to_latex() << "\n";
#endif

  auto combine_expr = [](const ExprPtr& exprA, const ExprPtr& exprB) {
    if (exprA->empty()) return exprB->clone();
    if (exprB->empty()) return exprA->clone();
    auto combined = std::make_shared<Product>();
    combined->append(exprA->clone());
    combined->append(exprB->clone());
    return std::dynamic_pointer_cast<Expr>(combined);
  };

  auto factorForm1 = combine_expr(subnet1, left1);
  auto factorForm2 = combine_expr(subnet2, left2);

  return std::tuple(factorForm1, factorForm2);
}

// Get a tensor network of an expression: expr
// and a symmetry operator: symmetry_op  eg A, P
auto expr_to_tnet = [](const auto& expr, const auto& symmetry_op) {
  auto prod = std::make_shared<Product>();
  prod->append(1, symmetry_op->clone());
  prod->append(1, expr->at(0)->clone());
  auto tnet = TensorNetwork(*prod);
  return tnet;
};

// shallow copy tensors from network into a product form
// exclude tensors with labels ignoreByLabel
auto tnet_to_expr = [](const auto& tnet, std::wstring_view ignoreByLabel) {
  auto prod = std::make_shared<Product>();
  for (auto& expr : tnet.tensors()) {
    auto tnsr = std::dynamic_pointer_cast<Expr>(expr);
    if (tnsr->template as<Tensor>().label() != ignoreByLabel)
      prod->append(tnsr);
  }
  return prod;
};

ExprPtr fuse_pair(const ExprPtr& expr1, const ExprPtr& expr2,
                  const ExprPtr& symop) {
  // let's check if the input terms are fusable at all
  if (!(expr1->is<Product>() && expr2->is<Product>() &&
        expr1->size() == expr2->size() &&
        expr1->size() == 2 &&  // 2 since (AB..)(...) is expected
        expr1->at(0)->is<Product>() && expr2->at(0)->is<Product>())) {
    return nullptr;
  }

  auto commonTnet1 = expr_to_tnet(expr1, symop);
  auto commonTnet2 = expr_to_tnet(expr2, symop);

  auto namedIdx = TensorNetwork::named_indices_t{};

  auto phase1 = commonTnet1.canonicalize(
      TensorCanonicalizer::cardinal_tensor_labels(), false, &namedIdx);
  auto phase2 = commonTnet2.canonicalize(
      TensorCanonicalizer::cardinal_tensor_labels(), false, &namedIdx);

  if (!phase1) phase1 = std::make_shared<Constant>(Constant{1});

  if (!phase2) phase2 = std::make_shared<Constant>(Constant{1});

  auto fusedFact1 = tnet_to_expr(commonTnet1, symop->as<Tensor>().label());
  auto fusedFact2 = tnet_to_expr(commonTnet2, symop->as<Tensor>().label());

  if (evaluate::EvalTree(fusedFact1).hash_value() !=
      evaluate::EvalTree(fusedFact2).hash_value()) {
#ifdef DEBUG_PRINT
    std::wcout << "fused common factors"
                  R"(\\)"
                  "\n$"
               << fusedFact1->to_latex() << "$ \\quad $"
               << fusedFact2->to_latex() << "$\n";
#endif
    return nullptr;  // not fusable pair
  }

  auto remainder1 = expr1->at(1)->clone();
  auto remainder2 = expr2->at(1)->clone();

  // replace the indices in the remainders
  for (auto& expr : *remainder1) {
    auto& tnsr = expr->as<Tensor>();
    tnsr.transform_indices(commonTnet1.idxrepl());
  }
  for (auto& expr : *remainder2) {
    auto& tnsr = expr->as<Tensor>();
    tnsr.transform_indices(commonTnet2.idxrepl());
  }

  auto commonFact = std::make_shared<Product>();
  for (auto& expr : *fusedFact1)
    if (expr != symop) commonFact->append(1, expr);

  auto remainderSum = std::make_shared<Sum>();
  remainderSum->append(remainder1);
  remainderSum->append(remainder2);

  auto result = std::make_shared<Product>();

  result->append(commonFact);
  result->append(remainderSum);

  result->scale(phase1->as<Constant>().value());

  auto scalProd = phase1->as<Constant>().value() / result->scalar();
  auto scalSum = phase2->as<Constant>().value() / result->scalar();

  return result;
}

ExprPtr make_intermediate(const ExprPtr& expr1, const ExprPtr& expr2,
                          const ExprPtr& symOp, evaluate::Operation arithOp) {
  ExprPtr imed{nullptr};

  if (arithOp == evaluate::Operation::SUM) {
    auto imed = std::make_shared<Sum>();
    imed->append(expr1);
    imed->append(expr2);
  } else {
    auto imed = std::make_shared<Product>();
    imed->append(expr1);
    imed->append(expr2);
  }

  const auto imed_hash = evaluate::EvalTree(imed).hash_value();

  // auto tn = TensorNetwork(imed);

  auto& tnsr1 = expr1->as<Tensor>();
  auto& tnsr2 = expr2->as<Tensor>();

  container::set<Index> idx1, idx2;

  idx1.reserve(tnsr1.rank());
  idx2.reserve(tnsr2.rank());

  idx1.insert(tnsr1.const_braket().begin(), tnsr1.const_braket().end());
  idx2.insert(tnsr2.const_braket().begin(), tnsr2.const_braket().end());

  // canonicalization
  auto tn = expr_to_tnet(imed, symOp);
  auto namedIdx = TensorNetwork::named_indices_t{};
  auto phase = tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(),
                               false, &namedIdx);
  if (!phase) phase = std::make_shared<Constant>(Constant{1});

  auto canon_prod = tnet_to_expr(tn, symOp->as<Tensor>().label());

  // figure the product's bra and ket
  container::svector<Index> prod_bra{}, prod_ket{};

  // if 'item' exists in iterable, remove it
  // if 'item' doesn't exist in iterable, append it
  auto collect_unique = [](auto& iterable, const auto& items) {
    for (const auto item : items) {
      auto it = iterable.begin();
      for (; it != iterable.end(); ++it) {
        if (*it == item) iterable.erase(it);
        break;
      }
      if (it == iterable.end()) iterable.push_back(item);
    }
  };

  collect_unique(prod_bra, tnsr1.bra());
  collect_unique(prod_bra, tnsr2.bra());

  collect_unique(prod_ket, tnsr1.ket());
  collect_unique(prod_ket, tnsr2.ket());

  auto imed_tensor_ptr = std::make_shared<Tensor>(
      Tensor{L"I", prod_bra, prod_ket, Symmetry::nonsymm,
             BraKetSymmetry::nonsymm, ParticleSymmetry::nonsymm});

  // update hash value
  //
  auto& imed_tensor = imed_tensor_ptr->as<Tensor>();
  // imed_tensor->hash_value([&imed_hash](){return imed_hash;});

  return imed_tensor_ptr;
}

}  // namespace sequant::factorize
