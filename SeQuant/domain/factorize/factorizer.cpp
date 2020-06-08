#include "factorizer.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_network.hpp>
// IndexSpace type based hashing of tensors for ColorMatrix
#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <ios>
#include <tuple>

namespace sequant::factorize {

using pos_type = AdjacencyMatrix::pos_type;
using color_mat_type = AdjacencyMatrix::color_mat_type;

AdjacencyMatrix::AdjacencyMatrix(const ExprPtr& expr)
    : colorMatrix_(
          expr->size(),
          color_mat_type::value_type(
              expr->size(), color_mat_type::value_type::value_type{})) {
  for (auto ii = 0; ii < expr->size(); ++ii)
    for (auto jj = ii + 1; jj < expr->size(); ++jj) {
      // set color data
      if (are_connected(expr->at(ii), expr->at(jj))) {
        colorMatrix_[ii][jj] = colorMatrix_[jj][ii] =
            evaluate::EvalTree(
                std::make_shared<Product>(Product{expr->at(ii), expr->at(jj)}))
                .hash_value();

        /////
        /* auto t1 = std::make_shared<Tensor>(expr->at(ii)->as<Tensor>()); */
        /* auto t2 = std::make_shared<Tensor>(expr->at(jj)->as<Tensor>()); */

        /* auto tn = TensorNetwork(*std::make_shared<Product>(Product{t1, t2})); */

        /* tn.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false); */

        /* auto tensor = [&tn](auto idx) { */
        /*   return std::dynamic_pointer_cast<Tensor>( */
        /*       *(tn.tensors().begin() + idx)); */
        /* }; */

        /* auto vertex0 = tensor(0); */
        /* auto vertex1 = tensor(1); */

        /* colorMatrix_[ii][jj] = colorMatrix_[jj][ii] = */
        /*     evaluate::EvalTree( */
        /*         std::make_shared<Product>(Product{vertex0, vertex1})) */
        /*         .hash_value(); */
        /////
      }
    }
}

AdjacencyMatrix::AdjacencyMatrix(const container::svector<ExprPtr>& tensors)
    : AdjacencyMatrix(
          std::make_shared<Product>(1, tensors.begin(), tensors.end())) {}

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

std::tuple<ExprPtr, ExprPtr> factorize_pair(const ExprPtr& exprA,
                                            const ExprPtr& exprB) {
  // lambda to get the common tensors in container_t1 and container_t2
  auto common_tensors = [](const auto& container_t1, const auto& container_t2) {
    container::set<ExprPtr> common_t1, common_t2;

    for (const auto& t1 : container_t1)
      for (const auto& t2 : container_t2) {
        //
        // NOTE: As an example: t_{i j}^{a b} has the same
        // hash value as t_{a b}^{i j}. To hash such expressions
        // differently, use EvalTree(expr, false).
        //
        if (evaluate::EvalTree(t1).hash_value() ==
            evaluate::EvalTree(t2).hash_value()) {
          common_t1.insert(t1);
          common_t2.insert(t2);
        }
      }

    // convert set to vector
    auto exprptr_vec = [](const auto& container) {
      container::svector<ExprPtr> result;
      result.reserve(container.size());
      for (auto&& xpr : container) result.push_back(xpr);
      return result;
    };

    return std::make_tuple(exprptr_vec(common_t1), exprptr_vec(common_t2));
  };  // lambda common_tensors

  auto [tensorsA, tensorsB] = common_tensors(*exprA, *exprB);

  // get positions in a container where different 'kind' of tensors begin
  // container: {f_ov, f_ov, t_oo, t_oovv, g_oovv, g_oovv}
  // output:    {(0, 2), (2, 3), (3, 4), (4, 6)}
  const auto parts_indices = [](const auto& container) {
    using evaluate::EvalTree;  // for hashing tensors by their types

    container::svector<pos_type> indices;
    indices.push_back(0);
    for (auto ii = 0; ii < container.size();) {
      auto lead_hvalue = EvalTree(container.at(ii)).hash_value();
      auto jj = ii + 1;
      for (; jj < container.size(); ++jj) {
        auto trail_hvalue = EvalTree(container.at(jj)).hash_value();
        if (lead_hvalue != trail_hvalue) break;
      }
      ii = jj;
      indices.push_back(ii);
    }
    container::svector<std::tuple<pos_type, pos_type>> result;
    for (auto ii = 0; ii < indices.size() - 1; ++ii)
      result.push_back(std::make_tuple(indices.at(ii), indices.at(ii + 1)));
    return result;
  };

  // get pairs of tensors from each of the two containers, container1 and
  // container2 where a pair from 1 has the same color of connectivity as the
  // corresponding pair from 2.
  //
  // returned value is a tuple of two vectors each vector has tuple(s) of two
  // position indices
  // eg.
  //     ({(0,4), (1,3)}, {(1,7), (3, 8)})
  //
  //     tensors at positions 0 and 4 from container1 are equivalently
  //     contracted as tensors at positions 1 and 7 from container2
  //
  //     tensors at positions 1 and 3 from container1 are equivalently
  //     contracted as tensors at positions 3 and 8 from container2
  //
  // container1 and container2 maybe different in size, however, they must have
  // same number of kinds of tensors
  // eg.
  // {t_ov, t_ov, t_ov, g_oovv} and {t_ov, g_oovv} are good
  // {t_ov, t_ov, t_ov, g_oovv, f_ov} and {t_ov, g_oovv} are not good
  //
  auto common_pairs = [&parts_indices](const auto& container1,
                                       const auto& container2) {
    // given @c pos position of an element in a container and
    // given @c parts partition vector of the container
    // return the index in the @c parts container in which the @c pos falls.
    // eg.
    // parts: {(0, 2), (2, 5), (5, 6)}
    // pos: 3
    // return: 1
    // since 3 belongs to (2, 5) partition whose index in the
    // partition vector is 1
    const auto part_idx = [](pos_type pos, const auto& parts) {
      for (size_t ii = 0; ii < parts.size(); ++ii) {
        if ((std::get<0>(parts.at(ii)) <= pos) &&
            (std::get<1>(parts.at(ii)) > pos))
          return ii;
      }
      return parts.size();
    };

    auto parts1 = parts_indices(container1);
    auto parts2 = parts_indices(container2);
    // assert(parts1.size() == parts2.size());

    auto adjMat1 = AdjacencyMatrix(container1);
    auto adjMat2 = AdjacencyMatrix(container2);

    container::svector<std::tuple<pos_type, pos_type>> commonPair1, commonPair2;
    for (auto p1 = 0; p1 < adjMat1.num_verts(); ++p1) {
      for (auto pp1 = p1 + 1; pp1 < adjMat1.num_verts(); ++pp1) {
        const auto [low1, up1] = parts2.at(part_idx(p1, parts1));
        const auto [low2, up2] = parts2.at(part_idx(pp1, parts1));
        for (auto p2 = low1; p2 < up1; ++p2)
          for (auto pp2 = low2; pp2 < up2; ++pp2) {
            if (auto color = adjMat1.color(p1, pp1);
                (color != color_mat_type::value_type::value_type{}) &&
                (color == adjMat2.color(p2, pp2))) {
              commonPair1.push_back(std::make_tuple(p1, pp1));
              commonPair2.push_back(std::make_tuple(p2, pp2));
            }
          }
      }
    }

    return std::make_tuple(commonPair1, commonPair2);
  };

  auto [commonPairsA, commonPairsB] = common_pairs(tensorsA, tensorsB);
  assert(commonPairsA.size() == commonPairsB.size());

  // keep track of indices that make a walk/common subnet
  container::svector<container::set<pos_type>> processedNetsA, processedNetsB;

  // keep track of tuple's positions that have been consumed in a walk
  container::svector<bool> processedTuples(commonPairsA.size(), false);

  for (auto ii = 0; ii < commonPairsA.size(); ++ii) {
    if (processedTuples[ii]) continue;
    decltype(processedNetsA)::value_type runningNetA, runningNetB;

    runningNetA.insert(std::get<0>(commonPairsA.at(ii)));
    runningNetA.insert(std::get<1>(commonPairsA.at(ii)));

    runningNetB.insert(std::get<0>(commonPairsB.at(ii)));
    runningNetB.insert(std::get<1>(commonPairsB.at(ii)));

    processedTuples[ii] = true;

    for (auto jj = ii + 1; jj < commonPairsA.size(); ++jj) {
      if (processedTuples[jj]) continue;
      auto [idxA1, idxA2] = commonPairsA.at(jj);
      auto [idxB1, idxB2] = commonPairsB.at(jj);
      if (runningNetA.contains(idxA1) || runningNetA.contains(idxA2)) {
        assert(runningNetB.contains(idxB1) || runningNetB.contains(idxB2));
        runningNetA.insert(idxA1);
        runningNetA.insert(idxA2);

        runningNetB.insert(idxB1);
        runningNetB.insert(idxB2);

        processedTuples[jj] = true;
      }
    }
    processedNetsA.push_back(runningNetA);
    processedNetsB.push_back(runningNetB);
  }

  /* std::wcout << "positions in netA\n"; */
  /* for (const auto& pos : processedNetsA) { */
  /*   std::wcout << "("; */
  /*   for (auto p : pos) std::wcout << p << "  "; */

  /*   std::wcout << ")\n"; */
  /* } */

  /* std::wcout << "positions in netB\n"; */
  /* for (const auto& pos : processedNetsB) { */
  /*   std::wcout << "("; */
  /*   for (auto p : pos) std::wcout << p << "  "; */

  /*   std::wcout << ")\n"; */
  /* } */

  auto form_product = [](const auto& nets, const auto& tensors) {
    auto result = std::make_shared<Product>();
    for (const auto& n : nets) {
      auto prod = std::make_shared<Product>();
      for (auto nn : n) {
        prod->append(tensors.at(nn));
      }
      result->append(prod);
    }
    return result;
  };

  auto subnetA = form_product(processedNetsA, tensorsA);
  auto subnetB = form_product(processedNetsB, tensorsB);
  /* std::wcout << "subnetA = " << subnetA->to_latex() << "\n"; */
  /* std::wcout << "subnetB = " << subnetB->to_latex() << "\n"; */

  decltype(tensorsA) unfactoredA, unfactoredB;
  for (const auto& tnsr : *exprA)
    if (!tensor_exists(subnetA, tnsr)) unfactoredA.push_back(tnsr);
  for (const auto& tnsr : *exprB)
    if (!tensor_exists(subnetB, tnsr)) unfactoredB.push_back(tnsr);

  subnetA->append(
      1, std::make_shared<Product>(1, unfactoredA.begin(), unfactoredA.end()));
  subnetA->scale(exprA->as<Product>().scalar());

  subnetB->append(
      1, std::make_shared<Product>(1, unfactoredB.begin(), unfactoredB.end()));
  subnetB->scale(exprB->as<Product>().scalar());

  return std::tuple(subnetA, subnetB);
}

}  // namespace sequant::factorize
