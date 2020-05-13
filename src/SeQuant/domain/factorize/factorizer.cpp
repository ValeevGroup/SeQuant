#include "factorizer.hpp"

#include <SeQuant/core/container.hpp>

#include <tuple>

namespace sequant::factorize {

std::tuple<container::svector<size_t>, container::svector<size_t>>
largest_common_subnet(const ExprPtr& exprA, const ExprPtr& exprB) {
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
    return std::make_tuple(common_t1, common_t2);
  };  // lambda common_tensors

  auto [common_t1, common_t2] = common_tensors(*exprA, *exprB);

  // canonicalize the common tensors
  auto tn_t1 = TensorNetwork(common_t1);
  auto tn_t2 = TensorNetwork(common_t2);
  tn_t1.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);
  tn_t2.canonicalize(TensorCanonicalizer::cardinal_tensor_labels(), false);

  // collect common tensors after canonicalization
  auto tensors = [](const auto& container) {
    container::svector<ExprPtr> tensors;
    for (const auto& expr : container.tensors())
      tensors.push_back(std::dynamic_pointer_cast<Tensor>(expr));
    return tensors;
  };

  auto tensorsA = tensors(tn_t1);
  auto tensorsB = tensors(tn_t2);

  //
  // Now creating an adjacency matrix of common tensors
  // based on their connectivity.
  //
  auto are_connected = [](const ExprPtr& t1, const ExprPtr& t2) {
    auto tnsr1 = t1->as<Tensor>();
    auto tnsr2 = t2->as<Tensor>();
    // iterate through the bra and the ket labels of tnsr1 and tnsr2
    // if any index label is common, they are connected.
    for (const auto& idx1 : tnsr1.const_braket())
      for (const auto& idx2 : tnsr2.const_braket()) {
        if (idx1.label() == idx2.label()) return true;
      }
    return false;
  };
  //
  const auto adj_mat = [&are_connected](const auto& container) {
    // initialize the matrix
    container::svector<container::svector<bool>> mat(
        container.size(), container::svector<bool>(container.size(), false));
    for (auto ii = 0; ii < container.size() - 1; ++ii)
      for (auto jj = ii + 1; jj < container.size(); ++jj) {
        mat[ii][jj] = are_connected(container.at(ii), container.at(jj));
        mat[jj][ii] = mat[ii][jj];
      }

    return mat;
  };
  // creating connectivity matrix
  const auto connect_mat = [&are_connected](const auto& container) {
    using hash_type = evaluate::HashType;
    // initialize the matrix
    container::svector<container::svector<hash_type>> mat(
        container.size(), container::svector<hash_type>(container.size(), 0));
    for (auto ii = 0; ii < container.size() - 1; ++ii)
      for (auto jj = ii + 1; jj < container.size(); ++jj) {
        auto prod = std::make_shared<Product>(
            Product{container.at(ii), container.at(jj)});
        mat[ii][jj] = evaluate::EvalTree(prod).hash_value();
        mat[jj][ii] = mat[ii][jj];
      }

    return mat;
  };

  // get positions in a container where different kind of tensors begin
  // container: {f_ov, f_ov, t_oo, t_oovv, g_oovv, g_oovv}
  // output:    {(0, 2), (2, 3), (3, 4), (4, 6)}
  const auto parts_indices = [](const auto& container){
    using pos_type = size_t;
    using evaluate::EvalTree;// for hashing tensors by their types

    container::svector<pos_type> indices;
    indices.push_back(0);
    for (auto ii = 0; ii < container.size(); ) {
      auto lead_hvalue = EvalTree(container.at(ii)).hash_value();
      auto jj = ii + 1;
      for (; jj < container.size(); ++jj) {
        auto trail_hvalue = EvalTree(container.at(jj)).hash_value();
        if (lead_hvalue != trail_hvalue) break;
      }
      ii = jj;
      indices.push_back(ii);
    }
    container::svector<std::tuple<pos_type, pos_type>>  result;
    for (auto ii =0; ii < indices.size()-1; ++ii)
      result.push_back(std::make_tuple(indices.at(ii), indices.at(ii+1)));
    return result;
  };

  auto print_vec = [](const auto& container){
    for (const auto& idx: container)
      std::wcout << std::boolalpha << "  " << idx;
    std::wcout << std::endl;
  };

  auto print_mat = [&print_vec](const auto& mat, const auto& header,
                      const auto& sidebar) {
    std::wcout << "    ";
    for (const auto& head : header) std::wcout << head->to_latex() << " ";
    std::wcout << "\n";
    auto ii = 0;
    for (const auto& side : sidebar) {
      std::wcout << side->to_latex() << " ";
      print_vec(mat[ii]);
      ++ii;
    }
  };

  print_mat(adj_mat(tensorsA), tensorsA, tensorsA);
  print_mat(adj_mat(tensorsB), tensorsB, tensorsB);
  print_mat(connect_mat(tensorsA), tensorsA, tensorsA);
  print_mat(connect_mat(tensorsB), tensorsB, tensorsB);

  return std::make_tuple(container::svector<size_t>(),
                         container::svector<size_t>());
}

}  // namespace sequant::factorize
