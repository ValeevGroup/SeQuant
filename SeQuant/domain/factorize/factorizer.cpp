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
        container.size(), container::svector<bool>(container.size()));
    auto aa = 0;
    for (const auto& t1 : container) {
      auto bb = 0;
      for (const auto& t2 : container) {
        mat[aa][bb] = (aa == bb) ? false : are_connected(t1, t2);
        ++bb;
      }
      ++aa;
    }
    return mat;
  };

  const auto connect_mat = [&are_connected](const auto& container) {
    using hash_type = evaluate::HashType;
    // initialize the matrix
    container::svector<container::svector<hash_type>> mat(
        container.size(), container::svector<hash_type>(container.size()));
    auto aa = 0;
    for (const auto& t1 : container) {
      auto bb = 0;
      for (const auto& t2 : container) {
        auto hvalue = 0;
        if ((aa != bb) && are_connected(t1, t2))
          hvalue =
              evaluate::EvalTree(std::make_shared<Product>(Product{t1, t2}))
                  .hash_value();
        mat[aa][bb] = hvalue;
        ++bb;
      }
      ++aa;
    }
    return mat;
  };

  auto print_mat = [](const auto& mat, const auto& header,
                          const auto& sidebar) {
    std::wcout << "    ";
    for (const auto& head : header) std::wcout << head->to_latex() << " ";
    std::wcout << "\n";
    auto ii = 0;
    for (const auto& side : sidebar) {
      std::wcout << side->to_latex() << " ";
      for (const auto& entry : mat[ii])
        std::wcout << std::boolalpha << " " << entry << " ";
      std::wcout << "\n";
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
