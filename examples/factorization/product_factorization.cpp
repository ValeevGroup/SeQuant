//
// Created by Bimal Gaudel on 12/22/19.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>
#include <SeQuant/domain/factorize/path_tree.hpp>
#include "../sequant_setup.hpp"

#include <memory>
#include <limits>
#include <string>
#include <iostream>

int main() {
  // global sequant setup...
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  Logger::get_instance().wick_stats = false;
  /* auto cc_r = cceqvec{ 2, 2 }(true, true, true, true); */

  using sequant::Tensor;
  using sequant::factorize::factorize_product;
  using sequant::factorize::detail::ContractionCostResult;
  using sequant::factorize::detail::PathCostResult;
  using sequant::factorize::detail::ContractionCostCounter;
  using sequant::factorize::detail::optimal_path;
  using sequant::factorize::detail::PathTree;

  /* using PathPtr = std::shared_ptr<PathTree>; */
  /* using vec_path_ptr = sequant::container::svector<PathPtr>; */

  std::wcout << "Creating a product..\n";
  std::shared_ptr<Tensor> t1, t2, t3, t4;
  t1 = std::make_shared<Tensor>(
      Tensor{L"t", {L"i_1", L"i_2"}, {L"a_1", L"a_2"}});
  t2 = std::make_shared<Tensor>(
      Tensor{L"t", {L"i_3", L"i_4"}, {L"a_1", L"a_3"}});
  t3 = std::make_shared<Tensor>(
      Tensor{L"t", {L"a_2", L"a_3"}, {L"a_4", L"a_5"}});
  t4 = std::make_shared<Tensor>(
      Tensor{L"t", {L"i_1", L"i_2"}, {L"i_3", L"a_6"}});

  auto prod1 = sequant::Product{};
  prod1.append(1, t1);
  prod1.append(1, t2);
  prod1.append(1, t3);
  prod1.append(1, t4);

  std::wcout << prod1.to_latex();
  std::wcout << "\nno. of factors = " << prod1.factors().size() << "\n";
  std::wcout << "\nSetting up a flops counter object with nocc = 10 and nvirt "
                "= 20..\n";
  auto counter = ContractionCostCounter{10, 20};

  auto initiate_path = [&](sequant::Product prod) {
    auto result = sequant::container::svector<std::shared_ptr<PathTree>>(
        prod.factors().size());
    for (ulong i = 0; i < prod.factors().size(); ++i)
      result[i] = std::make_shared<PathTree>(PathTree{i});
    return result;
  };

  auto initial_path = initiate_path(prod1);

  auto running_cost = std::make_shared<PathCostResult>();

  optimal_path(initial_path, prod1, counter, running_cost);

   std::wcout << "\nThe optimal path must be "
             << running_cost->path->print_tree()
             << " " << running_cost->flops << "\n";

  // testing path_to_product
  auto tree0 = std::make_shared<PathTree>(0);
  auto tree1 = std::make_shared<PathTree>(1);
  auto tree2 = std::make_shared<PathTree>(2);
  auto tree3 = std::make_shared<PathTree>(3);
  tree0->add_child(tree3);
  tree2->add_child(tree1);
  tree0->add_child(tree2);
  std::wcout << "Printing the tree..\n" << tree0->print_tree() << "\n";
  auto new_prod = sequant::factorize::detail::path_to_product(tree0, prod1);
  std::wcout << "New product in latex..\n" << new_prod->to_latex() << "\n";

  std::wcout << "\nProduct about to be factorized..\n" << prod1.to_latex()
             << "\nAfter factorization the product looks like..\n";
  auto result = factorize_product(prod1, counter);
  std::wcout << result->to_latex() << "\n";
  return 0;
}
