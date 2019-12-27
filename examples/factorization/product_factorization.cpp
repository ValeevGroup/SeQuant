//
// Created by Bimal Gaudel on 12/22/19.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>
#include <SeQuant/domain/factorize/path_tree.hpp>
#include "../sequant_setup.hpp"

#include <iostream>
#include <memory>
#include <string>

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
  using sequant::factorize::detail::CostCounter;
  using sequant::factorize::detail::path_scanner;
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

  std::wcout << "\nSetting up a flops counter object with nocc = 10 and nvirt "
                "= 20..\n";
  auto counter = CostCounter{10, 20};

  auto initiate_path = [&](sequant::Product prod) {
    auto result = sequant::container::svector<std::shared_ptr<PathTree>>(
        prod.factors().size());
    for (ulong i = 0; i < prod.factors().size(); ++i)
      result[i] = std::make_shared<PathTree>(PathTree{i});
    return result;
  };

  auto initial_path = initiate_path(prod1);

  path_scanner(initial_path, prod1, counter);

  return 0;
}
