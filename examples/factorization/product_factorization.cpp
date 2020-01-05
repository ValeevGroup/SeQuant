//
// Created by Bimal Gaudel on 12/22/19.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>
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
  using ispace_pair = std::pair<sequant::IndexSpace::Type, size_t>;
  using ispace_map = sequant::container::map<ispace_pair::first_type, ispace_pair::second_type>;

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
  prod1.append(2.0, t1);
  prod1.append(1.5, t2);
  prod1.append(1.5, t3);
  prod1.append(1.5, t4);

  // std::wcout << prod1.to_latex();
  std::wcout << "\nno. of factors = " << prod1.factors().size() << "\n";
  std::wcout << "\nSetting up a map with nocc = 10 and nvirt = 20..\n";

  auto counter_map = std::make_shared<ispace_map>(ispace_map{});

  counter_map->insert(ispace_pair{sequant::IndexSpace::active_occupied, 10});
  counter_map->insert(ispace_pair{sequant::IndexSpace::active_unoccupied, 20});

  std::wcout << "\nProduct about to be factorized..\n" << prod1.to_latex()
             << "\nAfter factorization the product looks like..\n";
  auto result = factorize_product(prod1, counter_map);
  std::wcout << result->to_latex() << "\n";
  return 0;
}
