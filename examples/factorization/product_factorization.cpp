//
// Created by Bimal Gaudel on 12/22/19.
//

#include <btas/btas.h>
#include <btas/tensorview.h>
#include <btas/tensor_func.h>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>
#include "../contract/interpret/contract.hpp"

#include "../sequant_setup.hpp"

#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include <functional>  // for std::bind
#include <random>      // for std::mt19937
#include <chrono>      // for seeding

// factorize an Expr
// Note: antisymmetrization Expr's will be gone from the returned result
sequant::ExprPtr factorize_expr(const ExprPtr& expr_ptr,
    const std::shared_ptr<sequant::container::map<sequant::IndexSpace::Type, size_t>>& ispace_map,
    bool factorize=true){
  using sequant::ExprPtr;
  using sequant::SumPtr;
  using sequant::ProductPtr;
  using sequant::Sum;
  using sequant::Product;
  using sequant::Tensor;

  if (expr_ptr->is<Product>()){
    // form a new product by omitting antisymmetrization tensors
    // and return the result of the sequant::factorize::product_factorize on it

    ProductPtr new_prod(new Product{});
    auto& prod = expr_ptr->as<Product>();
    for (auto &&fac: prod) {
      // CATCH: no further checking if the factors are non-Tensor type
      auto& tnsr = fac->as<Tensor>();
      // skip antisym. tensors
      if (tnsr.label() != L"A") new_prod->append(1, fac);
    }
    new_prod->scale(prod.scalar());

    if (factorize)
      return sequant::factorize::factorize_product(new_prod->as<Product>(), ispace_map);
    else return new_prod;

  } else if (expr_ptr->is<Sum>()){

    SumPtr new_sum(new Sum{});
    auto& sum = expr_ptr->as<Sum>();
    for (auto &&sumand: sum)
      new_sum->append(factorize_expr(sumand, ispace_map, factorize));
    return new_sum;
  } else {
    return expr_ptr;
  }
}

int main() {
  // global sequant setup...
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wostream::sync_with_stdio(true);
  std::wostream::sync_with_stdio(true);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wostream::sync_with_stdio(true);
  std::wostream::sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  Logger::get_instance().wick_stats = false;
  // CCSD equations
  auto cc_r = cceqvec{2, 2}(true, true, true, true);

  // factorization and evaluation of the CC equations
  using sequant::factorize::factorize_product;
  using ispace_pair = std::pair<sequant::IndexSpace::Type, size_t>;
  using ispace_map = sequant::container::map<ispace_pair::first_type,
                                             ispace_pair::second_type>;

  size_t nocc = 5, nvirt = 20;
  std::wcout << "\nSetting up a map with nocc = " << nocc << " and nvirt = " << nvirt << "..\n";
  auto counter_map = std::make_shared<ispace_map>(ispace_map{});
  counter_map->insert(ispace_pair{sequant::IndexSpace::active_occupied, nocc});
  counter_map->insert(
      ispace_pair{sequant::IndexSpace::active_unoccupied, nvirt});

  // initializing tensors present in the CCSD equations
  btas::Tensor<double> Fock_oo (nocc, nocc),
                       Fock_ov (nocc, nvirt),
                       Fock_vv (nvirt, nvirt),
                       G_oooo  (nocc, nocc, nocc, nocc),
                       G_vvvv  (nvirt, nvirt, nvirt, nvirt),
                       G_ovvv  (nocc, nvirt, nvirt, nvirt),
                       G_ooov  (nocc, nocc, nocc, nvirt),
                       G_oovv  (nocc, nocc, nvirt, nvirt),
                       G_ovov  (nocc, nvirt, nocc, nvirt),
                       T_ov    (nocc, nvirt),
                       T_oovv  (nocc, nocc, nvirt, nvirt);
  // fill random data to tensors
  auto randfill_tensor = [&](btas::Tensor<double>& tnsr){
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rgen(seed);
    auto dist = std::uniform_real_distribution<double>{-1.0, 1.0};
    tnsr.generate(std::bind(dist, rgen));
  };

   randfill_tensor(Fock_oo);
   randfill_tensor(Fock_ov);
   randfill_tensor(Fock_vv);
   randfill_tensor(G_oooo );
   randfill_tensor(G_vvvv );
   randfill_tensor(G_ovvv );
   randfill_tensor(G_ooov );
   randfill_tensor(G_oovv );
   randfill_tensor(G_ovov );
   randfill_tensor(T_ov   );
   randfill_tensor(T_oovv );

   /* std::wcout << L"\nnorm(Fock_oo) = " << std::sqrt(btas::dot(Fock_oo,Fock_oo)) */
   /*            << L"\nnorm(Fock_ov) = " << std::sqrt(btas::dot(Fock_ov,Fock_ov)) */
   /*            << L"\nnorm(Fock_vv) = " << std::sqrt(btas::dot(Fock_vv,Fock_vv)) */
   /*            << L"\nnorm(G_oooo ) = " << std::sqrt(btas::dot(G_oooo,G_oooo)) */
   /*            << L"\nnorm(G_vvvv ) = " << std::sqrt(btas::dot(G_vvvv,G_vvvv)) */
   /*            << L"\nnorm(G_ovvv ) = " << std::sqrt(btas::dot(G_ovvv,G_ovvv)) */
   /*            << L"\nnorm(G_ooov ) = " << std::sqrt(btas::dot(G_ooov,G_ooov)) */
   /*            << L"\nnorm(G_oovv ) = " << std::sqrt(btas::dot(G_oovv,G_oovv)) */
   /*            << L"\nnorm(G_ovov ) = " << std::sqrt(btas::dot(G_ovov,G_ovov)) */
   /*            << L"\nnorm(T_ov   ) = " << std::sqrt(btas::dot(T_ov,T_ov)) */
   /*            << L"\nnorm(T_oovv ) = " << std::sqrt(btas::dot(T_oovv,T_oovv)); */

   // for the evaluation of the expressions
   // a map is needed
   using str_to_tensor_pair = std::pair<std::wstring, btas::Tensor<double> const *>;
   using str_to_tensor_map = std::map<str_to_tensor_pair::first_type,
                                       str_to_tensor_pair::second_type>;

   str_to_tensor_map btensor_map;
   btensor_map.insert(str_to_tensor_pair(L"f_oo",   &Fock_oo));
   btensor_map.insert(str_to_tensor_pair(L"f_ov",   &Fock_ov));
   btensor_map.insert(str_to_tensor_pair(L"f_vv",   &Fock_vv));
   btensor_map.insert(str_to_tensor_pair(L"g_oooo", &G_oooo));
   btensor_map.insert(str_to_tensor_pair(L"g_vvvv", &G_vvvv));
   btensor_map.insert(str_to_tensor_pair(L"g_ovvv", &G_ovvv));
   btensor_map.insert(str_to_tensor_pair(L"g_ooov", &G_ooov));
   btensor_map.insert(str_to_tensor_pair(L"g_oovv", &G_oovv));
   btensor_map.insert(str_to_tensor_pair(L"g_ovov", &G_ovov));
   btensor_map.insert(str_to_tensor_pair(L"t_ov",   &T_ov));
   btensor_map.insert(str_to_tensor_pair(L"t_oovv", &T_oovv));

   // factorization and evaluation
   auto& expr_to_factorize = cc_r[2];
   auto unfactorized_expr  = factorize_expr(expr_to_factorize, counter_map, false);
   auto factorized_expr    = factorize_expr(expr_to_factorize, counter_map, true);

   std::wcout << "\nUnfactorized..\n"
              << unfactorized_expr->to_latex()
              << "\n"
              << "\nFactorized..\n"
              << factorized_expr->to_latex()
              << "\n";
  // to confirm that we are not working
  // on the same Expr
  if (*unfactorized_expr != *factorized_expr)
       std::wcout << "\nunfactorized and factorized Epr are not the same.. which is good:)\n";
  else std::wcout << "\nunfactorized and factorized Expr are the same.. time to debug:(\n";

  using std::chrono::duration_cast;
  using std::chrono::microseconds;

  auto tstart            = std::chrono::high_resolution_clock::now();
  auto unfactorized_eval = sequant::interpret::eval_equation(unfactorized_expr, btensor_map);
  auto tstop             = std::chrono::high_resolution_clock::now();
  auto duration          = duration_cast<microseconds>(tstop - tstart).count();
  std::wcout << "time(unfactorized_eval) = " << duration << " microseconds.\n";

  tstart                 = std::chrono::high_resolution_clock::now();
  auto factorized_eval   = sequant::interpret::eval_equation(factorized_expr, btensor_map);
  tstop                  = std::chrono::high_resolution_clock::now();
  duration               = duration_cast<microseconds>(tstop - tstart).count();
  std::wcout << "time(factorized_eval) = " << duration << " microseconds.\n";

  std::wcout << "\nnorm(unfac) = "
    << std::sqrt(btas::dot(unfactorized_eval.tensor(), unfactorized_eval.tensor()))
    << "\n";
  std::wcout << "\nnorm(fac) = "
    << std::sqrt(btas::dot(factorized_eval.tensor(), factorized_eval.tensor()))
    << "\n";
  return 0;
}