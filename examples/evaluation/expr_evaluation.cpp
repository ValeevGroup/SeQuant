#include <SeQuant/domain/evaluate/eval_expr.hpp>
#include <SeQuant/domain/evaluate/eval_fwd.hpp>
#include <SeQuant/domain/evaluate/eval_tensor.hpp>
#include <SeQuant/domain/factorize/factorizer.hpp>
#include "../contract/interpret/contract.hpp"
#include "../sequant_setup.hpp"

#include <btas/btas.h>
#include <btas/tensorview.h>

#include <complex>
#include <iostream>

#include <chrono>      // for seeding
#include <functional>  // for std::bind
#include <random>      // for std::mt19937

void print_evtensor(const sequant::evaluate::EvTensorPtr&);
void print_leaf_hashes(const sequant::evaluate::EvTensorPtr&);

using hash_count_pair = std::pair<sequant::evaluate::hash_type, std::size_t>;
using hash_count_map =
    container::map<hash_count_pair::first_type, hash_count_pair::second_type>;

void count_hashes(const sequant::evaluate::EvTensorPtr& evt_ptr,
                  hash_count_map& counts) {
  if (evt_ptr->is_leaf()) return;

  counts[evt_ptr->get_hash_value()] += 1;

  count_hashes(evt_ptr->left_tensor(), counts);
  count_hashes(evt_ptr->right_tensor(), counts);
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

  using ispace_pair = std::pair<sequant::IndexSpace::Type, size_t>;
  using ispace_map = sequant::container::map<ispace_pair::first_type,
                                             ispace_pair::second_type>;
  std::size_t nocc = 10, nvirt = 16;
  std::wcout << "\nSetting up a map with nocc = " << nocc
             << " and nvirt = " << nvirt << "..\n";
  auto fac_map = std::make_shared<ispace_map>(ispace_map{});
  fac_map->insert(ispace_pair{sequant::IndexSpace::active_occupied, nocc});
  fac_map->insert(ispace_pair{sequant::IndexSpace::active_unoccupied, nvirt});

  // CC equations
  auto cc_r = cceqvec{2, 2}(true, true, true, true);
  auto expr_to_factorize = cc_r[2];
  /* auto expr_to_factorize = std::make_shared<sequant::Sum>(); */
  /* size_t counter = 0; */
  /* for (const auto& summand : cc_r[2]->as<Sum>()) { */
  /*   if (counter < 5) { */
  /*     expr_to_factorize->append(summand); */
  /*   } */
  /*   ++counter; */
  /* } */
  if (expr_to_factorize->is<Product>())
    std::wcout << "The expr to be factorized is a Prod.\n";
  else if (expr_to_factorize->is<Sum>())
    std::wcout << "The expr to be factorized is a Sum.\n";

  auto factorized_expr =
      sequant::factorize::factorize_expr(expr_to_factorize, fac_map, true);
  std::wcout << "Done factorization" << std::endl;
  // std::wcout << "factorize_expr->to_latex(): \n" <<
  // factorized_expr->to_latex() << "\n";

  using std::endl;
  using std::wcout;
  using std::chrono::duration_cast;
  using std::chrono::microseconds;
  using BTensor = btas::Tensor<double>;

  auto Fock_oo = std::make_shared<BTensor>(nocc, nocc);
  auto Fock_ov = std::make_shared<BTensor>(nocc, nvirt);
  auto Fock_vv = std::make_shared<BTensor>(nvirt, nvirt);
  auto G_oooo = std::make_shared<BTensor>(nocc, nocc, nocc, nocc);
  auto G_vvvv = std::make_shared<BTensor>(nvirt, nvirt, nvirt, nvirt);
  auto G_ovvv = std::make_shared<BTensor>(nocc, nvirt, nvirt, nvirt);
  auto G_ooov = std::make_shared<BTensor>(nocc, nocc, nocc, nvirt);
  auto G_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
  auto G_ovov = std::make_shared<BTensor>(nocc, nvirt, nocc, nvirt);
  auto T_ov = std::make_shared<BTensor>(nocc, nvirt);
  auto T_oovv = std::make_shared<BTensor>(nocc, nocc, nvirt, nvirt);
  /* auto T_ooovvv = */
  /*     std::make_shared<BTensor>(nocc, nocc, nocc, nvirt, nvirt, nvirt); */

  // fill random data to tensors
  auto randfill_tensor = [&](BTensor& tnsr) {
    auto seed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rgen(seed);
    auto dist = std::uniform_real_distribution<double>{-1.0, 1.0};
    tnsr.generate(std::bind(dist, rgen));
  };

  randfill_tensor(*Fock_oo);
  randfill_tensor(*Fock_ov);
  randfill_tensor(*Fock_vv);
  randfill_tensor(*G_oooo);
  randfill_tensor(*G_vvvv);
  randfill_tensor(*G_ovvv);
  randfill_tensor(*G_ooov);
  randfill_tensor(*G_oovv);
  randfill_tensor(*G_ovov);
  randfill_tensor(*T_ov);
  randfill_tensor(*T_oovv);
  /* randfill_tensor(T_ooovvv); */

  std::wstring label_f = L"f";
  std::wstring label_g = L"g";
  std::wstring label_t = L"t";

  using namespace sequant::evaluate;

  auto generate_hash_map_entry = [](std::wstring& label, std::string&& spaces,
                                    std::shared_ptr<BTensor>& tnsr_ptr) {
    auto space_vector = container::svector<IndexSpace::Type>{};
    for (char sp : spaces) {
      if (sp == 'o')
        space_vector.push_back(IndexSpace::active_occupied);
      else if (sp == 'v')
        space_vector.push_back(IndexSpace::active_unoccupied);
      else
        throw std::logic_error(
            "Use 'o' for active_occupied and 'v' for active_unoccupied "
            "IndexSpace!");
    }
    auto hval = sequant::evaluate::DataTensorSpecs(label, space_vector)
                    .get_hash_value();
    return hash_to_dtensor_map<BTensor>::value_type(hval, tnsr_ptr);
  };

  auto hash_map_ptr = std::make_shared<hash_to_dtensor_map<BTensor>>();

  hash_map_ptr->insert(generate_hash_map_entry(label_f, "oo", Fock_oo));
  hash_map_ptr->insert(generate_hash_map_entry(label_f, "ov", Fock_ov));
  hash_map_ptr->insert(generate_hash_map_entry(label_f, "vv", Fock_vv));
  hash_map_ptr->insert(generate_hash_map_entry(label_g, "oooo", G_oooo));
  hash_map_ptr->insert(generate_hash_map_entry(label_g, "vvvv", G_vvvv));
  hash_map_ptr->insert(generate_hash_map_entry(label_g, "ovvv", G_ovvv));
  hash_map_ptr->insert(generate_hash_map_entry(label_g, "oovv", G_oovv));
  hash_map_ptr->insert(generate_hash_map_entry(label_g, "ooov", G_ooov));
  hash_map_ptr->insert(generate_hash_map_entry(label_g, "ovov", G_ovov));
  hash_map_ptr->insert(generate_hash_map_entry(label_t, "ov", T_ov));
  hash_map_ptr->insert(generate_hash_map_entry(label_t, "oovv", T_oovv));
  // hash_map_ptr->insert(generate_hash_map_entry(label_t, "ooovvv", T_ooovvv));
  //
  // CHECKS OUT
  // for (const auto& item: *hash_map_ptr)
  //   std::wcout << item.first << " " << btas::dot(*item.second, *item.second)
  //   << "\n";
  // std::wcout << "\n\n";
  // for (const auto& item: context.leaf_map())
  //     std::wcout /* << item.first */ << " " << btas::dot(*item.second,
  //     *item.second) << "\n";
  // print_leaf_hashes(evt_ptr);

  std::wcout << "Initialized tensors for new evaluation method.." << std::endl;

  // create a EvalTensor from the factorized expression
  auto evt_ptr = std::make_shared<EvalTensor>(factorized_expr);
  // print_evtensor(evt_ptr);

  // creating a context to evaluate on
  auto hash_counts = container::map<hash_type, std::size_t>{};
  fill_hash_counts(evt_ptr, hash_counts);
  auto context = EvalContext<BTensor>(*hash_map_ptr, hash_counts);

  // for (const auto& item: context.imed_counts())
  //   std::wcout << "hash: " << item.first << "    counts: " << item.second <<
  //   "\n";

#ifdef SEQUANT_HAS_BTAS
  evt_ptr->fill_btas_indices();
  std::wcout << "Filled btas indices..entering evaluation phase.." << std::endl;
  auto tstart = std::chrono::high_resolution_clock::now();
  //
  auto result = eval_evtensor<BTensor>(evt_ptr, context);
  auto top_scalar = evt_ptr->get_scalar().real();
  if (top_scalar != 1.0) {
    btas::scal(top_scalar, result);
  }
  auto tstop = std::chrono::high_resolution_clock::now();
   std::wcout << "norm(result) = " << std::sqrt(btas::dot(result, result))
              << std::endl;
  auto duration = duration_cast<microseconds>(tstop - tstart).count();
  std::wcout << "time(new_method) = " << duration << " microseconds.\n";

  using str_to_tensor_pair =
      std::pair<std::wstring, std::shared_ptr<BTensor>>;
  using str_to_tensor_map =
      std::map<str_to_tensor_pair::first_type, str_to_tensor_pair::second_type>;

  // for older evaluations
  str_to_tensor_map btensor_map;
  btensor_map.insert(str_to_tensor_pair(L"f_oo", Fock_oo));
  btensor_map.insert(str_to_tensor_pair(L"f_ov", Fock_ov));
  btensor_map.insert(str_to_tensor_pair(L"f_vv", Fock_vv));
  btensor_map.insert(str_to_tensor_pair(L"g_oooo", G_oooo));
  btensor_map.insert(str_to_tensor_pair(L"g_vvvv", G_vvvv));
  btensor_map.insert(str_to_tensor_pair(L"g_ovvv", G_ovvv));
  btensor_map.insert(str_to_tensor_pair(L"g_ooov", G_ooov));
  btensor_map.insert(str_to_tensor_pair(L"g_oovv", G_oovv));
  btensor_map.insert(str_to_tensor_pair(L"g_ovov", G_ovov));
  btensor_map.insert(str_to_tensor_pair(L"t_ov", T_ov));
  btensor_map.insert(str_to_tensor_pair(L"t_oovv", T_oovv));

  // for (const auto& item: btensor_map)
  //  std::wcout << "norm = " << btas::dot(*item.second, *item.second) << "\n";

  std::wcout << "\nevaluation using the old method.." << std::endl;
  tstart = std::chrono::high_resolution_clock::now();
  auto result2 = interpret::eval_equation(factorized_expr, btensor_map);
  tstop = std::chrono::high_resolution_clock::now();

  std::wcout << "norm(result2) = "
             << std::sqrt(btas::dot(result2.tensor(), result2.tensor()))
             << "\n";
  duration = duration_cast<microseconds>(tstop - tstart).count();
  std::wcout << "time(old_method) = " << duration << " microseconds.\n";

#endif

  return 0;
}

void print_evtensor(const sequant::evaluate::EvTensorPtr& evt_ptr) {
  std::wcout << "Scalar: " << evt_ptr->get_scalar() << "    Indices: ";
  for (const auto& idx : evt_ptr->indices()) std::wcout << idx << " ";
  std::wcout << "   Hash: " << evt_ptr->get_hash_value() << " Operation:  ";
  if (evt_ptr->get_op() == sequant::evaluate::EvalTensor::Operation::Sum)
    std::wcout << " Sum";
  else if (evt_ptr->get_op() ==
           sequant::evaluate::EvalTensor::Operation::Product)
    std::wcout << " Product";
  else if (evt_ptr->get_op() == sequant::evaluate::EvalTensor::Operation::Eval)
    std::wcout << " Eval";
  else
    throw std::logic_error(" Invalid Operation!");
  std::wcout << std::endl;
#ifdef SEQUANT_HAS_BTAS
  std::wcout << "Btas indices: ";
  for (const auto& idx : evt_ptr->btas_indices()) std::wcout << idx << " ";
  std::wcout << "\n";
#endif

  if (evt_ptr->is_leaf()) return;

  print_evtensor(evt_ptr->left_tensor());
  print_evtensor(evt_ptr->right_tensor());
}

void print_leaf_hashes(const sequant::evaluate::EvTensorPtr& evt_ptr) {
  if (evt_ptr->is_leaf()) {
    std::wcout << "leaf hash = " << evt_ptr->get_hash_value() << "\n";
    return;
  }
  print_leaf_hashes(evt_ptr->left_tensor());
  print_leaf_hashes(evt_ptr->right_tensor());
}
