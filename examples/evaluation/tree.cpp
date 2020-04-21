#include "../sequant_setup.hpp"

#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <tiledarray.h>

#include <cstdlib>
#include <iostream>

using namespace sequant;
using namespace sequant::evaluate;

// get a sequant Tensor made out of specs
// specs -> {label, b1, ..., b(n/2), k1, ..., k(n/2)}
// eg. {"g", "i_1", "i_2", "a_1", "a_2"}
auto make_tensor_expr =
    [](const sequant::container::svector<std::string>& specs) {
      // only equal bra-ket ranks are expected
      assert((specs.size() > 2) && (specs.size() % 2 != 0));
      std::wstring label = std::wstring(specs[0].begin(), specs[0].end());
      sequant::container::svector<sequant::Index, 4> bra_indices, ket_indices;
      for (auto i = 1; i < specs.size(); ++i) {
        if (i <= specs.size() / 2)
          bra_indices.push_back(
              sequant::Index(std::wstring(specs[i].begin(), specs[i].end())));
        else
          ket_indices.push_back(
              sequant::Index(std::wstring(specs[i].begin(), specs[i].end())));
      }
      return std::make_shared<sequant::Tensor>(label, bra_indices, ket_indices);
    };

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
  using DTensorType = TA::TArrayD;
  using ContextMapType =
      sequant::container::map<HashType, std::shared_ptr<DTensorType>>;
  // initialize MADWorld
  std::wcout << "initializing MADWorld..." << std::endl;
  int argc = 1;
  char* args[]{};
  char** argv = {args};
  auto& world = TA::initialize(argc, argv);

  const size_t nocc = 10;
  const size_t nvirt = 20;

  std::wcout << "initializing random filled tensors..." << std::endl;
  TA::TiledRange tr_ov{{0, nocc}, {0, nvirt}};
  TA::TiledRange tr_oovv{{0, nocc}, {0, nocc}, {0, nvirt}, {0, nvirt}};

  auto T_ov = std::make_shared<DTensorType>(world, tr_ov);
  auto T_oovv = std::make_shared<DTensorType>(world, tr_oovv);
  auto G_oovv = std::make_shared<DTensorType>(world, tr_oovv);

  T_ov->fill_random();
  T_oovv->fill_random();
  G_oovv->fill_random();

  std::wcout << "creating sequant tensors..." << std::endl;
  auto t = make_tensor_expr({"t", "i_1", "i_2", "a_1", "a_2"});
  auto g = make_tensor_expr({"g", "i_1", "i_2", "a_1", "a_2"});
  auto A = make_tensor_expr({"A", "i_1", "i_2", "a_1", "a_2"});

  std::wcout << "building context..." << std::endl;
  ContextMapType context;
  context.insert(ContextMapType::value_type(EvalTree(t).hash_value(), T_oovv));
  context.insert(ContextMapType::value_type(EvalTree(g).hash_value(), G_oovv));

  /*for (const auto& item : context) std::wcout << item.first << std::endl; */

  std::wcout << "creating sequant expression..." << std::endl;
  auto expr = std::make_shared<Product>(Product({A, t}));

  std::wcout << "generating tree from expression..." << std::endl;
  auto tree = EvalTree(expr);

  std::wcout << "performing manual antisymmetrization..." << std::endl;
  DTensorType manual_result;
  manual_result("i,j,a,b") = (*T_oovv)("i,j,a,b") - (*T_oovv)("i,j,b,a") +
                             (*T_oovv)("j,i,b,a") - (*T_oovv)("j,i,a,b");
  auto manual_norm =
      std::sqrt(manual_result("0,1,2,3").dot(manual_result("0,1,2,3")));

  std::wcout << "calling intrusive antisymmetrization..." << std::endl;

  auto intrusive_call_result = EvalTree::_antisymmetrize(*T_oovv, 2, 2);

  std::wcout << "intrusive_call_result.trange() = "
             << intrusive_call_result.trange().rank() << std::endl;

  // auto intrusive_call_norm = std::sqrt(
  // intrusive_call_result("0,1,2,3").dot(intrusive_call_result("0,1,2,3")));

  // std::wcout << "evaluating tree..." << std::endl;
  // auto eval_result = tree.evaluate(context);
  // auto eval_norm =
  // std::sqrt(eval_result("0,1,2,3").dot(eval_result("0,1,2,3")));

  std::wcout << "manual_norm = " << manual_norm;
  // std::wcout << "\nintrusive_call_norm = " << intrusive_call_norm;
  // std::wcout << "\neval_norm = " << eval_norm;

  std::wcout << std::endl;

  TA::finalize();
  return EXIT_SUCCESS;
}
