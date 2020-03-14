#include "../../src/SeQuant/domain/evaluate/eval_tensor.hpp"
#include "../../src/SeQuant/domain/evaluate/eval_tensor_builder.hpp"
#include "../sequant_setup.hpp"

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
  // CC equations
  auto cc_r = cceqvec{2, 2}(true, true, true, true);
  auto expr = cc_r[1];
  // std::wcout << expr->to_latex() << "\nexpr->to_latex()\n";
  using namespace sequant::evaluate;
  auto builder = EvalTensorBuilder{};
  auto tree = builder.build_tree(expr);
  std::wcout << "Digraph:\n";
  std::wcout << tree->to_digraph() << std::endl;

  /* auto evt_visitor = [](const EvalTensor& evt) { */
  /*     for (const auto& id: evt.get_indices()) { */
  /*         std::wcout << id << " "; */
  /*     } */
  /*     std::wcout << std::endl; */
  /* }; */
  /* tree->visit(evt_visitor); */
  return 0;
}
