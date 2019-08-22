#include "ctensor_test.hpp"
#include "contract.hpp"

int main(void) {
  // global setup...
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;

  sequant::mbpt::set_default_convention();

  using std::wcout;
  using std::endl;
  using std::uniform_real_distribution;
  using std::mt19937;
  using std::vector;
  using std::copy;
  using std::back_inserter;

  auto dist = uniform_real_distribution<double>{-1.0, 1.0};
  const size_t dim = 3;

  auto bg1 = BTensor(dim, dim, dim, dim);
  auto bg2 = BTensor(dim, dim, dim, dim);
  auto bt1 = BTensor(dim, dim);
  auto bt2 = BTensor(dim, dim, dim, dim);

  bg1.generate(bind(dist, mt19937(1)));
  bg2.generate(bind(dist, mt19937(2)));
  bt1.generate(bind(dist, mt19937(3)));
  bt2.generate(bind(dist, mt19937(4)));

  auto sg1 = STensor{L"g", {L"i_3", L"a_1"}, {L"i_1", L"i_2"}};
  auto st2 = STensor{L"t", {L"a_1", L"a_2"}, {L"i_2", L"i_3"}};

  auto t1 = CTensor{sg1};
  t1.link_btensor(bg1);
  print_ctensor(t1, "T1");

  auto t2 = CTensor{st2};
  t2.link_btensor(bt2);
  print_ctensor(t2, "T2");
  return 0;
}
