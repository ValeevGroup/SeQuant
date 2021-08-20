#include "catch.hpp"

#include <tiledarray.h>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/eval/eval_ta.hpp>

#include <cstdlib>
#include <string>
#include <vector>

template <typename Tensor_t>
class rand_tensor_yield {
 private:
  TA::World& world_;
  size_t const nocc_;
  size_t const nvirt_;
  std::map<std::wstring, Tensor_t> label_to_tnsr_;

  Tensor_t make_rand_tensor(sequant::Tensor const& tnsr) const {
    using sequant::IndexSpace;
    using ranges::views::transform;

    assert(ranges::all_of(tnsr.const_braket(), [](auto const& idx){
      return idx.space() == IndexSpace::active_occupied
             || idx.space() == IndexSpace::active_unoccupied;
    }) && "Unsupported IndexSpace type found while generating tensor.");

    auto trange_vec = tnsr.const_braket()
                      | transform([this](auto const& idx){
                        return TA::TiledRange1{0,
                            idx.space() == IndexSpace::active_occupied
                                                     ? nocc_ : nvirt_};})
                      | ranges::to_vector;

    TA::TArrayD result{world_,
                       TA::TiledRange{trange_vec.begin(), trange_vec.end()}};
    result.template init_elements([](auto const&){
      return static_cast<double>(std::rand()) / RAND_MAX;
    });
    return result;
  }

  static std::wstring tensor_to_label(sequant::Tensor const& tnsr) {
    auto res  = std::wstring{};
         res += tnsr.label();
         res += L"_";
    for (auto const& idx: tnsr.const_braket())
      if (idx.space() == sequant::IndexSpace::active_occupied)
        res += L"o";
      else {
        assert(idx.space() == sequant::IndexSpace::active_unoccupied
               && "unsupported IndexSpace type encountered");
        res += L"v";
      }
    return res;
  }
 public:
  rand_tensor_yield(TA::World& world,
                    size_t noccupied,
                    size_t nvirtual):
                                       world_{world},
                                       nocc_{noccupied},
                                       nvirt_{nvirtual} {}

  Tensor_t const& operator()(sequant::Tensor const& tnsr) {
    std::wstring const label = tensor_to_label(tnsr);
    if (auto&& found = label_to_tnsr_.find(label);
        found != label_to_tnsr_.end())
      return found->second;
    auto&& success = label_to_tnsr_.template emplace(label,
                                                      make_rand_tensor(tnsr));
    assert(success.second && "couldn't store tensor!");
    return success.first->second;
  }

  ///
  /// \param label eg. t_vvoo, f_ov
  /// \return const ref to Tensor_t type tensor
  /// \note The tensor should be already present in the yielder cache
  ///       otherwise throws assertion error. To avoid that use the other
  ///       overload of operator() that takes sequant::Tensor const&
  Tensor_t const& operator()(std::wstring_view label) const {
    if (auto&& found = label_to_tnsr_.find(label.data());
        found != label_to_tnsr_.end())
      return found->second;
    else
      assert(false && "attempted access of non-existent tensor!");
  }
};

TEST_CASE("TEST_EVAL_USING_TA", "[eval]") {
  using ranges::views::transform;
  using TA::TArrayD;
  using sequant::parse_expr_asymm;
  using sequant::to_eval_node;
  using sequant::eval::ta::eval;
  using sequant::eval::ta::eval_antisymm;
  using sequant::eval::ta::eval_symm;

  // tnsr is assumed to be single-tiled
  auto norm = [](TArrayD const& tnsr) { return tnsr.find(0).get().norm(); };

  std::srand(2021);

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto& world = TA::get_default_world();

  const size_t nocc = 2, nvirt = 20;
  auto yield = rand_tensor_yield<TArrayD>{world, nocc, nvirt};

  // nominal(empty) cache manager
  auto manager = sequant::eval::CacheManager<TArrayD const>{{}, {}};

  auto eval_bnode = [&yield, &manager](sequant::ExprPtr const& expr) {
      return eval(to_eval_node(expr), yield, manager);
  };

  auto eval_bnode_symm = [&yield, &manager](sequant::ExprPtr const& expr) {
      return eval_symm(to_eval_node(expr), yield, manager);
  };

  auto eval_bnode_antisymm = [&yield, &manager](sequant::ExprPtr const& expr) {
      return eval_antisymm(to_eval_node(expr), yield, manager);
  };

  SECTION("summation") {
    auto expr1 = parse_expr_asymm(L"t_{a1}^{i1} + f_{i1}^{a1}");
    auto sum1_eval = eval_bnode(expr1);

    auto sum1_man = TArrayD{};
    sum1_man("0,1") = yield(L"t_vo")("0,1") + yield(L"f_ov")("1,0");

    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 = parse_expr_asymm(L"2 * t_{a1}^{i1} + 1.5 * f_{i1}^{a1}");
    auto sum2_eval = eval_bnode(expr2);

    auto sum2_man = TArrayD{};
    sum2_man("0,1") = 2 * yield(L"t_vo")("0,1") + 1.5 * yield(L"f_ov")("1,0");

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 =
        parse_expr_asymm(L"1/2.0 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
    auto prod1_eval = eval_bnode(expr1);

    TArrayD prod1_man{};
    prod1_man("i4,a1,a4,i1") =
        1 / 2.0 * yield(L"g_oovv")("i2,i4,a2,a4") * yield(L"t_vvoo")("a1,a2,i1,i2");

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_expr_asymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}");
    auto prod2_eval = eval_bnode(expr2);

    auto prod2_man = TArrayD{};
    prod2_man("a1,a2,i1,i2") = -1 / 4.0 * yield(L"g_oovv")("i3,i4,a3,a4") *
                               yield(L"t_vvoo")("a2,a4,i1,i2") * yield(L"t_vvoo")("a1,a3,i3,i4");

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
  }

  SECTION("sum and product") {
    auto expr1 = parse_expr_asymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2} ");
    auto eval1 = eval_bnode(expr1);

    auto man1 = TArrayD{};
    man1("a1,a2,i1,i2") = -1.0 / 4 * yield(L"g_oovv")("i3,i4,a3,a4") *
                              yield(L"t_vvoo")("a2,a4,i1,i2") * yield(L"t_vvoo")("a1,a3,i3,i4") +
                          1.0 / 16 * yield(L"g_oovv")("i3,i4,a3,a4") *
                              yield(L"t_vvoo")("a1,a2,i3,i4") * yield(L"t_vvoo")("a3,a4,i1,i2");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Antisymmetrization") {
    auto expr1 = parse_expr_asymm(L"0.5 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_bnode_antisymm(expr1);

    auto man1 = TArrayD{};
    man1("0,1,2,3") = yield(L"g_oovv")("0,1,2,3") - yield(L"g_oovv")("1,0,2,3") +
                      yield(L"g_oovv")("1,0,3,2") - yield(L"g_oovv")("0,1,3,2");

    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Symmetrization") {
    auto expr1 = parse_expr_asymm(L"0.5 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_bnode_symm(expr1);

    auto man1 = TArrayD{};
    man1("0,1,2,3") = yield(L"g_oovv")("0,1,2,3") + yield(L"g_oovv")("1,0,3,2");
    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Debug") {
    auto expr1 =
        parse_expr_asymm(L"1/4 * g{i2,i3;a2,a3} * t{a1,a2,a3;i1,i2,i3}");

    // auto const node1 = sequant::to_eval_node(expr1);
    // std::wcout << node1.tikz<std::wstring>(
    //                   [](auto const& n) {
    //                     return L"$" + n->tensor().to_latex() + L"$";
    //                   },
    //                   [](auto const&) { return L"circle,draw"; })
    //            << std::endl;

    auto expr2 = parse_expr_asymm(
        L"1/4 g{i3,a1;a3,a4} t{a2,a3,a4;i1,i2,i3}"
         "+ 1/4 f{i3;a3} t{a1,a2,a3;i1,i2,i3}"
         "+ 1/4 g{i3,i4;i1,a3} t{a1,a2,a3;i2,i3,i4}"
         "+ 1/4 g{i3,i4;a3,a4} t{a3;i1} t{a1,a2,a4;i2,i3,i4}"
         "+ 1/4 g{i3,i4;a3,a4} t{a3;i3} t{a1,a2,a4;i1,i2,i4}"
         "+ 1/4 g{i3,i4;a3,a4} t{a1;i3} t{a2,a3,a4;i1,i2,i4}"
        );

    auto expr3 = parse_expr_asymm(
        L"-1/4 g{i4,a1;i1,i2} t{a2,a3;i3,i4}"
         "-1/4 g{a1,a2;i1,a4} t{a3,a4;i2,i3}"
        );

    auto eval1 = eval_bnode(expr1);
    auto eval2 = eval_bnode(expr2);
    auto eval3 = eval_bnode(expr3);

    auto man1 = TArrayD{};
    auto man2 = TArrayD{};
    auto man3 = TArrayD{};

    man1("i1,a1") =
        1. / 4. * yield(L"g_oovv")("i2,i3,a2,a3") * yield(L"t_vvvooo")("a1,a2,a3,i1,i2,i3");

    man2("i1,i2,a1,a2") =   1./4 * yield(L"g_ovvv")("i3,a1,a3,a4") * yield(L"t_vvvooo")("a2,a3,a4,i1,i2,i3")
                          + 1./4 * yield(L"f_ov")("i3,a3") * yield(L"t_vvvooo")("a1,a2,a3,i1,i2,i3")
                          + 1./4 * yield(L"g_ooov")("i3,i4,i1,a3") * yield(L"t_vvvooo")("a1,a2,a3,i2,i3,i4")
                          + 1./4 * yield(L"g_oovv")("i3,i4,a3,a4") * yield(L"t_vo")("a3,i1") * yield(L"t_vvvooo")("a1,a2,a4,i2,i3,i4")
                          + 1./4 * yield(L"g_oovv")("i3,i4,a3,a4") * yield(L"t_vo")("a3,i3") * yield(L"t_vvvooo")("a1,a2,a4,i1,i2,i4")
                          + 1./4 * yield(L"g_oovv")("i3,i4,a3,a4") * yield(L"t_vo")("a1,i3") * yield(L"t_vvvooo")("a2,a3,a4,i1,i2,i4");
    man3("i1,i2,i3,a1,a2,a3") = -1./4 * yield(L"g_ovoo")("i4,a1,i1,i2") * yield(L"t_vvoo")("a2,a3,i3,i4")
                                     -1./4 * yield(L"g_vvov")("a1,a2,i1,a4") * yield(L"t_vvoo")("a3,a4,i2,i3");
    REQUIRE(norm(man1) == Approx(norm(eval1)));
    REQUIRE(norm(man2) == Approx(norm(eval2)));
    REQUIRE(norm(man3) == Approx(norm(eval3)));
  }
}
