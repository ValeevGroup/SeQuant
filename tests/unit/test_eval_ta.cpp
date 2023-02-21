#include "catch.hpp"

#include <TiledArray/expressions/contraction_helpers.h>
#include <tiledarray.h>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/eval/eval_ta.hpp>
#include <boost/regex.hpp>

#include <cstdlib>
#include <string>
#include <vector>

auto index_label_list(std::string const& str) {
  return ranges::views::split(str, ',') | ranges::to<std::vector<std::string>>;
}

auto tensor_to_key(sequant::Tensor const& tnsr) {
  static auto const idx_rgx = boost::wregex{L"([ia])([↑↓])?(_?\\d+)"};
  auto formatter = [](boost::wsmatch mo) -> std::wstring {
    return (mo[1].str() == L"i" ? L"o" : L"v") + mo[2].str();
  };

  auto const tnsr_deparsed = sequant::deparse_expr(tnsr.clone(), false);
  return boost::regex_replace(tnsr_deparsed, idx_rgx, formatter);
}

auto tensor_to_key(std::wstring_view spec) {
  return tensor_to_key(sequant::parse_expr(spec, sequant::Symmetry::nonsymm)
                           ->as<sequant::Tensor>());
}

template <typename Tensor_t>
class rand_tensor_yield {
 private:
  TA::World& world_;
  size_t const nocc_;
  size_t const nvirt_;
  std::map<std::wstring, Tensor_t> label_to_tnsr_;

 public:
  Tensor_t make_rand_tensor(sequant::Tensor const& tnsr) const {
    using ranges::views::transform;
    using sequant::IndexSpace;

    assert(ranges::all_of(tnsr.const_braket(),
                          [](auto const& idx) {
                            return idx.space() == IndexSpace::active_occupied ||
                                   idx.space() == IndexSpace::active_unoccupied;
                          }) &&
           "Unsupported IndexSpace type found while generating tensor.");

    auto trange_vec =
        tnsr.const_braket() |
        transform([no = nocc_, nv = nvirt_](auto const& idx) {
          return TA::TiledRange1{
              0, idx.space() == IndexSpace::active_occupied ? no : nv};
        }) |
        ranges::to_vector;

    TA::TArrayD result{world_,
                       TA::TiledRange{trange_vec.begin(), trange_vec.end()}};
    result.template init_elements([](auto const&) {
      return static_cast<double>(std::rand()) / RAND_MAX;
    });
    return result;
  }

  rand_tensor_yield(TA::World& world, size_t noccupied, size_t nvirtual)
      : world_{world}, nocc_{noccupied}, nvirt_{nvirtual} {}

  Tensor_t const& operator()(sequant::Tensor const& tnsr) {
    std::wstring const label = tensor_to_key(tnsr);
    if (auto&& found = label_to_tnsr_.find(label);
        found != label_to_tnsr_.end())
      return found->second;
    auto&& success =
        label_to_tnsr_.template emplace(label, make_rand_tensor(tnsr));
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
    auto&& found = label_to_tnsr_.find(label.data());
    if (found == label_to_tnsr_.end())
      found = label_to_tnsr_.find(tensor_to_key(label));
    if (found == label_to_tnsr_.end())
      throw std::runtime_error{"attempted access of non-existent tensor!"};
    return found->second;
  }
};

auto print_node_ta = [](sequant::ExprPtr const& expr) {
  auto node = sequant::eval::to_eval_node_ta(expr);
  std::cout << node.tikz<std::string>(
                   [](auto&& n) { return "$" + n->annot() + "$"; },
                   [](auto&& n) {
                     return "label={left:" +
                            sequant::to_string(n->scalar().to_latex()) + "$}";
                   })
            << std::endl;
};

TEST_CASE("TEST_EVAL_USING_TA", "[eval]") {
  using ranges::views::transform;
  using sequant::to_eval_node;
  using sequant::eval::eval;
  using sequant::eval::eval_antisymm;
  using sequant::eval::eval_symm;
  using sequant::eval::to_eval_node_ta;
  using TA::TArrayD;
  auto parse_expr_antisymm = [](auto const& xpr) {
    return parse_expr(xpr, sequant::Symmetry::antisymm);
  };

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

  auto eval_bnode = [&yield, &manager](sequant::ExprPtr const& expr,
                                       std::string const& target_labels) {
    return eval(to_eval_node_ta(expr), index_label_list(target_labels), yield,
                manager);
  };

  auto eval_bnode_symm = [&yield, &manager](sequant::ExprPtr const& expr,
                                            std::string const& target_labels) {
    return eval_symm(to_eval_node_ta(expr), index_label_list(target_labels),
                     yield, manager);
  };

  auto eval_bnode_antisymm = [&yield, &manager](
                                 sequant::ExprPtr const& expr,
                                 std::string const& target_labels) {
    return eval_antisymm(to_eval_node_ta(expr), index_label_list(target_labels),
                         yield, manager);
  };

  SECTION("summation") {
    auto expr1 = parse_expr_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");
    auto sum1_eval = eval_bnode(expr1, "i_1,a_1");

    auto sum1_man = TArrayD{};
    sum1_man("i1,a1") =
        yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 = parse_expr_antisymm(L"2 * t_{a1}^{i1} + 1.5 * f_{i1}^{a1}");
    auto sum2_eval = eval_bnode(expr2, "i_1,a_1");

    auto sum2_man = TArrayD{};
    sum2_man("i1,a1") =
        2 * yield(L"t{a1;i1}")("a1,i1") + 1.5 * yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 =
        parse_expr_antisymm(L"1/2.0 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
    auto prod1_eval = eval_bnode(expr1, "i_4,a_1,a_4,i_1");

    TArrayD prod1_man{};
    prod1_man("i4,a1,a4,i1") = 1 / 2.0 *
                               yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                               yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_expr_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}");
    auto prod2_eval = eval_bnode(expr2, "a_1,a_2,i_1,i_2");

    auto prod2_man = TArrayD{};
    prod2_man("a1,a2,i1,i2") = -1 / 4.0 *
                               yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                               yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                               yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
  }

  SECTION("sum and product") {
    auto expr1 = parse_expr_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2} ");
    auto eval1 = eval_bnode(expr1, "a_1,a_2,i_1,i_2");

    auto man1 = TArrayD{};
    man1("a1,a2,i1,i2") = -1.0 / 4 * yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                              yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                              yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4") +
                          1.0 / 16 * yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                              yield(L"t{a1,a2;i3,i4}")("a1,a2,i3,i4") *
                              yield(L"t{a3,a4;i1,i2}")("a3,a4,i1,i2");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Antisymmetrization") {
    auto expr1 = parse_expr_antisymm(L"0.5 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_bnode_antisymm(expr1, "i_1,i_2,a_1,a_2");

    auto man1 = TArrayD{};
    man1("0,1,2,3") = yield(L"g{i1,i2;a1,a2}")("0,1,2,3") -
                      yield(L"g{i1,i2;a1,a2}")("1,0,2,3") +
                      yield(L"g{i1,i2;a1,a2}")("1,0,3,2") -
                      yield(L"g{i1,i2;a1,a2}")("0,1,3,2");

    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Symmetrization") {
    auto expr1 = parse_expr_antisymm(L"0.5 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_bnode_symm(expr1, "i_1,i_2,a_1,a_2");

    auto man1 = TArrayD{};
    man1("0,1,2,3") = yield(L"g{i1,i2;a1,a2}")("0,1,2,3") +
                      yield(L"g{i1,i2;a1,a2}")("1,0,3,2");
    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }
}