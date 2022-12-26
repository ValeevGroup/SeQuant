#include "catch.hpp"

#include <TiledArray/expressions/contraction_helpers.h>
#include <tiledarray.h>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
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
  static auto const idx_rgx = boost::wregex{L"([ia])([⁺⁻])?(_?\\d+)"};
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

template <size_t NOCC, size_t NVIRT>
class rand_tensor_of_tensor_yield {
 public:
  using DA_tot_t = TA::DistArray<TA::Tensor<TA::Tensor<double>>>;
  using DA_t = TA::TArrayD;
  using tot_result_t = sequant::eval::tot_result_t<DA_tot_t, DA_t>;

 private:
  TA::World& world_;
  rand_tensor_yield<DA_t> rand_da;
  std::map<std::wstring, tot_result_t> label_to_tnsr_;

  DA_tot_t make_rand_tot(std::string const& annot) const {
    using ranges::views::transform;
    auto split_pos = annot.find(';');

    auto outer_lbls_str = annot.substr(0, split_pos);
    auto inner_lbls_str = annot.substr(split_pos + 1);
    auto outer_lbls = index_label_list(outer_lbls_str);
    auto inner_lbls = index_label_list(inner_lbls_str);

    auto make_trange1 = [](std::string const& str) {
      assert(str[0] == 'i' || str[0] == 'a');
      return TA::TiledRange1{0, str[0] == 'i' ? NOCC : NVIRT};
    };

    auto make_ubounds = [](std::string const& str) {
      assert(str[0] == 'i' || str[0] == 'a');
      return str[0] == 'i' ? NOCC : NVIRT;
    };

    auto outer_tr1s = outer_lbls | transform(make_trange1) | ranges::to_vector;
    auto inner_ubounds =
        inner_lbls | transform(make_ubounds) | ranges::to_vector;
    auto inner_lbounds = ranges::views::repeat(size_t{0});
    auto inner_bounds =
        ranges::views::zip(inner_lbounds, inner_ubounds) | ranges::to_vector;

    DA_tot_t result{world_,
                    TA::TiledRange{outer_tr1s.begin(), outer_tr1s.end()}};

    result.init_elements([&inner_bounds, this](auto&&) {
      auto inner_tensor = TA::Tensor<double>{TA::Range(inner_bounds)};
      ranges::generate(inner_tensor, []() {
        return static_cast<double>(std::rand()) / RAND_MAX;
      });
      return inner_tensor;
    });

    assert(result.is_initialized() && "Returning uninitialized tensor!");
    TA::get_default_world().gop.fence();
    return result;
  }

  tot_result_t make_rand_tensor(sequant::Tensor const& tnsr) const {
    using namespace sequant;

    auto annot = to_eval_node(tnsr.clone())->annot();
    auto split_pos = annot.find(';');

    return split_pos == std::string::npos
               ? tot_result_t{rand_da.make_rand_tensor(tnsr)}
               : tot_result_t{make_rand_tot(annot)};
  }

 public:
  rand_tensor_of_tensor_yield(TA::World& world)
      : world_{world}, rand_da{world, NOCC, NVIRT} {}

  tot_result_t const& operator()(sequant::Tensor const& tnsr) {
    std::wstring const label = tensor_to_key(tnsr);
    if (auto&& found = label_to_tnsr_.find(label);
        found != label_to_tnsr_.end()) {
      return found->second;
    } else {
      auto&& success = label_to_tnsr_.emplace(label, make_rand_tensor(tnsr));
      assert(success.second && "Couldn't store tensor!");
      return success.first->second;
    }
  }

  tot_result_t const& operator()(std::wstring_view label) const {
    auto&& found = label_to_tnsr_.find(label.data());
    if (found == label_to_tnsr_.end())
      found = label_to_tnsr_.find(tensor_to_key(label));
    if (found == label_to_tnsr_.end())
      throw std::runtime_error{"attempted access of non-existent tensor!"};
    return found->second;
  }
};

TEST_CASE("TEST_EVAL_USING_TA", "[eval]") {
  using ranges::views::transform;
  using sequant::to_eval_node;
  using sequant::eval::eval;
  using sequant::eval::eval_antisymm;
  using sequant::eval::eval_symm;
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
    return eval(to_eval_node(expr), index_label_list(target_labels), yield,
                manager);
  };

  auto eval_bnode_symm = [&yield, &manager](sequant::ExprPtr const& expr,
                                            std::string const& target_labels) {
    return eval_symm(to_eval_node(expr), index_label_list(target_labels), yield,
                     manager);
  };

  auto eval_bnode_antisymm = [&yield, &manager](
                                 sequant::ExprPtr const& expr,
                                 std::string const& target_labels) {
    return eval_antisymm(to_eval_node(expr), index_label_list(target_labels),
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

TEST_CASE("TEST_EVAL_TOT_USING_TA", "[eval_tot]") {
  using ranges::views::transform;
  using sequant::to_eval_node;
  using sequant::eval::eval;
  using sequant::eval::eval_antisymm;
  using sequant::eval::eval_symm;
  using TA::TArrayD;

  std::srand(2021);

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  SECTION("Tensor-of-tensor") {
    using namespace sequant;

    std::wstring_view const pno_mp2_expr =
        L"g{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}"  //
        " +  f{a_1<i_1,i_2>;a_3<i_1,i_2>}"       //
        " * t{a_3<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}"
        " - f{i_3;i_1}"  //
        " * t{a_3<i_3,i_2>,a_4<i_3,i_2>;i_3,i_2}"
        " * O{a_1<i_1,i_2>;a_3<i_3,i_2>}"
        " * O{a_2<i_1,i_2>;a_4<i_3,i_2>}";

    auto expr = parse_expr(pno_mp2_expr, Symmetry::nonsymm);

    auto& world = TA::get_default_world();
    auto tensor_yield = rand_tensor_of_tensor_yield<3, 5>{world};

    using tot_result_t = decltype(tensor_yield)::tot_result_t;
    auto cman = eval::CacheManager<tot_result_t>{{}, {}};
    using DA_tot = decltype(tensor_yield)::DA_tot_t;
    using DA = decltype(tensor_yield)::DA_t;

    auto verify_equal_tot_tile = [](auto const& lhs, auto const& rhs) {
      REQUIRE(lhs.range() == rhs.range());
      for (auto&& [t1, t2] : ranges::views::zip(lhs, rhs)) {
        for (auto&& [e1, e2] : ranges::views::zip(t1, t2))
          if (e1 != 0.0)  // Catch2 flagged (0.0 == Approx(0.0) to be false.
            REQUIRE(e1 == Approx(e2));
      }
    };

    auto const node = to_eval_node(expr);

    auto eval_result = eval::eval_tot(
        node, std::vector<std::string>{"i_1", "i_2"},
        std::vector<std::string>{"a_1", "a_2"}, tensor_yield, cman);

    DA_tot man_result{}, temp{}, lhs{}, rhs{};

    // first term evaluated manually
    temp = std::get<DA_tot>(tensor_yield(L"g{v<o,o>,v<o,o>;o,o}"));
    man_result("i1,i2;a1,a2") = temp("i1,i2;a1,a2");

    // second term evaluated manually
    lhs = std::get<DA_tot>(tensor_yield(L"f{v<o,o>;v<o,o>}"));
    rhs = std::get<DA_tot>(tensor_yield(L"t{v<o,o>,v<o,o>;o,o}"));

    TA::expressions::einsum(temp("i1,i2;a1,a2"), lhs("i1,i2;a1,a3"),
                            rhs("i1,i2;a3,a2"));
    man_result("i1,i2;a1,a2") += temp("i1,i2;a1,a2");

    // third term evaluated manually
    // using rhs from above as the T2 amplitude
    DA_tot temp2;
    TA::expressions::einsum(temp2("i2,i1,i3;a4,a3"),
                            std::get<DA>(tensor_yield(L"f{o;o}"))("i3,i1"),
                            rhs("i3,i2;a3,a4"));
    DA_tot temp3;
    TA::expressions::einsum(
        temp3("i2,i1,i3;a4,a1"), temp2("i2,i1,i3;a4,a3"),
        std::get<DA_tot>(tensor_yield(L"O{v<o,o>;v<o,o>}"))("i1,i2,i3;a1,a3"));
    TA::expressions::einsum(
        temp("i1,i2;a1,a2"), temp3("i2,i1,i3;a4,a1"),
        std::get<DA_tot>(tensor_yield(L"O{v<o,o>;v<o,o>}"))("i1,i2,i3;a2,a4"));
    man_result("i1,i2;a1,a2") -= temp("i1,i2;a1,a2");

    TA::get_default_world().gop.fence();

    verify_equal_tot_tile(eval_result.find(0).get(), man_result.find(0).get());

    auto eval_result_symm_tot = eval_symm_tot(
        node, container::vector<std::string>{"i_1", "i_2"},
        container::vector<std::string>{"a_1", "a_2"}, tensor_yield, cman);

    auto man_result_symm_tot = eval_result;
    man_result_symm_tot("0,1;2,3") -=
        eval_result("0,1;2,3");  // zero initial result
    man_result_symm_tot("0,1;2,3") += eval_result("0,1;2,3");
    man_result_symm_tot("0,1;2,3") += eval_result("1,0;3,2");

    TA::get_default_world().gop.fence();
    verify_equal_tot_tile(eval_result_symm_tot.find(0).get(),
                          man_result_symm_tot.find(0).get());

    auto eval_result_antisymm_tot = eval_antisymm_tot(
        node, container::vector<std::string>{"i_1", "i_2"},
        container::vector<std::string>{"a_1", "a_2"}, tensor_yield, cman);

    auto man_result_antisymm_tot = eval_result;
    man_result_antisymm_tot("0,1;2,3") -=
        eval_result("0,1;2,3");  // zero initial result
    man_result_antisymm_tot("0,1;2,3") += eval_result("0,1;2,3");

    man_result_antisymm_tot("0,1;2,3") -= eval_result("1,0;2,3");
    man_result_antisymm_tot("0,1;2,3") += eval_result("1,0;3,2");
    man_result_antisymm_tot("0,1;2,3") -= eval_result("0,1;3,2");

    TA::get_default_world().gop.fence();
    verify_equal_tot_tile(eval_result_antisymm_tot.find(0).get(),
                          man_result_antisymm_tot.find(0).get());
  }
}