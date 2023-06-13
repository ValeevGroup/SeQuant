#include "catch.hpp"

#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/eval.hpp>

#include <tiledarray.h>
#include <boost/regex.hpp>

#include <cstdlib>
#include <string>
#include <vector>

namespace {
auto eval_node(sequant::ExprPtr const& expr) {
  return sequant::to_eval_node<sequant::eval::EvalExprTA>(expr);
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
  std::map<std::wstring, sequant::eval::ERPtr> label_to_tnsr_;

 public:
  [[nodiscard]] Tensor_t make_rand_tensor(sequant::Tensor const& tnsr) const {
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

  sequant::eval::ERPtr operator()(sequant::Tensor const& tnsr) {
    using result_t = sequant::eval::EvalTensorTA<TA::TArrayD>;
    std::wstring const label = tensor_to_key(tnsr);
    if (auto&& found = label_to_tnsr_.find(label);
        found != label_to_tnsr_.end()) {
      //      std::cout << "label = [" << sequant::to_string(label)
      //                << "] FOUND in cache. Returning.." << std::endl;
      return found->second;
    }
    Tensor_t t = make_rand_tensor(tnsr);
    auto success = label_to_tnsr_.emplace(
        label, sequant::eval::eval_result<result_t>(std::move(t)));
    assert(success.second && "couldn't store tensor!");
    //    std::cout << "label = [" << sequant::to_string(label)
    //              << "] NotFound in cache. Creating.." << std::endl;
    return success.first->second;
  }

  template <typename T,
            typename = std::enable_if_t<sequant::eval::IsEvaluable<T>>>
  sequant::eval::ERPtr operator()(T const& node) {
    using namespace sequant::eval;
    if (node->result_type() == sequant::ResultType::Tensor) {
      assert(node->expr()->template is<sequant::Tensor>());
      return (*this)(node->as_tensor());
    }

    using result_t = EvalConstant<double>;

    assert(node->expr()->template is<sequant::Constant>());
    return eval_result<result_t>(node->as_constant().value().real());
  }

  ///
  /// \param label eg. t_vvoo, f_ov
  /// \return const ref to Tensor_t type tensor
  /// \note The tensor should be already present in the yielder cache
  ///       otherwise throws assertion error. To avoid that use the other
  ///       overload of operator() that takes sequant::Tensor const&
  sequant::eval::ERPtr operator()(std::wstring_view label) const {
    auto&& found = label_to_tnsr_.find(label.data());
    if (found == label_to_tnsr_.end())
      found = label_to_tnsr_.find(tensor_to_key(label));
    if (found == label_to_tnsr_.end())
      throw std::runtime_error{"attempted access of non-existent tensor!"};
    return found->second;
  }
};
}  // namespace

TEST_CASE("TEST_EVAL_USING_TA", "[eval]") {
  using ranges::views::transform;
  using sequant::eval::EvalExprTA;
  using sequant::eval::evaluate;
  using sequant::eval::evaluate_antisymm;
  using sequant::eval::evaluate_symm;

  using TA::TArrayD;
  auto parse_antisymm = [](auto const& xpr) {
    return parse_expr(xpr, sequant::Symmetry::antisymm);
  };

  // tnsr is assumed to be single-tiled
  auto norm = [](TArrayD const& tnsr) { return TA::norm2(tnsr); };

  std::srand(2021);

  auto& world = TA::get_default_world();

  const size_t nocc = 2, nvirt = 20;
  auto yield_ = rand_tensor_yield<TArrayD>{world, nocc, nvirt};
  auto yield = [&yield_](std::wstring_view lbl) -> TA::TArrayD const& {
    return yield_(lbl)->get<TA::TArrayD>();
  };

  auto eval = [&yield_](sequant::ExprPtr const& expr,
                        std::string const& target_labels) {
    return evaluate(eval_node(expr), target_labels, yield_)->get<TA::TArrayD>();
  };

  auto eval_symm = [&yield_](sequant::ExprPtr const& expr,
                             std::string const& target_labels) {
    return evaluate_symm(eval_node(expr), target_labels, yield_)
        ->get<TA::TArrayD>();
  };

  auto eval_antisymm = [&yield_](sequant::ExprPtr const& expr,
                                 std::string const& target_labels) {
    return evaluate_antisymm(eval_node(expr), target_labels, yield_)
        ->get<TA::TArrayD>();
  };

  SECTION("summation") {
    auto expr1 = parse_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");

    auto sum1_eval = eval(expr1, "i_1,a_1");

    auto sum1_man = TArrayD{};
    sum1_man("i1,a1") =
        yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 = parse_antisymm(L"2 * t_{a1}^{i1} + 1.5 * f_{i1}^{a1}");

    auto sum2_eval = eval(expr2, "i_1,a_1");

    auto sum2_man = TArrayD{};
    sum2_man("i1,a1") =
        2 * yield(L"t{a1;i1}")("a1,i1") + 1.5 * yield(L"f{i1;a1}")("i1,a1");

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 =
        parse_antisymm(L"1/2.0 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
    auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

    TArrayD prod1_man{};
    prod1_man("i4,a1,a4,i1") = 1 / 2.0 *
                               yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                               yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, i4}");
    auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

    auto prod2_man = TArrayD{};
    prod2_man("a1,a2,i1,i2") = -1 / 4.0 *
                               yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                               yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                               yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
  }

  SECTION("sum and product") {
    auto expr1 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");
    auto eval1 = eval(expr1, "a_1,a_2,i_1,i_2");

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
    auto expr1 = parse_antisymm(L"0.5 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");

    auto man1 = TArrayD{};
    man1("0,1,2,3") = yield(L"g{i1,i2;a1,a2}")("0,1,2,3") -
                      yield(L"g{i1,i2;a1,a2}")("1,0,2,3") +
                      yield(L"g{i1,i2;a1,a2}")("1,0,3,2") -
                      yield(L"g{i1,i2;a1,a2}")("0,1,3,2");

    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Symmetrization") {
    auto expr1 = parse_antisymm(L"0.5 * g_{i1, i2}^{a1, a2}");
    auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");

    auto man1 = TArrayD{};
    man1("0,1,2,3") = yield(L"g{i1,i2;a1,a2}")("0,1,2,3") +
                      yield(L"g{i1,i2;a1,a2}")("1,0,3,2");
    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Others") {
    using namespace std::string_literals;
    auto expr1 = parse_antisymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

    auto eval1 = evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)
                     ->get<TA::TArrayD>();

    auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
      return eval_node(x);
    }) | ranges::to_vector;

    auto eval2 =
        evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TA::TArrayD>();

    REQUIRE(norm(eval1) == Approx(norm(eval2)));
  }
}